#ifndef PRESCALEPROVIDER_H
#define PRESCALEPROVIDER_H


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/split.hpp>
#include <iostream>
#include <string>

//Name: PrescaleProvider
//Simple class to provide prescale information L1 and HLT paths for a given run and lumi number
//It reads the data from three specially constructed .json files which store all the information
//about the L1 and HLT prescales and simplifies accessing this information
//The information is scraped from WBM, infact the PS tables are very clearly the PS tables from WBM
//usage:
//   PrescaleProvider psProv("path to json + first part of the name")
//   L1 prescale : psProv.l1Prescale(seedname,runnr,lumi-section)
//   HLT prescale : psProv.hltPrescale(pathname,runnr,lumi-section)
//   if a negative number is returned, it means there was an error
//
//   in both cases it automatically strips HLT version information, so anything after
//   the last "_v" is removed. Note, it expects the _v in the name. 
//   So HLT_Photon200_v3 will convert to HLT_Photon200_v and will match all versions
//   HLT_Photon200 will match nothing, it has to be HLT_Photon200_v or HLT_Photon200_vX
//
//limitations:
//   the json data is scraped from WBM so we rely on that being accurate
//   it has been noticed that sometimes for the first few lumi sections or last lumi sections 
//   the prescale column isnt in WBM, here we guess that its the same as the nearest lumi sections
//   which limited testing appears to support
//   appears to be limited to 2016 data

//technical details:
//   There are two sub classes :
//     PrescaleTbl: actually holds the L1/HLT prescale table, basically the path name + the prescales 
//                  for each column
//     RunInfo: stores the L1/HLT prescale tables for each run as well as which prescale column CMS
//              was in for a give lumisection range
//   Dont expect users to use these classes directly and instead use them via PrescaleProvider which
//   fills all the information and retrives it easily
//
//Author: Sam Harper (RAL), Dec 2017

namespace psprov{
  std::string rmTrigVersionFromName(std::string pathName){
    size_t versionIndex = pathName.rfind("_v");
    if(versionIndex!=std::string::npos) pathName.erase(versionIndex+2);
    return pathName;
  }
}

class PrescaleTbl {
public:
  struct Entry {
    std::string name;
    std::vector<unsigned int> prescales;
    bool operator<(const Entry& rhs)const{return name<rhs.name;}
    bool operator<(const std::string& rhs)const{return name<rhs;}
    friend bool operator<(const std::string& lhs,const Entry& rhs){return lhs<rhs.name;}
  };
  enum class TrigType{L1,HLT};
  
private:
  std::vector<Entry> tblData_;
  std::string menuName_;

public:
  PrescaleTbl(){}
  ~PrescaleTbl(){}
  
  void fill(const std::string& keyName,const boost::property_tree::ptree& psTree,TrigType trigType){
    for(const auto& psTbl : psTree){
      //const auto& psTbl = psTree.get_child(keyName);
      if(psTbl.first==keyName){
	fill(psTbl.second,trigType);
	break;
      }
    }
  }
  void fill(const boost::property_tree::ptree& psTbl,TrigType trigType){
    for(const auto& pathEntry : psTbl){ 
      std::vector<std::string> pathData;
      for(const auto& entry : pathEntry.second){
	const auto& val = entry.second.get_value<std::string>();
	pathData.push_back(val);
      }
      if(pathData[0]=="n") continue;

      //For HLT:
      //trigger number (meaningless) is entry 0
      //trigger name is entry 1
      //prescales are entries 2 to second to last
      //L1 seeds are last entry
      
      //For L1
      //trigger number (meaningless) is entry 0
      //trigger name is entry 1
      //prescales are entries 2 to last

      size_t lastIndxOffset = 0;
      if(trigType==TrigType::HLT) lastIndxOffset = 1;
      Entry tblEntry; 
      tblEntry.name = psprov::rmTrigVersionFromName(pathData[1]);
      for(size_t indx=2;indx+lastIndxOffset<pathData.size();indx++){
	tblEntry.prescales.push_back(std::stoi(pathData[indx]));
      }
      tblData_.emplace_back(std::move(tblEntry));
    } 
    std::sort(tblData_.begin(),tblData_.end());
  }
  
  int prescale(const std::string& path,const int psColumn)const{
    if(psColumn<0) return -1;
    const size_t psColUnsigned = static_cast<size_t>(psColumn);
    auto res = std::equal_range(tblData_.begin(),tblData_.end(),path);
    if(res.first==res.second) return -1;//path not found
    else{
      if(distance(res.first,res.second)!=1){
	std::cout <<"PrescaleTbl::prescale error "<<path<<" has multiple entries, this should not be possible"<<std::endl;
      }
      if(psColUnsigned<res.first->prescales.size()){
	return res.first->prescales[psColUnsigned];
      }else{
	std::cout <<"PrescaleTbl::prescale error "<<path<<" column "<<psColUnsigned<<" nr set "<<res.first->prescales.size()<<std::endl;
	return -2;
      }
    } 
  }
};




class RunInfo{
private:
  struct PSColData {
    int psCol;
    unsigned int minLumiSec;
    unsigned int maxLumiSec;
    
    PSColData(int iPSCol,unsigned int iMinLumiSec,unsigned int iMaxLumiSec):
      psCol(iPSCol),minLumiSec(iMinLumiSec),maxLumiSec(iMaxLumiSec){}
    bool isInRange(unsigned int lumiSec)const{
      return lumiSec>=minLumiSec && lumiSec<=maxLumiSec;
    }
    
  };

  int runnr_;
  std::string l1Menu_;
  std::string hltMenu_;
  std::string triggerMode_;
  std::string cmsswVersion_;
  std::vector<PSColData> psColData_;
 
  PrescaleTbl l1PSTbl_;
  PrescaleTbl hltPSTbl_;


public:
  explicit RunInfo(int runnr):runnr_(runnr){}
  bool operator<(const RunInfo& rhs)const{return runnr_<rhs.runnr_;}
  bool operator<(const int rhs)const{return runnr_<rhs;}
  friend bool operator<(const int lhs,const RunInfo& rhs){return lhs<rhs.runnr_;}

  int psColumn(unsigned int lumiSec)const{
    for(const auto& entry : psColData_){
      if(entry.isInRange(lumiSec)) return entry.psCol;
    }
    //okay we didnt find it,there is some WBM limitations sometimes 
    //the last and first limis are not taken into account
    //so we take our best guess
    //the guess is very good for the last lumis, probably okay for first lumis
    unsigned int minLumi=std::numeric_limits<int>::max();
    unsigned int maxLumi=0;
    for(const auto& entry : psColData_){
      minLumi = std::min(minLumi,entry.minLumiSec);
      maxLumi = std::max(maxLumi,entry.maxLumiSec);
    }
    if(lumiSec==maxLumi+1) return psColumn(maxLumi);
    if(lumiSec>0 && lumiSec<minLumi) return psColumn(minLumi);

    return -1; //didnt find it, return invalid column
  }
  int l1PrescaleFromLumiSec(const std::string& path,unsigned int lumiSec)const{
    return l1PrescaleFromPSColumn(path,psColumn(lumiSec));
  }
  int l1PrescaleFromPSColumn(const std::string& path,int psColumn)const{
    return l1PSTbl_.prescale(path,psColumn);
  }
  int hltPrescaleFromLumiSec(const std::string& path,unsigned int lumiSec)const{
    return hltPrescaleFromPSColumn(path,psColumn(lumiSec));
  }
  int hltPrescaleFromPSColumn(const std::string& path,int psColumn)const{
    return hltPSTbl_.prescale(path,psColumn);
  }
  void fill(const boost::property_tree::ptree& runInfoTree,
	    const boost::property_tree::ptree& l1PSTree,
	    const boost::property_tree::ptree& hltPSTree){
    
    const auto& runData = runInfoTree.get_child(std::to_string(runnr_));
    l1Menu_ = runData.get<std::string>("l1_menu");
    hltMenu_ = runData.get<std::string>("hlt_menu");    
    triggerMode_ = runData.get<std::string>("trig_mode");
    cmsswVersion_ = runData.get<std::string>("cmssw_version");
    
    const auto& psColData = runData.get_child("ps_cols");
    for(auto& entry : psColData){
      int colNr = std::stoi(entry.first);
      for(auto& lumiPair : entry.second){
	unsigned int minLumiSec = lumiPair.second.begin()->second.get_value<int>();
	unsigned int maxLumiSec = (++lumiPair.second.begin())->second.get_value<int>();
	//	std::cout <<"run "<<runnr_<<" colnr "<<colNr<<" min lumi "<<minLumiSec<<" max Lumi "<<maxLumiSec<<std::endl;
	psColData_.push_back(PSColData(colNr,minLumiSec,maxLumiSec));
      }
    }
    l1PSTbl_.fill(triggerMode_,l1PSTree,PrescaleTbl::TrigType::L1);
    hltPSTbl_.fill(hltMenu_,hltPSTree,PrescaleTbl::TrigType::HLT);
  }
  
  std::ostream& print(std::ostream& out)const{
    out <<runnr_<<" "<<l1Menu_<<" "<<hltMenu_<<" "<<triggerMode_;
    return out;
  }
  friend std::ostream& operator<<(std::ostream& out,const RunInfo& runInfo){
    return runInfo.print(out);
  }
    
};


class PrescaleProvider {
private:
  std::vector<RunInfo> runInfo_;
public:  
  PrescaleProvider(const std::string& jsonBaseName){
    boost::property_tree::ptree runInfoTree,l1PSTree,hltPSTree;
    boost::property_tree::read_json(jsonBaseName+"_runInfo.json",runInfoTree);
    boost::property_tree::read_json(jsonBaseName+"_l1prescales.json",l1PSTree);
    boost::property_tree::read_json(jsonBaseName+"_hltprescales.json",hltPSTree);
    for(auto& entry : runInfoTree){
      //std::cout <<"entry "<<entry.first<<" : "<<entry.second.get<std::string>("hlt_menu")<<std::endl; 
      if(entry.second.get<std::string>("hlt_menu").find("/cdaq/physics")==0){
	runInfo_.emplace_back(RunInfo(std::stoi(entry.first)));
      }
    }
    for(auto& info : runInfo_) info.fill(runInfoTree,l1PSTree,hltPSTree);
  }
  
  const RunInfo* getRunInfo(int runnr)const{
    auto res = std::equal_range(runInfo_.begin(),runInfo_.end(),runnr);
    if(res.first==res.second){
      //std::cout <<"PrescaleProvider::prescale : warning run "<<runnr<<" not found "<<std::endl;
      return nullptr; 
    }else{
      if(std::distance(res.first,res.second)>1) std::cout <<"PrescaleProvider::prescale : warning run "<<runnr<<" has multiple entries, this should not be possible "<<std::endl;
      return &(*res.first);
    }
  }

  int l1Prescale(const std::string& l1Seed,int runnr,int lumiSec)const{
    auto runInfo = getRunInfo(runnr);
    // std::cout <<"run info "<<*runInfo<<std::endl;
    if(runInfo) return runInfo->l1PrescaleFromLumiSec(psprov::rmTrigVersionFromName(l1Seed),lumiSec);
    else return -1;
  }
  int hltPrescale(const std::string& hltPath,int runnr,int lumiSec)const{
    auto runInfo = getRunInfo(runnr);
    if(runInfo) return runInfo->hltPrescaleFromLumiSec(psprov::rmTrigVersionFromName(hltPath),lumiSec);
    else return -1;
  }
    
};

#endif
