#ifndef VIDCutCodes_h
#define VIDCutCodes_h

template<typename BitsDef>
class VIDCutCodes : public BitsDef  {
public:
  
  enum PassCondition{REQUIRE,IGNORE};
  
public:
  VIDCutCodes()=delete;
  ~VIDCutCodes()=delete;
  

  static bool pass(unsigned int vidBitmap,size_t bitNr,PassCondition passCond=REQUIRE){
    return pass(vidBitmap,std::vector<size_t>{bitNr},passCond);
  }
  static bool pass(unsigned int vidBitmap,const std::vector<size_t>& bitNrs,PassCondition passCond=REQUIRE){
    if(passCond==REQUIRE){
      for(auto bitNr : bitNrs){
	unsigned int mask = 0x1<<bitNr;
	if( (vidBitmap&mask)==0 ) return false;
      }
      return true;
    }else{
      unsigned int bitMask = ( ~mask(bitNrs) ) & BitsDef::kFullMask;
      return (bitMask & vidBitmap) == bitMask;
    }
  } 
  
  static unsigned int mask(size_t bitNr){return (0x1 << bitNr);}
  static unsigned int mask(const std::vector<size_t>& bitNrs){
    unsigned int val=0x0;
    for(auto bitNr : bitNrs) val|=mask(bitNr);
    return val;
  }
    
};


#endif
