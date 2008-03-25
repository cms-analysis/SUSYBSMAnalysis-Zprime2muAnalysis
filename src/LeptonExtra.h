#ifndef LEPTONEXTRA_H
#define LEPTONEXTRA_H

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/Utilities/interface/Exception.h"

// details about the number of reconstruction levels stored
const int NUM_REC_LEVELS = 4;
const int MAX_LEVELS = 8;
enum RECLEVEL { lgen, l1, l2, l3, lgmr, ltk, lfms, lpmr, lbest };
const std::string str_level[MAX_LEVELS+1] = {
  "Gen", " L1", " L2", " L3", "GMR", "Tracker-only", "TPFMS", "PMR", "OPT"
};
const std::string str_level_short[MAX_LEVELS+1] = {
  "GN", "L1", "L2", "L3", "GR", "TK", "FS", "PR", "BS"
};

inline void checkRecLevel(const int level, const char* name) {
  if (level < 0 || level > MAX_LEVELS)
    throw cms::Exception(name)
      << "invalid level " << level << " in " << name << "!\n";
}

#include <iostream>

namespace zp2mu {
  class LeptonExtra {
  public:
    LeptonExtra() {
      init_();
    }
      
    LeptonExtra(int iden, int rec) {
      init_(iden, rec);
      valid_ = 1;
    }
    
    LeptonExtra(const LeptonExtra& rhs) {
      copy_(rhs);
    }

    LeptonExtra& operator=(const LeptonExtra& rhs) {
      copy_(rhs);
      return *this;
    }

    void setValid(const int valid) {
      valid_ = valid;
    }

    int valid() const {
      return valid_;
    }

    void setId(const int id) {
      id_ = id;
    }

    int id() const {
      return id_;
    }
    
    void setRecLevel(const int recLevel) {
      recLevel_ = recLevel;
    }

    int recLevel() const {
      return recLevel_;
    }
    
    void setSeedIndex(const int index) {
      seedIndex_ = index;
    }

    int seedIndex() const {
      return seedIndex_;
    }
    
    void setSeedMatch(const int index, const int level) {
      checkRecLevel(level, "setSeedMatch");
      seedMatch_[level] = index;
    }

    int seedMatch(const int level) const {
      checkRecLevel(level, "seedMatch");
      return seedMatch_[level];
    }

    void setClosestMatch(const int index, const int level) {
      checkRecLevel(level, "setClosestMatch");
      closestMatch_[level] = index;
    }

    int closestMatch(const int level) const {
      checkRecLevel(level, "setClosestMatch");
      return closestMatch_[level];
    }

    void setClosestPhoton(const reco::Particle::LorentzVector& ph) {
      photon_ = ph;
    }

    const reco::Particle::LorentzVector& closestPhoton() const {
      return photon_;
    }

  private:
    void init_(int iden=-1, int rec=-1) {
      valid_ = 0;
      id_ = iden;
      recLevel_ = rec;
      seedIndex_ = -999;
      for (int irec = 0; irec < MAX_LEVELS; irec++) {
	setClosestMatch(irec == recLevel_ ? id_ : -999, irec);
	setSeedMatch(irec >= l3 && irec == recLevel_ ? id_ : -999, irec);
      }
    }

    void copy_(const LeptonExtra& rhs) {
      valid_ = rhs.valid();
      id_ = rhs.id();
      recLevel_ = rhs.recLevel();
      seedIndex_ = rhs.seedIndex();
      for (int i = 0; i < MAX_LEVELS; i++) {
	setClosestMatch(rhs.closestMatch(i), i);
	setSeedMatch(rhs.seedMatch(i), i);
      }
      photon_ = rhs.closestPhoton();
    }
  
    int valid_;
    int id_;
    int recLevel_;
    int seedIndex_;
    int seedMatch_[MAX_LEVELS];
    int closestMatch_[MAX_LEVELS];
    reco::Particle::LorentzVector photon_;
  };
}

#endif
