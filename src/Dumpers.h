#ifndef Zp2mu_Dumps_h
#define Zp2mu_Dumps_h

#include <iosfwd>

namespace edm {
  class TriggerNames;
}

namespace reco {
  class GenParticle;
  class HitPattern;
  class Track;
}

namespace pat {
  class CompositeCandidate;
  class Electron;
  class Muon;
}

int mlprintf(const char* category, const char* fmt, ...);

std::ostream& operator<<(std::ostream& out, const reco::GenParticle& gen);
std::ostream& operator<<(std::ostream& out, const reco::HitPattern& hp);
std::ostream& operator<<(std::ostream& out, const reco::Track& tk);
std::ostream& operator<<(std::ostream& out, const pat::Muon& mu);
std::ostream& operator<<(std::ostream& out, const pat::Electron& el);
std::ostream& operator<<(std::ostream& out, const reco::CandidateBaseRef& cbr);
std::ostream& operator<<(std::ostream& out, const pat::CompositeCandidate& dil);
std::ostream& operator<<(std::ostream& out, const edm::TriggerNames& t);

#endif
