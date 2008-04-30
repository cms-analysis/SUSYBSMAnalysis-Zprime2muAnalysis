#ifndef TOCONCRETE_H
#define TOCONCRETE_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

// Functions to cast base types of Candidates to concrete derived types
// e.g. to a reco::Muon.

template <typename T>
inline const T& toConcrete(const reco::Candidate& cand) {
  return *dynamic_cast<const T*>(&cand);
}

template <typename T>
inline const T& toConcrete(const reco::CandidateRef& cand) {
  return *dynamic_cast<const T*>(&*cand);
}

template <typename T>
inline const T& toConcrete(const reco::CandidateBaseRef& cand) {
  return *dynamic_cast<const T*>(&*cand);
}

template <typename T>
inline const T* toConcretePtr(const reco::Candidate& cand) {
  return dynamic_cast<const T*>(&cand);
}

template <typename T>
inline const T* toConcretePtr(const reco::CandidateRef& cand) {
  return dynamic_cast<const T*>(&*cand);
}

template <typename T>
inline const T* toConcretePtr(const reco::CandidateBaseRef& cand) {
  return dynamic_cast<const T*>(&*cand);
}

#endif
