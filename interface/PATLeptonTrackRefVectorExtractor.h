#ifndef TauAnalysis_RecoTools_PATLeptonTrackRefVectorExtractor_h
#define TauAnalysis_RecoTools_PATLeptonTrackRefVectorExtractor_h

/** \class PATLeptonTrackRefVectorExtractor
 *
 * Auxiliary class to encapsulate the different methods 
 * for accessing the tracks of pat::Electrons and pat::Muons 
 * and "signal" tracks of pat::Taus
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: PATLeptonTrackRefVectorExtractor.h,v 1.4 2010/11/20 17:07:05 veelken Exp $
 *
 */

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/Candidate/interface/Candidate.h"

template <typename T>
class PATLeptonTrackRefVectorExtractor
{
 public:
  std::vector<reco::TrackBaseRef> operator()(const T& lepton) const
  {
    assert (0);
  }
};

// add template specialization for pat::(GSF)Electrons
template <>
class PATLeptonTrackRefVectorExtractor<pat::Electron>
{
 public:
  std::vector<reco::TrackBaseRef> operator()(const pat::Electron& electron) const
  {
    std::vector<reco::TrackBaseRef> retVal;
    retVal.push_back(reco::TrackBaseRef(electron.gsfTrack()));
    return retVal;
  }
};

// add template specialization for pat::Muons
template <>
class PATLeptonTrackRefVectorExtractor<pat::Muon>
{
 public:
  std::vector<reco::TrackBaseRef> operator()(const pat::Muon& muon) const
  {
    std::vector<reco::TrackBaseRef> retVal;
    retVal.push_back(reco::TrackBaseRef(muon.track()));
    return retVal;
  }
};

// add template specialization for pat::Taus
template <>
class PATLeptonTrackRefVectorExtractor<pat::Tau>
{
 public:
  std::vector<reco::TrackBaseRef> operator()(const pat::Tau& tau) const
  {
    std::vector<reco::TrackBaseRef> retVal;
    if ( tau.isPFTau() ) {
      const reco::PFCandidateRefVector& pfChargedHadrons = tau.signalPFChargedHadrCands();
      for ( reco::PFCandidateRefVector::const_iterator pfChargedHadron = pfChargedHadrons.begin();
	    pfChargedHadron != pfChargedHadrons.end(); ++pfChargedHadron ) {
	retVal.push_back(reco::TrackBaseRef((*pfChargedHadron)->trackRef()));
      }
    } else { 
      const reco::TrackRefVector& signalTracks = tau.signalTracks();
      for ( reco::TrackRefVector::const_iterator signalTrack = signalTracks.begin();
	    signalTrack != signalTracks.end(); ++signalTrack ) {
	retVal.push_back(reco::TrackBaseRef(*signalTrack));
      }
    }
    return retVal;
  }
};

#endif
