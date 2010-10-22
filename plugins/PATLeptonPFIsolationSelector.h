#ifndef TauAnalysis_RecoTools_PATLeptonPFIsolationSelector_h
#define TauAnalysis_RecoTools_PATLeptonPFIsolationSelector_h

/** \class PATLeptonPFIsolationSelector
 *
 * Flexible selection of pat::Leptons based on charged hadrons, 
 * neutral hadrons and photons/pi0s reconstructed by particle-flow algorithm;
 * a pat::Electron, pat::Muon or pat::Tau passed the selection
 * if the sum of deposits of the speficied type and with Pt > ptMin 
 * in an annulus within dRvetoCone (inner cone) and dRisoCone (outer cone)
 * is more than sumPtMin and less than sumPtMax 
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATLeptonPFIsolationSelector.h,v 1.1 2010/10/14 09:01:25 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

#include <vector>
#include <string>

template <class T>
class PATLeptonPFIsolationSelector
{
 public:
  typedef std::vector<T> collection;

  explicit PATLeptonPFIsolationSelector(const edm::ParameterSet&);
  ~PATLeptonPFIsolationSelector();

  typename std::vector<const T*>::const_iterator begin() const { return selected_.begin(); }
  typename std::vector<const T*>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<collection>&, const edm::Event&, const edm::EventSetup&);

  size_t size() const { return selected_.size(); }

 private:
  std::vector<const T*> selected_;

  struct pfIsoConfigType
  {
    pfIsoConfigType(reco::PFCandidate::ParticleType pfParticleType, const edm::ParameterSet& cfg)
      : pfParticleType_(pfParticleType),
	ptMin_(cfg.getParameter<double>("ptMin")),
	dRvetoCone_(cfg.getParameter<double>("dRvetoCone")),
	dRisoCone_(cfg.getParameter<double>("dRisoCone"))
    {
      dEtaVeto_ = cfg.exists("dEtaVeto") ? cfg.getParameter<double>("dEtaVeto") : -1.;
      dPhiVeto_ = cfg.exists("dPhiVeto") ? cfg.getParameter<double>("dPhiVeto") : -1.;
    }
    ~pfIsoConfigType() {}

    double compSumPt(const reco::PFCandidateCollection& pfCandidates, const reco::Particle::LorentzVector& patLeptonP4)
    {
      double sumPt = 0.;

      for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates.begin();
	    pfCandidate != pfCandidates.end(); ++pfCandidate ) {
	if ( pfCandidate->particleId() != pfParticleType_ ) continue;

	if ( pfCandidate->pt() <  ptMin_ ) continue;

	double dR = deltaR(pfCandidate->p4(), patLeptonP4);
	if ( dR < dRvetoCone_ || dR > dRisoCone_ ) continue;

	if ( TMath::Abs(pfCandidate->eta() - patLeptonP4.eta()) < dEtaVeto_ ) continue;
	if ( TMath::Abs(pfCandidate->phi() - patLeptonP4.phi()) < dPhiVeto_ ) continue;

	sumPt += pfCandidate->pt();
      }    

      return sumPt;
    }

    reco::PFCandidate::ParticleType pfParticleType_;

    double ptMin_;

    double dEtaVeto_;
    double dPhiVeto_;

    double dRvetoCone_;
    double dRisoCone_;
  };

  edm::InputTag pfCandidateSrc_;

  pfIsoConfigType* pfChargedHadronIso_;
  pfIsoConfigType* pfNeutralHadronIso_;
  pfIsoConfigType* pfPhotonIso_;

  double sumPtMin_;
  double sumPtMax_;
  
  enum { kAbsoluteIso, kRelativeIso };
  int sumPtMethod_;

  int cfgError_;
};

#endif
