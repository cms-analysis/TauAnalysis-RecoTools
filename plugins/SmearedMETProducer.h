#ifndef TauAnalysis_RecoTools_SmearedMETProducer_h
#define TauAnalysis_RecoTools_SmearedMETProducer_h

/** \class SmearedMETProducer
 *
 * Recompute pat::MET according to systematic shifts/smearing applied to collections 
 * of pat::Electrons, pat::Muons, pat::Taus and pat::Jets
 * 
 * \author Christian Veelken, UC Davis;
 *         Manuel Zeise, Karlsruhe University & KIT;
 *         Michail Bachtis, University of Wisconsin
 *
 * \version $Revision: 1.3 $
 *
 * $Id: SmearedMETProducer.h,v 1.3 2010/12/12 08:52:26 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <TRandom3.h>

class SmearedMETProducer : public edm::EDProducer 
{
 public:
  explicit SmearedMETProducer(const edm::ParameterSet&);
  ~SmearedMETProducer();
    
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&); 

  std::string moduleLabel_;
      
//--- configuration parameters

  // "original" MET collection
  edm::InputTag src_;

  struct smearedParticleType
  {
    smearedParticleType(const edm::ParameterSet& cfg)
      : srcOriginal_(cfg.getParameter<edm::InputTag>("srcOriginal")),
	srcSmeared_(cfg.getParameter<edm::InputTag>("srcSmeared"))
    {
      smearByResolutionUncertainty_ = cfg.exists("smearByResolutionUncertainty") ? 
	cfg.getParameter<double>("smearByResolutionUncertainty") : 0.;
    }
    ~smearedParticleType() {}
    edm::InputTag srcOriginal_;
    edm::InputTag srcSmeared_;
    double smearByResolutionUncertainty_;
  };

  // collections of smeared pat::Electrons, pat::Muons, pat::Taus and pat::Jets
  // the systematic shifts/smearing of which is to enter the MET recomputation
  std::vector<smearedParticleType> smearedParticleCollections_;

  TRandom rnd_;
};

#endif


