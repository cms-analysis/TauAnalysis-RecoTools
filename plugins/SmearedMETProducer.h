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
 * \version $Revision: 1.2 $
 *
 * $Id: SmearedMETProducer.h,v 1.2 2010/11/08 09:23:06 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

class SmearedMETProducer : public edm::EDProducer 
{
 public:
  explicit SmearedMETProducer(const edm::ParameterSet&);
  ~SmearedMETProducer();
    
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&); 
      
//--- configuration parameters

  // "original" MET collection
  edm::InputTag src_;

  struct smearedParticleType
  {
    smearedParticleType(const edm::ParameterSet& cfg)
      : srcOriginal_(cfg.getParameter<edm::InputTag>("srcOriginal")),
	srcSmeared_(cfg.getParameter<edm::InputTag>("srcSmeared"))
    {}
    ~smearedParticleType() {}
    edm::InputTag srcOriginal_;
    edm::InputTag srcSmeared_;
  };

  // collections of smeared pat::Electrons, pat::Muons, pat::Taus and pat::Jets
  // the systematic shifts/smearing of which is to enter the MET recomputation
  std::vector<smearedParticleType> smearedParticleCollections_;
};

#endif


