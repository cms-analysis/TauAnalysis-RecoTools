#ifndef TauAnalysis_RecoTools_SmearedMETProducer_h
#define TauAnalysis_RecoTools_SmearedMETProducer_h

/** \class SmearedMETProducer
 *
 * Recompute pat::MET according to systematic shifts/smearing applied to collections 
 * of pat::Electrons, pat::Muons, pat::Taus and pat::Jets
 * (with the option to apply additional shift/extra smearing to pat::MET
 *  to account for particles at high |eta|/outside detector acceptance
 *  and not entering collections of PAT objects)
 * 
 * \author Christian Veelken, UC Davis;
 *         Manuel Zeise, Karlsruhe University & KIT;
 *         Michail Bachtis, University of Wisconsin
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SmearedMETProducer.h,v 1.1 2010/10/14 09:01:25 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <TFormula.h>
#include <TRandom3.h>

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

  // additional shifts/extra smearing to take into account systematic effects 
  // not accounted for by shifting/smearing PAT objects
  //
  // NOTE: the additional shifts/extra smearing is computed separately for parallel and perpendicular directions;
  //       the parallel direction is taken from the first entry in the collection of particles
  //       specified by the configuration parameter srcDirection
  //      (e.g. direction of Z boson in Z --> tau+ tau-/Z --> mu+ mu- events)
  // 
  struct extraCorrType
  {
    extraCorrType(const edm::ParameterSet& cfg)
      : srcDirection_(cfg.getParameter<edm::InputTag>("srcDirection")),
	objValExtractorParallel_(0),
	extraCorrParallel_(0),
	objValExtractorPerpendicular_(0),
	extraCorrPerpendicular_(0)	      
    {
      if ( cfg.exists("parallel") ) {
	edm::ParameterSet cfgParallel = cfg.getParameter<edm::ParameterSet>("parallel");
	initialize(cfgParallel, objValExtractorParallel_, extraCorrParallel_);     
      }

      if ( cfg.exists("perpendicular") ) {
	edm::ParameterSet cfgPerpendicular = cfg.getParameter<edm::ParameterSet>("perpendicular");
	initialize(cfgPerpendicular, objValExtractorPerpendicular_, extraCorrPerpendicular_);
      }
    }
    ~extraCorrType()
    {
      delete objValExtractorParallel_;
      delete extraCorrParallel_;
      delete objValExtractorPerpendicular_;
      delete extraCorrPerpendicular_;
    }
    void initialize(const edm::ParameterSet& cfg, ObjValExtractorBase*& objValExtractor, TFormula*& extraCorr)
    {
      edm::ParameterSet cfgObjValExtractor = cfg.getParameter<edm::ParameterSet>("parametrization");
      std::string pluginTypeObjValExtractor = cfgObjValExtractor.getParameter<std::string>("pluginType");
      objValExtractor = ObjValExtractorPluginFactory::get()->create(pluginTypeObjValExtractor, cfgObjValExtractor);

      std::string extraCorr_string = cfg.getParameter<std::string>("correction");
      extraCorr = new TFormula("perpendicular", extraCorr_string.data());
    }      
    edm::InputTag srcDirection_; 
    ObjValExtractorBase* objValExtractorParallel_;
    TFormula* extraCorrParallel_;
    ObjValExtractorBase* objValExtractorPerpendicular_;
    TFormula* extraCorrPerpendicular_;
  };
    
  extraCorrType* extraShift_;
  extraCorrType* extraSmearing_;
  
  TRandom3 rnd_;
};

#endif


