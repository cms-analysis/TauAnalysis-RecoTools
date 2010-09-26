#include "TauAnalysis/RecoTools/plugins/PATMuonVbTfSelector.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/Core/interface/eventAuxFunctions.h"

#include <TMath.h>

#include <string>

//--- define default muon selection criteria defined by VBTF,
//    as documented in CMS AN-10-264
//   (to be used in case selection criteria are not specified explicitely)
double defaultMaxIPxy_              =  0.2; // cm
double defaultMaxChi2red_           = 10.;  // chi^2/nDoF
unsigned defaultMinTrackerHits_     = 10;
unsigned defaultMinPixelHits_       =  1; 
unsigned defaultMinMuonChamberHits_ =  1;
unsigned defaultMinMuonSegments_    =  2;

template<typename T>
T getConfigurationParameter(const edm::ParameterSet& cfg, const std::string& cfgParName, T* cfgParDefaultValue = 0)
{
  if ( cfg.exists(cfgParName) || !cfgParDefaultValue ) return cfg.getParameter<T>(cfgParName);
  else return (*cfgParDefaultValue);
}

PATMuonVbTfSelectorImp::PATMuonVbTfSelectorImp(const edm::ParameterSet& cfg) 
{
  srcBeamSpot_        = cfg.getParameter<edm::InputTag>("beamSpotSource");

  maxIPxy_            = getConfigurationParameter<double>  (cfg, "maxIPxy",            &defaultMaxIPxy_);
  maxChi2red_         = getConfigurationParameter<double>  (cfg, "maxChi2red",         &defaultMaxChi2red_);
  minTrackerHits_     = getConfigurationParameter<unsigned>(cfg, "minTrackerHits",     &defaultMinTrackerHits_);
  minPixelHits_       = getConfigurationParameter<unsigned>(cfg, "minPixelHits",       &defaultMinPixelHits_);
  minMuonChamberHits_ = getConfigurationParameter<unsigned>(cfg, "minMuonChamberHits", &defaultMinMuonChamberHits_);
  minMuonSegments_    = getConfigurationParameter<unsigned>(cfg, "minMuonSegments",    &defaultMinMuonSegments_);
}

PATMuonVbTfSelectorImp::~PATMuonVbTfSelectorImp() 
{
//--- nothing to be done yet...
}

void PATMuonVbTfSelectorImp::select(const edm::Handle<collection>& patMuonCollection,
				    edm::Event& evt, const edm::EventSetup& es) 
{
  selected_.clear();

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_, beamSpot);

  for ( collection::const_iterator patMuon = patMuonCollection->begin(); 
	patMuon != patMuonCollection->end(); ++patMuon ) {

    if ( !patMuon->isGlobalMuon()                                                   ) continue;
    if ( !patMuon->isTrackerMuon()                                                  ) continue;
    if ( !isValidRef(patMuon->globalTrack())                                        ) continue;
    if ( !isValidRef(patMuon->innerTrack())                                         ) continue;

    const reco::TrackRef& globalTrack = patMuon->globalTrack();
    if ( !(TMath::Abs(globalTrack->dxy(beamSpot->position())) < maxIPxy_)           ) continue;
    if ( !(globalTrack->normalizedChi2() < maxChi2red_)                             ) continue;

    const reco::HitPattern& innerTrackHitPattern = patMuon->innerTrack()->hitPattern();
    if ( !(innerTrackHitPattern.numberOfValidTrackerHits() >= (int)minTrackerHits_) ) continue;
    if ( !(innerTrackHitPattern.numberOfValidPixelHits() >= (int)minPixelHits_)     ) continue;
    
    if ( !(patMuon->numberOfChambers() >= (int)minMuonChamberHits_)                 ) continue;
    if ( !(patMuon->numberOfMatches() >= (int)minMuonSegments_)                     ) continue;

    selected_.push_back(&(*patMuon));
  }
}

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

typedef ObjectSelector<PATMuonVbTfSelectorImp> PATMuonVbTfSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATMuonVbTfSelector);
