

#ifndef TauAnalysis_RecoTools_PATElectronConversionFinder_h
#define TauAnalysis_RecoTools_PATElectronConversionFinder_h

// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

//Tracker tracks
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "Math/GenVector/VectorUtil.h"


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

class PATElectronConversionFinderImp
{
   public:
  
      typedef std::vector<pat::Electron> PATElectronCollection;
      typedef std::vector<pat::Electron> collection;

       PATElectronConversionFinderImp(const edm::ParameterSet& iConfig);
      ~PATElectronConversionFinderImp();
  
      std::vector<const pat::Electron*>::const_iterator begin() const { return selected_.begin(); }
      std::vector<const pat::Electron*>::const_iterator end() const { return selected_.end();}

      void select(const edm::Handle<PATElectronCollection>&,edm::Event & iEvent, const edm::EventSetup & iSetup);

   private:
      double GetDcaBetweenTracks(const edm::Event&, 
				 const edm::EventSetup&,
				 const MagneticField*,
				 const reco::Track&,
				 const reco::GsfTrack&);  
      
      
      std::vector<const pat::Electron*> selected_;

      edm::InputTag elecs_,tracks_;
      double cotTheta_,dca_,dRin_,eventWeight_;
      bool doPXBcut_;bool doHistos_;
      double nMax_;
      TH1F *hnTracksPass, *hPXBLayer, *hPXFDisk;
      TH2F *hDcotTheta, *hDca; 
};

#endif