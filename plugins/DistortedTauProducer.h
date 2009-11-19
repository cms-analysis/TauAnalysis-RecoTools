// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include <TFile.h>
#include <TH1F.h>


class DistortedTauProducer : public edm::EDProducer {
   public:
    explicit DistortedTauProducer(const edm::ParameterSet& iConfig);
    ~DistortedTauProducer();
    
   private:
    virtual void beginJob(const edm::EventSetup& es);
    virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup); 
    virtual void endJob();

    TH1F* integratedHisto(TH1F*);
    double sampleFromHisto(TH1F *, double);
      
      // ----------member data ---------------------------
      
      edm::InputTag src_;        //input Collection
      edm::FileInPath fileName_; //Name of the file with the histograms

      //Smear MC particle!
      bool smearMCParticle_;

      bool smearConstituents_;

      //Histogram Smearing
      bool smearFromPtHistogram_;
      bool smearFromEtaHistogram_;
      std::string ptHisto_;        //Histogram name in File
      std::string etaHisto_;        //Histogram name in File

      //Flat smearing
      double energyScale_;        //energy Scale 

      double hadronEnergyScale_;        //hadronEnergy Scale 
      double gammaEnergyScale_;        //Gamma Energy Scale 

      double deltaEta_;           //delta Eta
      double deltaPt_;            //delta Pt


      //Gaussian Smearing
      double gSigmaPt_;
      double gSigmaEta_;

      //Private data
      TFile *f;

      TH1F *ptHisto;
      TH1F *etaHisto;

      edm::Service<edm::RandomNumberGenerator> rng;
};


