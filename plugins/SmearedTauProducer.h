// system include files
#include <memory>
#include "TauAnalysis/RecoTools/plugins/SmearedParticleProducer.h"

class SmearedTauProducer : public SmearedParticleProducer<pat::Tau,GenJetRetriever<pat::Tau> > {
   public:
  explicit SmearedTauProducer(const edm::ParameterSet& iConfig):
    SmearedParticleProducer<pat::Tau,GenJetRetriever<pat::Tau> >(iConfig),
    smearConstituents_(iConfig.getParameter<bool>("smearConstituents")),  
    hadronEnergyScale_(iConfig.getParameter<double>("hadronEnergyScale")),
    gammaEnergyScale_(iConfig.getParameter<double>("gammaEnergyScale"))
    {

    }



  ~SmearedTauProducer() {}
    
   protected:
  virtual void beginJob(const edm::EventSetup& es) {}
  virtual void endJob() {}

  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;

    std::auto_ptr<std::vector<pat::Tau> > out(new std::vector<pat::Tau> );
    Handle<std::vector<pat::Tau> > srcH;
  
    if(iEvent.getByLabel(src_,srcH) &&srcH->size()>0) 
      for(unsigned int i=0;i<srcH->size();++i) {

	pat::Tau  object = srcH->at(i);
	smear(object);

        if(smearConstituents_) {
	
	  math::XYZTLorentzVector hadronLV;
	  PFCandidateRefVector hadrons = object.signalPFChargedHadrCands();
	  if(hadrons.size()>0)
	    for(unsigned int i=0;i<hadrons.size();++i)
	    hadronLV+=hadrons.at(i)->p4();
	  //apply hadron energy scale
	  
	  math::PtEtaPhiMLorentzVector vechES(hadronLV.pt()*hadronEnergyScale_,hadronLV.eta(),hadronLV.phi(),hadronLV.M());
	  math::XYZTLorentzVector cartVechES(vechES.px(),vechES.py(),vechES.pz(),vechES.energy());
	  hadronLV=cartVechES;
	  
	  math::XYZTLorentzVector gammaLV;
	  PFCandidateRefVector gammas = object.signalPFGammaCands();
	  if(gammas.size()>0)
	    for(unsigned int i=0;i<gammas.size();++i)
	      gammaLV+=gammas.at(i)->p4();
	  
	
	  math::PtEtaPhiMLorentzVector vecgES(gammaLV.pt()*gammaEnergyScale_,gammaLV.eta(),gammaLV.phi(),gammaLV.M());
	  math::XYZTLorentzVector cartVecgES(vecgES.px(),vecgES.py(),vecgES.pz(),vecgES.energy());
	  gammaLV=cartVecgES;
	  
	  object.setP4(gammaLV+hadronLV);
	}
	
	out->push_back(object);
      }

    iEvent.put(out);

  } 

  
  // ----------member data ---------------------------
  
  bool smearConstituents_;
  double hadronEnergyScale_;
  double gammaEnergyScale_;
};



 

