// -*- C++ -*-
//
// Package:    METSmearer
// Class:      METSmearer
// 
/**\class METSmearer METSmearer.cc TauAnalysis/RecoTools/plugins/METSmearer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Manuel Zeise
//         Created:  Mo Nov 16 10:12:40 CET 2009
// $Id: METSmearer.cc,v 1.1 2010/01/12 07:19:24 zeise Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/JetReco/interface/Jet.h>
#include <DataFormats/JetReco/interface/CaloJet.h>
#include <DataFormats/METReco/interface/CaloMET.h>
#include <DataFormats/METReco/interface/CaloMETCollection.h>


#include "DataFormats/TauReco/interface/CaloTauDiscriminator.h"

#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

#include <iostream>
#include <vector>
#include <numeric>

#include <CLHEP/Random/RandGauss.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/TauReco/interface/PFTau.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

//
// class declaration
//

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
class METSmearer : public edm::EDProducer {
   public:
      explicit METSmearer(const edm::ParameterSet&);
      ~METSmearer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      template <typename T> void productNotFound(edm::Event& iEvent, edm::InputTag input);
      //bool jetIsTau(edm::Event& iEvent, edm::Handle<reco::CaloTauCollection> tausHandle, reco::Particle::LorentzVector jetp4);
      
      void repairMETp4(reco::Particle::LorentzVector * p4);
      
			reco::Particle::LorentzVector processTaus(edm::Event& iEvent, edm::Handle< std::vector<TauClass> > tausHandle, std::vector<reco::Particle::LorentzVector> * alreadySmearedJets);
			reco::Particle::LorentzVector processJets(edm::Event& iEvent, edm::Handle< std::vector<JetClass> > jetsHandle, std::vector<reco::Particle::LorentzVector> * alreadySmearedJets);
      reco::Particle::LorentzVector processCaloTowers(edm::Event& iEvent, edm::Handle< std::vector<JetClass> > jetsHandle);      

      double getSmearingFactor(double pt, std::vector<double> smearingBins, std::vector<double> smearingMeans, std::vector<double> smearingWidths);

      // ----------member data ---------------------------
      edm::InputTag jetTag_, metTag_, tauTag_, towerTag_;
      
            
      bool produceSmearedJetCollection_, produceSmearedTauCollection_, produceSmearedTowerCollection_;
      
      std::vector<double> jetSmearingBins_, jetSmearingMeans_, jetSmearingWidths_;
      std::vector<double> tauSmearingBins_, tauSmearingMeans_, tauSmearingWidths_;
      std::vector<double> towerSmearingBins_, towerSmearingMeans_, towerSmearingWidths_;
            
      std::vector<edm::InputTag> tauDiscrTags_;      
      int minPositiveTauDiscrs_;
      
      bool useAlreadySmearedJets_;
			edm::InputTag smearedJetTag_;
      bool useAlreadySmearedTaus_;
			edm::InputTag smearedTauTag_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::METSmearer(const edm::ParameterSet& iConfig)
{
	jetTag_ 	= iConfig.getParameter< edm::InputTag > ("jetTag");
	metTag_ 	= iConfig.getParameter< edm::InputTag > ("metTag");
	tauTag_ 	= iConfig.getParameter< edm::InputTag > ("tauTag");
	towerTag_ = iConfig.getParameter< edm::InputTag > ("towerTag");
	
	useAlreadySmearedJets_ = iConfig.getParameter< bool > ("useAlreadySmearedJets");
	if (useAlreadySmearedJets_)
		smearedJetTag_ = iConfig.getParameter< edm::InputTag > ("smearedJetTag");
	
	useAlreadySmearedTaus_ = iConfig.getParameter< bool > ("useAlreadySmearedTaus");
	if (useAlreadySmearedTaus_)
		smearedTauTag_ = iConfig.getParameter< edm::InputTag > ("smearedTauTag");
	
	jetSmearingBins_ = iConfig.getParameter< std::vector< double > > ("jetSmearingBins");
	jetSmearingMeans_ = iConfig.getParameter< std::vector< double > > ("jetSmearingMeans");
	jetSmearingWidths_ = iConfig.getParameter< std::vector< double > > ("jetSmearingWidths");	
	
	if (jetSmearingBins_.size() != jetSmearingMeans_.size())
		throw("jetSmearingBins and jetSmearingMeans have a different number of entries!");
	if (jetSmearingBins_.size() != jetSmearingWidths_.size())
		throw("jetSmearingBins and jetSmearingWidths have a different number of entries!");

	tauSmearingBins_ = iConfig.getParameter< std::vector< double > > ("tauSmearingBins");
	tauSmearingMeans_ = iConfig.getParameter< std::vector< double > > ("tauSmearingMeans");
	tauSmearingWidths_ = iConfig.getParameter< std::vector< double > > ("tauSmearingWidths");	
	
	if (tauSmearingBins_.size() != tauSmearingMeans_.size())
		throw("tauSmearingBins and tauSmearingMeans have a different number of entries!");
	if (tauSmearingBins_.size() != tauSmearingWidths_.size())
		throw("tauSmearingBins and tauSmearingWidths have a different number of entries!");
		

	towerSmearingBins_ = iConfig.getParameter< std::vector< double > > ("towerSmearingBins");
	towerSmearingMeans_ = iConfig.getParameter< std::vector< double > > ("towerSmearingMeans");
	towerSmearingWidths_ = iConfig.getParameter< std::vector< double > > ("towerSmearingWidths");	
	
	if (towerSmearingBins_.size() != towerSmearingMeans_.size())
		throw("towerSmearingBins and towerSmearingMeans have a different number of entries!");
	if (towerSmearingBins_.size() != towerSmearingWidths_.size())
		throw("towerSmearingBins and towerSmearingWidths have a different number of entries!");
			
	tauDiscrTags_ = iConfig.getParameter< std::vector< edm::InputTag > > ("tauDiscrTags");
	minPositiveTauDiscrs_ = iConfig.getParameter< int > ("minPositiveTauDiscrs");	
	
	produceSmearedJetCollection_ = iConfig.getParameter< bool > ("produceSmearedJetCollection");
	produceSmearedTauCollection_ = iConfig.getParameter< bool > ("produceSmearedTauCollection");
	produceSmearedTowerCollection_ = iConfig.getParameter< bool > ("produceSmearedTowerCollection");
	
	produces<std::vector<METClass> >();
	
	if (produceSmearedTowerCollection_)
	{
		produces< std::vector<TowerClass> >();
	}
	produces<std::vector< JetClass > >();
	produces<std::vector< TauClass > >("taus");
}


template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::~METSmearer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
void METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
   
	edm::Handle< std::vector<METClass>  > metHandle;
	if (!iEvent.getByLabel(metTag_, metHandle))
		return productNotFound< std::vector<METClass> >(iEvent, metTag_);
	
	const METClass met = *(metHandle->begin());
	
	edm::Handle< std::vector< JetClass > > jetsHandle;
	if (!iEvent.getByLabel(jetTag_, jetsHandle))
		return productNotFound< std::vector< reco::CaloJet > >(iEvent, jetTag_);
	LogDebug("METSmearer") << "found " << jetsHandle->size() << " jets";

  edm::Handle< std::vector< TauClass > > tausHandle; 
  if (!iEvent.getByLabel(tauTag_,tausHandle))
      return productNotFound< std::vector<TauClass> >(iEvent, tauTag_);	
	LogDebug("METSmearer") << "found " << tausHandle->size() << " taus";
 	
 	reco::Particle::LorentzVector preSmearingFromJets;
	if (useAlreadySmearedJets_)
 	{
		edm::Handle< std::vector< JetClass > > smearedJetsHandle;
		if (!iEvent.getByLabel(smearedJetTag_, smearedJetsHandle))
			return productNotFound< std::vector< JetClass > >(iEvent, smearedJetTag_);

		for (typename std::vector< JetClass >::const_iterator it = jetsHandle->begin(); it != jetsHandle->end(); ++it)
	 		preSmearingFromJets -= it->p4();
		for (typename std::vector< JetClass >::const_iterator it = smearedJetsHandle->begin(); it != smearedJetsHandle->end(); ++it)
	 		preSmearingFromJets += it->p4();
 	}
	
	reco::Particle::LorentzVector preSmearingFromTaus;
 	if (useAlreadySmearedTaus_)
 	{
		edm::Handle< std::vector< TauClass > > tausHandle;
		if (!iEvent.getByLabel(tauTag_, tausHandle))
			return productNotFound< std::vector< TauClass > >(iEvent, tauTag_);
			
		edm::Handle< std::vector< TauClass > > smearedTausHandle;
		if (!iEvent.getByLabel(smearedTauTag_, smearedTausHandle))
			return productNotFound< std::vector< TauClass > >(iEvent, smearedTauTag_);

		for (typename std::vector< TauClass >::const_iterator it = tausHandle->begin(); it != tausHandle->end(); ++it)
	 		preSmearingFromTaus -= it->p4();
		for (typename std::vector< TauClass >::const_iterator it = smearedTausHandle->begin(); it != smearedTausHandle->end(); ++it)
	 		preSmearingFromTaus += it->p4();
 	}
 		
 	std::vector<reco::Particle::LorentzVector> alreadySmearedJets;
 	reco::Particle::LorentzVector smearingFromTaus = processTaus(iEvent, tausHandle, &alreadySmearedJets);
 	reco::Particle::LorentzVector smearingFromJets = processJets(iEvent, jetsHandle, &alreadySmearedJets);
 	
 	
 	reco::Particle::LorentzVector smearingFromTowers;

	smearingFromTowers = processCaloTowers(iEvent, jetsHandle);
 	
	
	repairMETp4(&preSmearingFromTaus);
	repairMETp4(&preSmearingFromJets);
	repairMETp4(&smearingFromTaus);
	repairMETp4(&smearingFromJets);
	repairMETp4(&smearingFromTowers);

 	std::cout << "delta p4 from tau smearing:   " << smearingFromTaus << "\n";
 	std::cout << "delta p4 from jet smearing:   " << smearingFromJets << "\n";
 	std::cout << "delta p4 from tower smearing: " << smearingFromTowers << "\n";
 	std::cout << "existent tau smearing: " << preSmearingFromTaus << "\n";
 	std::cout << "existent jet smearing: " <<preSmearingFromJets << "\n";
	 	
  std::cout << "---- old MET: "<< met.p4() << "\n";  
  
 	reco::Particle::LorentzVector finalMETp4 = met.p4() - smearingFromTaus - smearingFromJets - smearingFromTowers;
 	finalMETp4.SetE(sqrt(finalMETp4.px()*finalMETp4.px()+finalMETp4.py()*finalMETp4.py()));
  METClass newMET(met.getSpecific(),finalMETp4.Et(), finalMETp4, met.vertex());
  
  std::auto_ptr< std::vector<METClass> > newmetCollection (new std::vector<METClass>);
	newmetCollection->push_back(newMET);
	iEvent.put(newmetCollection);
	  
  std::cout << "---- new MET: "<< newMET.p4() << "\n\n"; 	
 	return;

}

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
void METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::repairMETp4(reco::Particle::LorentzVector * p4)
{
 	p4->SetPz(0);
	p4->SetE(sqrt(p4->px()*p4->px()+p4->py()*p4->py()));
}

// ------------ method called once each job just before starting event loop  ------------
template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
void METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
void METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::endJob() {
}

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
reco::Particle::LorentzVector METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::processJets(edm::Event& iEvent, edm::Handle< std::vector<JetClass> > jetsHandle, std::vector<reco::Particle::LorentzVector> * alreadySmearedJets)
{
	std::auto_ptr<std::vector< JetClass > > smearedJets(new std::vector< JetClass >());
	reco::Particle::LorentzVector smearingFromJets;

	for (typename std::vector< JetClass >::const_iterator it = jetsHandle->begin(); it != jetsHandle->end(); ++it)
  { 		
  	bool ignoreJet = false;
    double tempSmearing=0.;
    
    JetClass correspondingJet(*it);
   	for (std::vector< reco::Particle::LorentzVector >::const_iterator it2 = alreadySmearedJets->begin(); it2 != alreadySmearedJets->end() && !ignoreJet; ++it2)
   		if (deltaR(*it2,it->p4())<0.1)
				ignoreJet=true;
    if (!ignoreJet && !useAlreadySmearedJets_)
			tempSmearing = getSmearingFactor(it->pt(),jetSmearingBins_, jetSmearingMeans_, jetSmearingWidths_);

    smearingFromJets+=(correspondingJet.p4())*tempSmearing;
		correspondingJet.setP4((correspondingJet.p4())*(1.+tempSmearing));
		smearedJets->push_back(correspondingJet);    	    
	}
	
	iEvent.put(smearedJets);
	
	return smearingFromJets;	
}

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
reco::Particle::LorentzVector METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::processTaus(edm::Event& iEvent, edm::Handle< std::vector<TauClass> > tausHandle, std::vector<reco::Particle::LorentzVector> * alreadySmearedJets)
{
	std::auto_ptr<std::vector< TauClass > > smearedTaus(new std::vector< TauClass >());
	reco::Particle::LorentzVector smearingFromTaus;
	
	//std::vector<reco::Particle::LorentzVector> alreadySmearedJets;
	
  for (typename std::vector<TauClass>::size_type iPFTau=0;iPFTau<tausHandle->size();iPFTau++)
  { 
    //retrieve the persistent reference to a PFTau.  This reference serves as the lookup key in the map of  discriminant values
    edm::Ref< std::vector<TauClass> > thePFTau(tausHandle,iPFTau);
			
		int positiveTauDiscrs=0;  
    for (std::vector<edm::InputTag>::iterator discr_it = tauDiscrTags_.begin(); discr_it!=tauDiscrTags_.end(); discr_it++)
    {
			edm::Handle< TauDiscriminatorClass > temp_Handle;
 
			if (!iEvent.getByLabel((*discr_it),temp_Handle))
				productNotFound< edm::View<TauDiscriminatorClass> >(iEvent, (*discr_it));

			if (((*temp_Handle)[thePFTau])!=0)
				positiveTauDiscrs += 1;
    }
    LogDebug("METSmearer") << "the tau " << thePFTau->p4() << " has " << positiveTauDiscrs << " positive tau tags ("<< minPositiveTauDiscrs_ <<" required)";

    double tempSmearing=0.;
    TauClass correspondingTau(*thePFTau);
    if (positiveTauDiscrs>=minPositiveTauDiscrs_)
    {
    	alreadySmearedJets->push_back(thePFTau->p4());
			if (!useAlreadySmearedTaus_)
				tempSmearing = getSmearingFactor(thePFTau->pt(),tauSmearingBins_, tauSmearingMeans_, tauSmearingWidths_);
	    std::cout << "the tau " << thePFTau->p4() << " has " << positiveTauDiscrs << " positive tau tags ("<< minPositiveTauDiscrs_ <<" required)\n";				
    }
    smearingFromTaus+=(correspondingTau.p4())*tempSmearing;
		correspondingTau.setP4((correspondingTau.p4())*(1.+tempSmearing));
		if (tempSmearing!=0)
		  std::cout << "---smearing factor for this tau: " << tempSmearing << "\n";		
		smearedTaus->push_back(correspondingTau);    	    
	}
	
	iEvent.put(smearedTaus, "taus");
	
	return smearingFromTaus;
}

/*
bool METSmearer::jetIsTau(edm::Event& iEvent, edm::Handle< std::vector<CaloTau> > tausHandle, reco::Particle::LorentzVector jetp4)
{      
  for (reco::CaloTauCollection::size_type iPFTau=0;iPFTau<tausHandle->size();iPFTau++)
  { 
    //retrieve the persistent reference to a PFTau.  This reference serves as the lookup key in the map of  discriminant values
    reco::CaloTauRef thePFTau(tausHandle,iPFTau);
    
		if (jetp4 != thePFTau->p4())
			continue;
			
		int positiveTauDiscrs=0;  
    for (std::vector<edm::InputTag>::iterator discr_it = tauDiscrTags_.begin(); discr_it!=tauDiscrTags_.end(); discr_it++)
    {
			edm::Handle< reco::CaloTauDiscriminator > temp_Handle;
 
			if (!iEvent.getByLabel((*discr_it),temp_Handle))
			{
				productNotFound< edm::View<reco::CaloTauDiscriminator> >(iEvent, (*discr_it));
				return false;
			}

			if (((*temp_Handle)[thePFTau])!=0)
				positiveTauDiscrs += 1;

    }
    LogDebug("METSmearer") << "the jet " << thePFTau->p4() << " has " << positiveTauDiscrs << " positive tau tags ("<< minPositiveTauDiscrs_ <<" required)";
    return (positiveTauDiscrs>=minPositiveTauDiscrs_);
	}
  return false;
}
*/

template <>
reco::Particle::LorentzVector METSmearer<reco::PFTau, reco::PFJet, CaloTower, reco::PFTauDiscriminator, reco::PFMET>::processCaloTowers(edm::Event& iEvent, edm::Handle< std::vector<reco::PFJet> > jetsHandle)
{
	return reco::Particle::LorentzVector();
}

// note the specialised class defined before
template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
reco::Particle::LorentzVector METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::processCaloTowers(edm::Event& iEvent, edm::Handle< std::vector<JetClass> > jetsHandle)
{
	reco::Particle::LorentzVector caloTowerSmearing;
	
	std::vector< CaloTowerDetId > usedCaloTowerDetIds;
	for (std::vector< reco::CaloJet >::const_iterator it = jetsHandle->begin(); it != jetsHandle->end(); ++it)
	{
		std::vector< CaloTowerDetId > caloTowerIDs = it->getTowerIndices();
		for (std::vector< CaloTowerDetId >::iterator towerIt = caloTowerIDs.begin(); towerIt != caloTowerIDs.end() ; towerIt++)
			usedCaloTowerDetIds.push_back(*towerIt);
	}
	
	std::auto_ptr< std::vector<TowerClass> > prod(new std::vector<TowerClass>());
	
	edm::Handle<CaloTowerCollection> calotowers; 
	if(!iEvent.getByLabel(towerTag_,calotowers))
	{
		productNotFound< std::vector<TowerClass> >(iEvent, towerTag_);
		return caloTowerSmearing;
	}
	
	for (typename std::vector<TowerClass>::const_iterator it = calotowers->begin(); it!=calotowers->end(); it++)
	{
		double tempSmearing = 0.;
		if (std::find(usedCaloTowerDetIds.begin (), usedCaloTowerDetIds.end (), it->id()) == usedCaloTowerDetIds.end())
			tempSmearing = getSmearingFactor(it->pt(),towerSmearingBins_, towerSmearingMeans_, towerSmearingWidths_);
		
		if (produceSmearedTowerCollection_)
		{
			tempSmearing = getSmearingFactor(it->pt(),towerSmearingBins_, towerSmearingMeans_, towerSmearingWidths_); 
			CaloTower tempCaloTower((*it));
			tempCaloTower.setP4(it->p4()*(1.+tempSmearing));
			prod->push_back(tempCaloTower);
		}
		
		caloTowerSmearing+=(it->p4())*tempSmearing;
	}
	//std::cout << caloTowerSmearing << " = caloTowerSmearing\n";
	
	if (produceSmearedTowerCollection_)
		iEvent.put(prod);
		
	return caloTowerSmearing;
}

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
double METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::getSmearingFactor(double pt, std::vector<double> smearingBins, std::vector<double> smearingMeans, std::vector<double> smearingWidths)
{	
	unsigned int rangeSize = smearingBins.size();
	double rangeMin = *(smearingBins.begin());
	double rangeMax = *(smearingBins.end()-1);
	
	double mean = 0.;
	double width = 0.;	
	if (pt<=rangeMin)
	{
		mean = smearingMeans.at(0);
		width = smearingWidths.at(0);
	}
	else
		if (pt>=rangeMax)
		{
			mean = smearingMeans.at(rangeSize-1);
			width = smearingWidths.at(rangeSize-1);			
		}
		else
		{
			for (unsigned int i = 1; i < rangeSize; i++)
			{
				if (pt>=smearingBins.at(i-1))
				{
					double lowerBin = smearingBins.at(i-1);
					double upperBin = smearingBins.at(i);
					
					double lowerValMean = smearingMeans.at(i-1);					
					double upperValMean = smearingMeans.at(i);
					mean = (upperValMean - lowerValMean)/(upperBin-lowerBin)*(pt - lowerBin) + lowerValMean;

					double lowerValWidth = smearingWidths.at(i-1);					
					double upperValWidth = smearingWidths.at(i);					
					width = (upperValWidth - lowerValWidth)/(upperBin-lowerBin)*(pt - lowerBin) + lowerValWidth;
					
					break;
				}
			}
		}
	
	double rnd = CLHEP::RandGauss::shoot(mean, width);
	LogDebug("METSmearer") << "jet with pt=" << pt << " scaled by " << rnd << " as a result of a gauss("<< mean << "," << width <<") random number";

	return rnd;
}

template <class TauClass, class JetClass, class TowerClass, class TauDiscriminatorClass, class METClass>
template <typename T>
void METSmearer<TauClass, JetClass, TowerClass, TauDiscriminatorClass, METClass>::productNotFound(edm::Event& iEvent, edm::InputTag input)
{
	std::cout << "--------------------------------------------------------\n";
	std::cout << "Could not find product with:\n"<< input << "\n";
	std::vector< edm::Handle< T >  > allHandles; 
	iEvent.getManyByType(allHandles);
	typename std::vector< edm::Handle< T > >::iterator it;
	for (it = allHandles.begin(); it != allHandles.end(); it++)
	{
		std::cout << "module label: " << (*it).provenance()->moduleLabel() << "\n";
		std::cout << "productInstanceName: " << (*it).provenance()->productInstanceName() << "\n";
		std::cout << "processName:         " << (*it).provenance()->processName() << "\n\n";
	}
	std::cout << "--------------------------------------------------------\n";
}

//define this as a plug-in

#include "FWCore/Framework/interface/MakerMacros.h"

typedef METSmearer<reco::CaloTau, reco::CaloJet, CaloTower, reco::CaloTauDiscriminator, reco::CaloMET> TauAnalysisRecoToolsCaloMETSmearer;
typedef METSmearer<PFTau, PFJet, CaloTower, reco::PFTauDiscriminator, reco::PFMET> TauAnalysisRecoToolsPFMETSmearer;
typedef METSmearer<reco::PFTau, reco::CaloJet, CaloTower, reco::PFTauDiscriminator, reco::CaloMET> TauAnalysisRecoToolsPFTauCaloMETSmearer;

DEFINE_FWK_MODULE(TauAnalysisRecoToolsCaloMETSmearer);
DEFINE_FWK_MODULE(TauAnalysisRecoToolsPFMETSmearer);
DEFINE_FWK_MODULE(TauAnalysisRecoToolsPFTauCaloMETSmearer);


