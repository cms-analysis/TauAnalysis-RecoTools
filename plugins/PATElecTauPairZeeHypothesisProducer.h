#ifndef TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisProducer_h
#define TauAnalysis_RecoTools_PATElecTauPairZeeHypothesisProducer_h

/** \class PATElecTauPairZeeHypothesisProducer
 *
 * Produce data-formats providing information 
 * about compatibility of an elecron + tau-jet pair
 * with the hypothesis of being a pair of electron,
 * resulting from a Z --> e+ e- decay
 * 
 * \authors Christian Veelken
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATElecTauPairZeeHypothesisProducer.h,v 1.1 2009/07/17 14:21:33 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "TauAnalysis/RecoTools/interface/PATElecTauPairZeeHypothesisAlgorithm.h"

class PATElecTauPairZeeHypothesisProducer : public edm::EDProducer 
{
 public:

  explicit PATElecTauPairZeeHypothesisProducer(const edm::ParameterSet&);
  ~PATElecTauPairZeeHypothesisProducer();

  void produce(edm::Event&, const edm::EventSetup&);

 private:

  PATElecTauPairZeeHypothesisAlgorithm algorithm_;
  
  edm::InputTag srcElecTauPairs_;
};

#endif

