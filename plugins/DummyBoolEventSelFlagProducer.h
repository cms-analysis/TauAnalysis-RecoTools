#ifndef TauAnalysis_RecoTools_DummyBoolEventSelFlagProducer_h
#define TauAnalysis_RecoTools_DummyBoolEventSelFlagProducer

/** \class DummyBoolEventSelFlagProducer
 *
 * Produce boolean flag that is always set to true
 * (indicating that a whole path has been run up to the end, 
 *  passing all EDFilter modules which are in this path)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: DummyBoolEventSelFlagProducer.h,v 1.4 2010/04/28 14:51:15 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DummyBoolEventSelFlagProducer : public edm::EDProducer
{
 public:

  // constructor 
  explicit DummyBoolEventSelFlagProducer(const edm::ParameterSet& cfg);

  // destructor
  virtual ~DummyBoolEventSelFlagProducer();

 private:

  // method for adding boolean flag to the event
  void produce(edm::Event&, const edm::EventSetup&);
};

#endif    
