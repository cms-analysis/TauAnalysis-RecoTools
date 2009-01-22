#ifndef TauAnalysis_RecoTools_HistManagerBase_h
#define TauAnalysis_RecoTools_HistManagerBase_h

/** \class HistManagerBase
 *
 * Base-class for histogram booking and filling in physics analyses
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.25 $
 *
 * $Id: Cut.h,v 1.25 2008/02/04 10:44:27 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HistManagerBase
{
 public:
  // constructor 
  explicit HistManagerBase() {}
  
  // destructor
  virtual ~HistManagerBase() {}

  // methods for booking and filling of histograms
  virtual void bookHistograms(const edm::EventSetup&) = 0;
  virtual void fillHistograms(const edm::Event&, const edm::EventSetup&) = 0;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<HistManagerBase* (const edm::ParameterSet&)> HistManagerPluginFactory;

#endif  

