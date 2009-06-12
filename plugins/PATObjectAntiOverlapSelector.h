#ifndef TauAnalysis_RecoTools_PATObjectAntiOverlapSelector_h
#define TauAnalysis_RecoTools_PATObjectAntiOverlapSelector_h

/** \class PATObjectAntiOverlapSelector
 *
 * Remove particles overlapping with other particles,
 * in order to avoud double-counting
 * 
 * \author Alfredo Gurrola, Texas A&M;
 *  modified by Konstantinos A. Petridis,
 *              Christian Veelken
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATObjectAntiOverlapSelector.h,v 1.2 2009/02/16 16:00:31 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include <vector>

template <class T>
class PATObjectAntiOverlapSelector
{
 public:
  typedef std::vector<T> collection;

  explicit PATObjectAntiOverlapSelector(const edm::ParameterSet&);

  typename std::vector<const T*>::const_iterator begin() const { return selected_.begin(); }
  typename std::vector<const T*>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<collection>&, const edm::Event&, const edm::EventSetup&);

  size_t size() const { return selected_.size(); }

 private:
  std::vector<const T*> selected_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;

  double dRmin_;
};

#endif
