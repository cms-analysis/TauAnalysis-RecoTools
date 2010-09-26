
/** \class PATMuonVbTfSelector
 *
 * Selection of "good quality" muons according to criteria
 * defined by Vector Boson Task Force and documented in Analysis Note CMS AN-10-264
 *
 * \author Michail Bachtis,
 *         Christian Veelken
 *
 * \version $Revision: 1.21 $
 *
 * $Id: PATMuonVbTfSelector.h,v 1.21 2010/09/24 10:10:54 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include <vector>

class PATMuonVbTfSelectorImp
{
 public:
  typedef pat::MuonCollection collection;

  explicit PATMuonVbTfSelectorImp(const edm::ParameterSet&);
  ~PATMuonVbTfSelectorImp();

  std::vector<const pat::Muon*>::const_iterator begin() const { return selected_.begin(); }
  std::vector<const pat::Muon*>::const_iterator end()   const { return selected_.end();   }

  void select(const edm::Handle<collection>&, edm::Event&, const edm::EventSetup&);
    
  size_t size() const { return selected_.size(); }

 private:
  std::vector<const pat::Muon*> selected_;

//--- configuration parameters
  edm::InputTag srcBeamSpot_;

  double   maxIPxy_;            // max. transverse impact parameter of muon (global) track wrt. beam spot
  double   maxChi2red_;         // max. (normalized) chi^2 of global muon track fit per degree of freedom
  unsigned minTrackerHits_;     // min. number of hits in SiStrip + Pixel detectors
  unsigned minPixelHits_;       // min. number of hits in Pixel detector
  unsigned minMuonChamberHits_; // min. number of hits in Muon chambers
  unsigned minMuonSegments_;    // min. number of segments in Muon stations matched to inner track
};
