#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TMath.h"

using namespace pat;

const double pat::IsolatedTrack::getTrackIsolation (const reco::Track &track, const std::vector<reco::Track> &tracks, const double outerDeltaR, const double innerDeltaR) const{

  double sumPt = 0.0;

  for (const auto &t : tracks) {

    if (fabs( track.dz (t.vertex ())) > 3.0 * hypot (track.dzError (), t.dzError ()))
      continue;

    double dR = deltaR (track, t);
    if (dR < outerDeltaR && dR > innerDeltaR)
      sumPt += t.pt ();
  }

  return sumPt;
}
