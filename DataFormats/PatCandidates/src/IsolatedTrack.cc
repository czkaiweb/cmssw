#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

using namespace pat;

const double pat::IsolatedTrack::getTrackIsolation (const reco::Track &track, const vector<reco::Track> &tracks, const double outerDeltaR, const double innerDeltaR= 1.0e-12) {

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

const double pat::IsolatedTrack::caloTotDR5NoPU ( RhoType rhoType = All, CaloType caloType = Sum) {
// For reference, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_AN1
  double rho;
  switch (rhoType) {
    case All:
      rho = rhoPUCorr();
      break;
    case Calo:
      rho = rhoPUCorrCalo();
      break;
    case CentralCalo:
      rho = rhoPUCorrCentralCalo();
      break;
    default:
      throw cms::Exception("FatalError") << "Unknown or not implemented rho type requested, type:" << rhoType;
  }

  double rawCaloTot = 0.0;
  double dR = 0.5;
  int intDR = dR * 10.0;
  switch (caloType) {
    case Sum:
      rawCaloTot = caloDR5();
    case EM:
      rawCaloTot = caloEMDR5();
    case Had:
      rawCaloTot = caloHadDR5();
    default:
      throw cms::Exception("FatalError")<< "Unknown or not implemented calo type requested, type:" << caloType;
    }
  double caloCorr = rho * TMath::Pi() * dR * dR;  // Define effective area as pi*r^2, where r is radius of DeltaR cone.
  double caloTotNoPU = TMath::Max(0., rawCaloTot - caloCorr);
  return caloTotNoPU;
}
