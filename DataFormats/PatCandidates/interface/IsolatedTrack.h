#ifndef __DataFormats_PatCandidates_IsolatedTrack_h__
#define __DataFormats_PatCandidates_IsolatedTrack_h__

/*
  \class    pat::IsolatedTrack IsolatedTrack.h "DataFormats/PatCandidates/interface/IsolatedTrack.h"
  \brief Small class to store key info on isolated tracks
   pat::IsolatedTrack stores important info on isolated tracks. Draws from
   packedPFCandidates, lostTracks, and generalTracks.
  \author   Bennett Marsh
*/

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"

namespace pat {

    class IsolatedTrack : public reco::LeafCandidate {

      public:

        IsolatedTrack() :
          LeafCandidate(0, LorentzVector(0,0,0,0)),
          pfIsolationDR03_(pat::PFIsolation()),
          miniIsolation_(pat::PFIsolation()),
          matchedCaloJetEmEnergy_(0.), matchedCaloJetHadEnergy_(0.),
          pfLepOverlap_(false), pfNeutralSum_(0.),
          dz_(0.), dxy_(0.), dzError_(0.), dxyError_(0.), fromPV_(-1), trackQuality_(0),
          dEdxStrip_(0), dEdxPixel_(0), hitPattern_(reco::HitPattern()),
          crossedEcalStatus_(std::vector<uint16_t>()),
          crossedHcalStatus_(std::vector<uint32_t>()),
          deltaEta_(0), deltaPhi_(0),
          packedCandRef_(PackedCandidateRef()),
          nearestPFPackedCandRef_(PackedCandidateRef()),
          nearestLostTrackPackedCandRef_(PackedCandidateRef()),
          assocEMCaloDR05_           (-1.),
          assocHadCaloDR05_          (-1.),
          rhoPUCorr_               (-1.),
          rhoPUCorrCalo_           (-1.),
          rhoPUCorrCentralCalo_    (-1.),
          trackIsoDR05_            (-1.){}

        explicit IsolatedTrack(const PFIsolation &iso, const PFIsolation &miniiso, float caloJetEm, float caloJetHad,
                               bool pfLepOverlap, float pfNeutralSum,
                               const LorentzVector &p4, int charge, int id,
                               float dz, float dxy, float dzError, float dxyError,
                               const reco::HitPattern &hp, float dEdxS, float dEdxP, int fromPV, int tkQual,
                               const std::vector<uint16_t> &ecalst,
                               const std::vector<uint32_t> &hcalst, int dEta, int dPhi,
                               const PackedCandidateRef &pcref, const PackedCandidateRef &refToNearestPF, const PackedCandidateRef &refToNearestLostTrack,
                               const reco::Track &track,
                               const std::vector<reco::Track> &tracks):
          LeafCandidate(charge, p4, Point(0.,0.,0.), id),
          pfIsolationDR03_(iso), miniIsolation_(miniiso),
          matchedCaloJetEmEnergy_(caloJetEm), matchedCaloJetHadEnergy_(caloJetHad),
          pfLepOverlap_(pfLepOverlap), pfNeutralSum_(pfNeutralSum),
          dz_(dz), dxy_(dxy), dzError_(dzError), dxyError_(dxyError),
          fromPV_(fromPV), trackQuality_(tkQual), dEdxStrip_(dEdxS), dEdxPixel_(dEdxP),
          hitPattern_(hp),
          crossedEcalStatus_(ecalst), crossedHcalStatus_(hcalst),
          deltaEta_(dEta), deltaPhi_(dPhi),
          packedCandRef_(pcref),
          nearestPFPackedCandRef_(refToNearestPF),
          nearestLostTrackPackedCandRef_(refToNearestLostTrack),
          assocEMCaloDR05_           (-1.),
          assocHadCaloDR05_          (-1.),
          rhoPUCorr_               (-1.),
          rhoPUCorrCalo_           (-1.),
          rhoPUCorrCentralCalo_    (-1.),
          trackIsoDR05_            (getTrackIsolation (track, tracks, 0.5)){}

        ~IsolatedTrack() override {}

        const PFIsolation& pfIsolationDR03() const { return pfIsolationDR03_; }

        const PFIsolation& miniPFIsolation() const { return miniIsolation_; }

        float matchedCaloJetEmEnergy() const { return matchedCaloJetEmEnergy_; }
        float matchedCaloJetHadEnergy() const { return matchedCaloJetHadEnergy_; }

        bool pfLepOverlap() const { return pfLepOverlap_; }
        float pfNeutralSum() const { return pfNeutralSum_; }

        float dz() const { return dz_; }
        float dzError() const override { return dzError_; }
        float dxy() const { return dxy_; }
        float dxyError() const override { return dxyError_; }

        int fromPV() const { return fromPV_; }

        bool isHighPurityTrack() const 
          {  return (trackQuality_ & (1 << reco::TrackBase::highPurity)) >> reco::TrackBase::highPurity; }
        bool isTightTrack() const 
          {  return (trackQuality_ & (1 << reco::TrackBase::tight)) >> reco::TrackBase::tight; }
        bool isLooseTrack() const 
          {  return (trackQuality_ & (1 << reco::TrackBase::loose)) >> reco::TrackBase::loose; }

       const reco::HitPattern& hitPattern() const { return hitPattern_; }


       //! just the status code part of an EcalChannelStatusCode for all crossed Ecal cells
       const std::vector<uint16_t>& crossedEcalStatus() const { return crossedEcalStatus_; }
       //! just the status code part of an HcalChannelStatus for all crossed Hcal cells
       const std::vector<uint32_t>& crossedHcalStatus() const { return crossedHcalStatus_; }

       //! difference in eta/phi between initial traj and intersection w/ ecal
       //! Values are between +-0.5 with a precision of 0.002
       float deltaEta() const { return float(deltaEta_) / 500.f; }
       float deltaPhi() const { return float(deltaPhi_) / 500.f; }

       const PackedCandidateRef& packedCandRef() const { return packedCandRef_; }
       const PackedCandidateRef& nearestPFPackedCandRef() const { return nearestPFPackedCandRef_; }
       const PackedCandidateRef& nearestLostTrackPackedCandRef() const { return nearestPFPackedCandRef_; }

       enum RhoType { All, Calo, CentralCalo };
       enum CaloType { Sum, EM, Had };

       //////////////////////////////////////
       // Un-corrected (rho) calo energies
       //////////////////////////////////////

       // New calculation that uses all rec hits in dR < 0.5 cone.
       const float assocEMCaloDR05 ()  const { return this->assocEMCaloDR05_; };
       const float assocHadCaloDR05 () const { return this->assocHadCaloDR05_; };
       const float assocCaloDR05 ()    const { return this->assocEMCaloDR05_ + this->assocHadCaloDR05_; };

       //////////////////////////////////////
       // Rho-corrected calo energies
       //////////////////////////////////////

       // New calculation that uses all rec hits in dR < 0.5 cone.
       const float assocAllCaloDR05NoPU ()           const;
       const float assocAllEMCaloDR05NoPU ()         const;
       const float assocAllHadCaloDR05NoPU ()        const;
     
       const float assocCaloDR05NoPUCalo ()          const;
       const float assocCaloDR05NoPUCaloEm ()        const;
       const float assocCaloDR05NoPUCaloHad ()       const;

       const float caloDR05NoPUCentralCalo ()        const;
       const float caloDR05NoPUCentralCaloEm ()      const;
       const float caloDR05NoPUCentralCaloHad ()     const;

       //////////////////////////////////////
       // Set calo energies
       //////////////////////////////////////

       void set_assocEMCaloDR05 (double value) { assocEMCaloDR05_  = value; };
       void set_assocHadCaloDR05 (double value) { assocHadCaloDR05_ = value; };

       //////////////////////////////////////
       // Set rhos
       //////////////////////////////////////

       void set_rhoPUCorr  (double value) { rhoPUCorr_   = value; };
       void set_rhoPUCorrCalo         (double value) { rhoPUCorrCalo_   = value; };
       void set_rhoPUCorrCentralCalo  (double value) { rhoPUCorrCentralCalo_   = value; };

       //////////////////////////////////////
       // Set track isolations
       //////////////////////////////////////

       void set_trackIsoDR05 (double value) { trackIsoDR05_ = value; };

       void set_trackIsoNoPUDR05 (double value) { trackIsoNoPUDR05_ = value; };

       void set_trackIsoNoFakesDR05 (double value) { trackIsoNoFakesDR05_ = value; };

       void set_trackIsoNoPUNoFakesDR05 (double value) { trackIsoNoPUNoFakesDR05_ = value; };

       void set_dEdx_pixel (double value, 
                            double error, 
                            int nSaturatedMeasurements, 
                            unsigned int nMeasurements) { 
           dEdxPixel_ = value; 
           dEdxPixel_Error_ = error;
           dEdxPixel_numberOfSaturatedMeasurements_ = nSaturatedMeasurements;
           dEdxPixel_numberOfMeasurements_ = nMeasurements;
       };
       void set_dEdx_strip (double value, 
                           double error, 
                           int nSaturatedMeasurements, 
                           unsigned int nMeasurements) { 
           dEdxStrip_ = value; 
           dEdxStrip_Error_ = error;
           dEdxStrip_numberOfSaturatedMeasurements_ = nSaturatedMeasurements;
           dEdxStrip_numberOfMeasurements_ = nMeasurements;
       };

       //////////////////////////////////////
       // Get rhos
       //////////////////////////////////////

       const float rhoPUCorr ()            const { return this->rhoPUCorr_; };
       const float rhoPUCorrCalo ()        const { return this->rhoPUCorrCalo_; };
       const float rhoPUCorrCentralCalo () const { return this->rhoPUCorrCentralCalo_; };

       //////////////////////////////////////
       // Get track isolations
       //////////////////////////////////////

       const float trackIsoDR05 ()            const { return this->trackIsoDR05_; };

       float dEdxStrip()   const { return dEdxStrip_; }
       float dEdxStrip_Error                 () const { return this->dEdxStrip_Error_; };
       int dEdxStrip_nSaturatedMeasurements  () const { return this->dEdxStrip_numberOfSaturatedMeasurements_; };
       unsigned int dEdxStrip_nMeasurements  () const { return this->dEdxStrip_numberOfMeasurements_; };
       float dEdxPixel()  const { return dEdxPixel_; }
       float dEdxPixel_Error                 () const { return this->dEdxPixel_Error_; };
       int dEdxPixel_nSaturatedMeasurements  () const { return this->dEdxPixel_numberOfSaturatedMeasurements_; };
       unsigned int dEdxPixel_nMeasurements  () const { return this->dEdxPixel_numberOfMeasurements_; };

       // missing hits differentiated by location on track
       // re-implement these methods from osu::Track to provide a getter function when plotting osu::Track::matchedCandidateTrack()
       const unsigned char missingInnerHits_ ()  const { return this->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS); };
       const unsigned char missingMiddleHits_ () const { return this->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS); };
       const unsigned char missingOuterHits_ ()  const { return this->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS); };

     protected:
       PFIsolation pfIsolationDR03_;
       PFIsolation miniIsolation_;
       float matchedCaloJetEmEnergy_;  //energy of nearest calojet within a given dR;
       float matchedCaloJetHadEnergy_;
       bool pfLepOverlap_;
       float pfNeutralSum_;
       float dz_, dxy_, dzError_, dxyError_;
       int fromPV_;  //only stored for packedPFCandidates
       int trackQuality_;
       float dEdxStrip_, dEdxPixel_;  //in MeV/mm
       float dEdxStrip_Error_, dEdxPixel_Error_;
       int   dEdxStrip_numberOfSaturatedMeasurements_, dEdxPixel_numberOfSaturatedMeasurements_;
       unsigned int dEdxStrip_numberOfMeasurements_, dEdxPixel_numberOfMeasurements_;
   
       reco::HitPattern hitPattern_;
   
       std::vector<uint16_t> crossedEcalStatus_;
       std::vector<uint32_t> crossedHcalStatus_;
       int deltaEta_, deltaPhi_;
   
       PackedCandidateRef packedCandRef_;  // stored only for packedPFCands/lostTracks. NULL for generalTracks
       PackedCandidateRef nearestPFPackedCandRef_;
       PackedCandidateRef nearestLostTrackPackedCandRef_;
   
       float assocEMCaloDR05_;
       float assocHadCaloDR05_;
   
       float rhoPUCorr_;
       float rhoPUCorrCalo_;
       float rhoPUCorrCentralCalo_;
   
       float trackIsoDR05_;
       float trackIsoNoPUDR05_;
       float trackIsoNoFakesDR05_;
       float trackIsoNoPUNoFakesDR05_;
   
       const double getTrackIsolation (const reco::Track &, const std::vector<reco::Track> &, const double, const double = 1.0e-10) const;
   
       const double CaloTotDR05NoPU (RhoType, CaloType) const;

     };
   
     typedef std::vector<IsolatedTrack> IsolatedTrackCollection;

}  // namespace pat

#endif
