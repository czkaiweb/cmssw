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
    IsolatedTrack()
        : LeafCandidate(0, LorentzVector(0, 0, 0, 0)),
          pfIsolationDR03_(pat::PFIsolation()),
          miniIsolation_(pat::PFIsolation()),
          matchedCaloJetEmEnergy_(0.),
          matchedCaloJetHadEnergy_(0.),
          pfLepOverlap_(false),
          pfNeutralSum_(0.),
          dz_(0.),
          dxy_(0.),
          dzError_(0.),
          dxyError_(0.),
          fromPV_(-1),
          trackQuality_(0),
          dEdxStrip_(0),
          dEdxPixel_(0),
          hitPattern_(reco::HitPattern()),
          crossedEcalStatus_(std::vector<uint16_t>()),
          crossedHcalStatus_(std::vector<uint32_t>()),
          deltaEta_(0),
          deltaPhi_(0),
          packedCandRef_(PackedCandidateRef()),
          nearestPFPackedCandRef_(PackedCandidateRef()),
          nearestLostTrackPackedCandRef_(PackedCandidateRef()),
          assocEMCaloDR05_           (-1.),
          assocHadCaloDR05_          (-1.),
          rhoPUCorr_               (-1.),
          rhoPUCorrCalo_           (-1.),
          rhoPUCorrCentralCalo_    (-1.),
          trackIsoDR05_            (-1.),

    explicit IsolatedTrack(const PFIsolation& iso,
                           const PFIsolation& miniiso,
                           float caloJetEm,
                           float caloJetHad,
                           bool pfLepOverlap,
                           float pfNeutralSum,
                           const LorentzVector& p4,
                           int charge,
                           int id,
                           float dz,
                           float dxy,
                           float dzError,
                           float dxyError,
                           const reco::HitPattern& hp,
                           float dEdxS,
                           float dEdxP,
                           int fromPV,
                           int tkQual,
                           const std::vector<uint16_t>& ecalst,
                           const std::vector<uint32_t>& hcalst,
                           int dEta,
                           int dPhi,
                           const PackedCandidateRef& pcref,
                           const PackedCandidateRef& refToNearestPF,
                           const PackedCandidateRef& refToNearestLostTrack,
                           const reco::Track &track,
                           const vector<reco::Track> &tracks)
        : LeafCandidate(charge, p4, Point(0., 0., 0.), id),
          pfIsolationDR03_(iso),
          miniIsolation_(miniiso),
          matchedCaloJetEmEnergy_(caloJetEm),
          matchedCaloJetHadEnergy_(caloJetHad),
          pfLepOverlap_(pfLepOverlap),
          pfNeutralSum_(pfNeutralSum),
          dz_(dz),
          dxy_(dxy),
          dzError_(dzError),
          dxyError_(dxyError),
          fromPV_(fromPV),
          trackQuality_(tkQual),
          dEdxStrip_(dEdxS),
          dEdxPixel_(dEdxP),
          hitPattern_(hp),
          crossedEcalStatus_(ecalst),
          crossedHcalStatus_(hcalst),
          deltaEta_(dEta),
          deltaPhi_(dPhi),
          packedCandRef_(pcref),
          nearestPFPackedCandRef_(refToNearestPF),
          nearestLostTrackPackedCandRef_(refToNearestLostTrack),
          assocEMCaloDR05_           (-1.),
          assocHadCaloDR05_          (-1.),
          rhoPUCorr_               (-1.),
          rhoPUCorrCalo_           (-1.),
          rhoPUCorrCentralCalo_    (-1.),
          trackIsoDRp5_            (getTrackIsolation (track, tracks, 0.5)),
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

    bool isHighPurityTrack() const {
      return (trackQuality_ & (1 << reco::TrackBase::highPurity)) >> reco::TrackBase::highPurity;
    }
    bool isTightTrack() const { return (trackQuality_ & (1 << reco::TrackBase::tight)) >> reco::TrackBase::tight; }
    bool isLooseTrack() const { return (trackQuality_ & (1 << reco::TrackBase::loose)) >> reco::TrackBase::loose; }

    const reco::HitPattern& hitPattern() const { return hitPattern_; }

    float dEdxStrip()   const { return dEdxStrip_; }
    float dEdxStrip_Error                 () const { return this->dEdxStrip_Error_; };
    int dEdxStrip_nSaturatedMeasurements  () const { return this->dEdxStrip_numberOfSaturatedMeasurements_; };
    unsigned int dEdx_Strip_nMeasurements () const { return this->dEdxStrip_numberOfMeasurements_; };
    float dEdxPixel() const { return dEdxPixel_; }
    float dEdxPixel_Error                 () const { return this->dEdxPixel_Error_; };
    int dEdxPixel_nSaturatedMeasurements  () const { return this->dEdxPixel_numberOfSaturatedMeasurements_; };
    unsigned int dEdx_Pixel_nMeasurements () const { return this->dEdxPixel_numberOfMeasurements_; };

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
    const float assocAllCaloDR05NoPU ()                   const { return CaloTotDR05NoPU(All, Sum); };
    const float assocAllEMCaloDR05NoPU ()                 const { return CaloTotDR05NoPU(All, EM); };
    const float assocAllHadCaloDR05NoPU ()                const { return CaloTotDR05NoPU(All, Had); };
     
    const float assocCaloDR05NoPUCalo ()               const { return CaloTotDR05NoPU(Calo, Sum); };
    const float assocCaloDR5NoPUCaloEm ()             const { return CaloTotDR05NoPU(Calo, EM); };
    const float assocCaloDR5NoPUCaloHad ()            const { return CaloTotDR05NoPU(Calo, Had); };

    const float caloDR5NoPUCentralCalo ()        const { return CaloTotDR05NoPU(CentralCalo, Sum); };
    const float caloDR5NoPUCentralCaloJEm ()     const { return CaloTotDR05NoPU(CentralCalo, EM); };
    const float caloDR5NoPUCentralCaloHad ()     const { return CaloTotDR05NoPU(CentralCalo, Had); };

    //////////////////////////////////////
    // Set calo energies
    //////////////////////////////////////

    void set_caloEMDR5 (double value) { caloEMDR5_  = value; };
    void set_caloHadDR5(double value) { caloHadDR5_ = value; };

    //////////////////////////////////////
    // Set rhos
    //////////////////////////////////////

    void set_rhoPUCorr  (double value) { rhoPUCorr_   = value; };
    void set_rhoPUCorrCalo         (double value) { rhoPUCorrCalo_   = value; };
    void set_rhoPUCorrCentralCalo  (double value) { rhoPUCorrCentralCalo_   = value; };

    //////////////////////////////////////
    // Set track isolations
    //////////////////////////////////////

    void set_trackIsoDRp5 (double value) { trackIsoDRp5_ = value; };

    void set_trackIsoNoPUDRp5 (double value) { trackIsoNoPUDRp5_ = value; };

    void set_trackIsoNoFakesDRp5 (double value) { trackIsoNoFakesDRp5_ = value; };

    void set_trackIsoNoPUNoFakesDRp5 (double value) { trackIsoNoPUNoFakesDRp5_ = value; };

    void set_dEdx_pixel (double value, 
                         double error, 
                         int nSaturatedMeasurements, 
                         unsigned int nMeasurements) { 
        dEdx_pixel_ = value; 
        dEdxError_pixel_ = error;
        dEdx_numberOfSaturatedMeasurements_pixel_ = nSaturatedMeasurements;
        dEdx_numberOfMeasurements_pixel_ = nMeasurements;
    };
    void set_dEdx_strip (double value, 
                        double error, 
                        int nSaturatedMeasurements, 
                        unsigned int nMeasurements) { 
        dEdx_strip_ = value; 
        dEdxError_strip_ = error;
        dEdx_numberOfSaturatedMeasurements_strip_ = nSaturatedMeasurements;
        dEdx_numberOfMeasurements_strip_ = nMeasurements;
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

    const float trackIsoDRp5 ()            const { return this->trackIsoDRp5_; };

    const float dEdx_pixel                      () const { return this->dEdx_pixel_; };
    const float dEdxError_pixel                 () const { return this->dEdxError_pixel_; };
    const int dEdx_nSaturatedMeasurements_pixel () const { return this->dEdx_numberOfSaturatedMeasurements_pixel_; };
    const unsigned int dEdx_nMeasurements_pixel () const { return this->dEdx_numberOfMeasurements_pixel_; };

    const float dEdx_strip                      () const { return this->dEdx_strip_; };
    const float dEdxError_strip                 () const { return this->dEdxError_strip_; };
    const int dEdx_nSaturatedMeasurements_strip () const { return this->dEdx_numberOfSaturatedMeasurements_strip_; };
    const unsigned int dEdx_nMeasurements_strip () const { return this->dEdx_numberOfMeasurements_strip_; };

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

    float caloEMDR5_;
    float caloHadDR5_;

    float rhoPUCorr_;
    float rhoPUCorrCalo_;
    float rhoPUCorrCentralCalo_;

    float trackIsoDRp5_;

    float dEdx_pixel_;
    float dEdxError_pixel_;
    int dEdx_numberOfSaturatedMeasurements_pixel_;
    unsigned int dEdx_numberOfMeasurements_pixel_;

    float dEdx_strip_;
    float dEdxError_strip_;
    int dEdx_numberOfSaturatedMeasurements_strip_;
    unsigned int dEdx_numberOfMeasurements_strip_;

    const double getTrackIsolation (const reco::Track &, const vector<reco::Track> &, const bool, const bool, const double, const double = 1.0e-10) const;

    const double caloTotDR5NoPU (RhoType, CaloType) const;
  };

  typedef std::vector<IsolatedTrack> IsolatedTrackCollection;

}  // namespace pat

#endif
