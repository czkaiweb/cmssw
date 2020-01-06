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

#define INVALID_VALUE (numeric_limits<int>::min ())

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
          caloNewEMDRp5_           (INVALID_VALUE),
          caloNewHadDRp5_          (INVALID_VALUE),
          caloNewEMDRp3_           (INVALID_VALUE),
          caloNewHadDRp3_          (INVALID_VALUE),
          caloNewEMDRp2_           (INVALID_VALUE),
          caloNewHadDRp2_          (INVALID_VALUE),
          caloNewEMDRp1_           (INVALID_VALUE),
          caloNewHadDRp1_          (INVALID_VALUE),
          rhoPUCorr_               (INVALID_VALUE),
          rhoPUCorrCalo_           (INVALID_VALUE),
          rhoPUCorrCentralCalo_    (INVALID_VALUE),
          trackIsoDRp5_            (INVALID_VALUE),
          trackIsoDRp3_            (INVALID_VALUE),
          trackIsoDRp2_            (INVALID_VALUE),
          trackIsoDRp1_            (INVALID_VALUE),
          trackIsoNoPUDRp5_        (INVALID_VALUE),
          trackIsoNoPUDRp3_        (INVALID_VALUE),
          trackIsoNoPUDRp2_        (INVALID_VALUE),
	  trackIsoNoPUDRp1_        (INVALID_VALUE),
	  trackIsoNoFakesDRp5_     (INVALID_VALUE),
	  trackIsoNoFakesDRp3_     (INVALID_VALUE),
	  trackIsoNoFakesDRp2_     (INVALID_VALUE),
	  trackIsoNoFakesDRp1_     (INVALID_VALUE),
	  trackIsoNoPUNoFakesDRp5_ (INVALID_VALUE),
	  trackIsoNoPUNoFakesDRp3_ (INVALID_VALUE),
	  trackIsoNoPUNoFakesDRp2_ (INVALID_VALUE),
	  trackIsoNoPUNoFakesDRp1_ (INVALID_VALUE) {}

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
                           const PackedCandidateRef& refToNearestLostTrack)
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
          caloNewEMDRp5_           (INVALID_VALUE),
          caloNewHadDRp5_          (INVALID_VALUE),
          caloNewEMDRp3_           (INVALID_VALUE),
          caloNewHadDRp3_          (INVALID_VALUE),
          caloNewEMDRp2_           (INVALID_VALUE),
          caloNewHadDRp2_          (INVALID_VALUE),
          caloNewEMDRp1_           (INVALID_VALUE),
          caloNewHadDRp1_          (INVALID_VALUE),
          rhoPUCorr_               (INVALID_VALUE),
          rhoPUCorrCalo_           (INVALID_VALUE),
          rhoPUCorrCentralCalo_    (INVALID_VALUE),
          trackIsoDRp5_            (INVALID_VALUE),
          trackIsoDRp3_            (INVALID_VALUE),
          trackIsoDRp2_            (INVALID_VALUE),
          trackIsoDRp1_            (INVALID_VALUE),
          trackIsoNoPUDRp5_        (INVALID_VALUE),
          trackIsoNoPUDRp3_        (INVALID_VALUE),
          trackIsoNoPUDRp2_        (INVALID_VALUE),
          trackIsoNoPUDRp1_        (INVALID_VALUE),
          trackIsoNoFakesDRp5_     (INVALID_VALUE),
          trackIsoNoFakesDRp3_     (INVALID_VALUE),
          trackIsoNoFakesDRp2_     (INVALID_VALUE),
          trackIsoNoFakesDRp1_     (INVALID_VALUE),
          trackIsoNoPUNoFakesDRp5_ (INVALID_VALUE),
          trackIsoNoPUNoFakesDRp3_ (INVALID_VALUE),
          trackIsoNoPUNoFakesDRp2_ (INVALID_VALUE),
          trackIsoNoPUNoFakesDRp1_ (INVALID_VALUE) {}

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
    const float caloNewEMDRp5 ()  const { return this->caloNewEMDRp5_; };
    const float caloNewHadDRp5 () const { return this->caloNewHadDRp5_; };
    const float caloNewDRp5 ()    const { return this->caloNewEMDRp5_ + this->caloNewHadDRp5_; };

    // dR < 0.3
    const float caloNewEMDRp3 ()  const { return this->caloNewEMDRp3_; };
    const float caloNewHadDRp3 () const { return this->caloNewHadDRp3_; };
    const float caloNewDRp3 ()    const { return this->caloNewEMDRp3_ + this->caloNewHadDRp3_; };

    // dR < 0.2
    const float caloNewEMDRp2 ()  const { return this->caloNewEMDRp2_; };
    const float caloNewHadDRp2 () const { return this->caloNewHadDRp2_; };
    const float caloNewDRp2 ()    const { return this->caloNewEMDRp2_ + this->caloNewHadDRp2_; };

    // dR < 0.1
    const float caloNewEMDRp1 ()  const { return this->caloNewEMDRp1_; };
    const float caloNewHadDRp1 () const { return this->caloNewHadDRp1_; };
    const float caloNewDRp1 ()    const { return this->caloNewEMDRp1_ + this->caloNewHadDRp1_; };

    //////////////////////////////////////
    // Rho-corrected calo energies
    //////////////////////////////////////

    // New calculation that uses all rec hits in dR < 0.5 cone.
    const float caloNewNoPUDRp5 ()                   const { return caloTotNoPU(0.5, All, Sum); };
    const float caloNewNoPUDRp5JustEm ()             const { return caloTotNoPU(0.5, All, EM); };
    const float caloNewNoPUDRp5JustHad ()            const { return caloTotNoPU(0.5, All, Had); };
    
    const float caloNewNoPUDRp5Calo ()               const { return caloTotNoPU(0.5, Calo, Sum); };
    const float caloNewNoPUDRp5CaloJustEm ()         const { return caloTotNoPU(0.5, Calo, EM); };
    const float caloNewNoPUDRp5CaloJustHad ()        const { return caloTotNoPU(0.5, Calo, Had); };

    const float caloNewNoPUDRp5CentralCalo ()        const { return caloTotNoPU(0.5, CentralCalo, Sum); };
    const float caloNewNoPUDRp5CentralCaloJustEm ()  const { return caloTotNoPU(0.5, CentralCalo, EM); };
    const float caloNewNoPUDRp5CentralCaloJustHad () const { return caloTotNoPU(0.5, CentralCalo, Had); };

    // dR < 0.3
    const float caloNewNoPUDRp3 ()                   const { return caloTotNoPU(0.3, All, Sum); };
    const float caloNewNoPUDRp3JustEm ()             const { return caloTotNoPU(0.3, All, EM); };
    const float caloNewNoPUDRp3JustHad ()            const { return caloTotNoPU(0.3, All, Had); };
    
    const float caloNewNoPUDRp3Calo ()               const { return caloTotNoPU(0.3, Calo, Sum); };
    const float caloNewNoPUDRp3CaloJustEm ()         const { return caloTotNoPU(0.3, Calo, EM); };
    const float caloNewNoPUDRp3CaloJustHad ()        const { return caloTotNoPU(0.3, Calo, Had); };

    const float caloNewNoPUDRp3CentralCalo ()        const { return caloTotNoPU(0.3, CentralCalo, Sum); };
    const float caloNewNoPUDRp3CentralCaloJustEm ()  const { return caloTotNoPU(0.3, CentralCalo, EM); };
    const float caloNewNoPUDRp3CentralCaloJustHad () const { return caloTotNoPU(0.3, CentralCalo, Had); };

    // dR < 0.2
    const float caloNewNoPUDRp2 ()                   const { return caloTotNoPU(0.2, All, Sum); };
    const float caloNewNoPUDRp2JustEm ()             const { return caloTotNoPU(0.2, All, EM); };
    const float caloNewNoPUDRp2JustHad ()            const { return caloTotNoPU(0.2, All, Had); };
    
    const float caloNewNoPUDRp2Calo ()               const { return caloTotNoPU(0.2, Calo, Sum); };
    const float caloNewNoPUDRp2CaloJustEm ()         const { return caloTotNoPU(0.2, Calo, EM); };
    const float caloNewNoPUDRp2CaloJustHad ()        const { return caloTotNoPU(0.2, Calo, Had); };

    const float caloNewNoPUDRp2CentralCalo ()        const { return caloTotNoPU(0.2, CentralCalo, Sum); };
    const float caloNewNoPUDRp2CentralCaloJustEm ()  const { return caloTotNoPU(0.2, CentralCalo, EM); };
    const float caloNewNoPUDRp2CentralCaloJustHad () const { return caloTotNoPU(0.2, CentralCalo, Had); };

    // dR < 0.1
    const float caloNewNoPUDRp1 ()                   const { return caloTotNoPU(0.1, All, Sum); };
    const float caloNewNoPUDRp1JustEm ()             const { return caloTotNoPU(0.1, All, EM); };
    const float caloNewNoPUDRp1JustHad ()            const { return caloTotNoPU(0.1, All, Had); };
    
    const float caloNewNoPUDRp1Calo ()               const { return caloTotNoPU(0.1, Calo, Sum); };
    const float caloNewNoPUDRp1CaloJustEm ()         const { return caloTotNoPU(0.1, Calo, EM); };
    const float caloNewNoPUDRp1CaloJustHad ()        const { return caloTotNoPU(0.1, Calo, Had); };

    const float caloNewNoPUDRp1CentralCalo ()        const { return caloTotNoPU(0.1, CentralCalo, Sum); };
    const float caloNewNoPUDRp1CentralCaloJustEm ()  const { return caloTotNoPU(0.1, CentralCalo, EM); };
    const float caloNewNoPUDRp1CentralCaloJustHad () const { return caloTotNoPU(0.1, CentralCalo, Had); };

    //////////////////////////////////////
    // Set calo energies
    //////////////////////////////////////

    void set_caloNewEMDRp5 (double value) { caloNewEMDRp5_  = value; };
    void set_caloNewHadDRp5(double value) { caloNewHadDRp5_ = value; };

    void set_caloNewEMDRp3 (double value) { caloNewEMDRp3_  = value; };
    void set_caloNewHadDRp3(double value) { caloNewHadDRp3_ = value; };

    void set_caloNewEMDRp2 (double value) { caloNewEMDRp2_  = value; };
    void set_caloNewHadDRp2(double value) { caloNewHadDRp2_ = value; };

    void set_caloNewEMDRp1 (double value) { caloNewEMDRp1_  = value; };
    void set_caloNewHadDRp1(double value) { caloNewHadDRp1_ = value; };

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
    void set_trackIsoDRp3 (double value) { trackIsoDRp3_ = value; };
    void set_trackIsoDRp2 (double value) { trackIsoDRp2_ = value; };
    void set_trackIsoDRp1 (double value) { trackIsoDRp1_ = value; };

    void set_trackIsoNoPUDRp5 (double value) { trackIsoNoPUDRp5_ = value; };
    void set_trackIsoNoPUDRp3 (double value) { trackIsoNoPUDRp3_ = value; };
    void set_trackIsoNoPUDRp2 (double value) { trackIsoNoPUDRp2_ = value; };
    void set_trackIsoNoPUDRp1 (double value) { trackIsoNoPUDRp1_ = value; };

    void set_trackIsoNoFakesDRp5 (double value) { trackIsoNoFakesDRp5_ = value; };
    void set_trackIsoNoFakesDRp3 (double value) { trackIsoNoFakesDRp3_ = value; };
    void set_trackIsoNoFakesDRp2 (double value) { trackIsoNoFakesDRp2_ = value; };
    void set_trackIsoNoFakesDRp1 (double value) { trackIsoNoFakesDRp1_ = value; };

    void set_trackIsoNoPUNoFakesDRp5 (double value) { trackIsoNoPUNoFakesDRp5_ = value; };
    void set_trackIsoNoPUNoFakesDRp3 (double value) { trackIsoNoPUNoFakesDRp3_ = value; };
    void set_trackIsoNoPUNoFakesDRp2 (double value) { trackIsoNoPUNoFakesDRp2_ = value; };
    void set_trackIsoNoPUNoFakesDRp1 (double value) { trackIsoNoPUNoFakesDRp1_ = value; };

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
    const float trackIsoDRp3 ()            const { return this->trackIsoDRp3_; };
    const float trackIsoDRp2 ()            const { return this->trackIsoDRp2_; };
    const float trackIsoDRp1 ()            const { return this->trackIsoDRp1_; };

    const float trackIsoNoPUDRp5 ()        const { return this->trackIsoNoPUDRp5_; };
    const float trackIsoNoPUDRp3 ()        const { return this->trackIsoNoPUDRp3_; };
    const float trackIsoNoPUDRp2 ()        const { return this->trackIsoNoPUDRp2_; };
    const float trackIsoNoPUDRp1 ()        const { return this->trackIsoNoPUDRp1_; };
    
    const float trackIsoNoFakesDRp5 ()     const { return this->trackIsoNoFakesDRp5_; };
    const float trackIsoNoFakesDRp3 ()     const { return this->trackIsoNoFakesDRp3_; };
    const float trackIsoNoFakesDRp2 ()     const { return this->trackIsoNoFakesDRp2_; };
    const float trackIsoNoFakesDRp1 ()     const { return this->trackIsoNoFakesDRp1_; };
    
    const float trackIsoNoPUNoFakesDRp5 () const { return this->trackIsoNoPUNoFakesDRp5_; };
    const float trackIsoNoPUNoFakesDRp3 () const { return this->trackIsoNoPUNoFakesDRp3_; };
    const float trackIsoNoPUNoFakesDRp2 () const { return this->trackIsoNoPUNoFakesDRp2_; };
    const float trackIsoNoPUNoFakesDRp1 () const { return this->trackIsoNoPUNoFakesDRp1_; };

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

    float caloNewEMDRp5_;
    float caloNewHadDRp5_;

    float caloNewEMDRp3_;
    float caloNewHadDRp3_;

    float caloNewEMDRp2_;
    float caloNewHadDRp2_;

    float caloNewEMDRp1_;
    float caloNewHadDRp1_;

    float rhoPUCorr_;
    float rhoPUCorrCalo_;
    float rhoPUCorrCentralCalo_;

    float trackIsoDRp5_;
    float trackIsoDRp3_;
    float trackIsoDRp2_;
    float trackIsoDRp1_;

    float trackIsoNoPUDRp5_;
    float trackIsoNoPUDRp3_;
    float trackIsoNoPUDRp2_;
    float trackIsoNoPUDRp1_;

    float trackIsoNoFakesDRp5_;
    float trackIsoNoFakesDRp3_;
    float trackIsoNoFakesDRp2_;
    float trackIsoNoFakesDRp1_;

    float trackIsoNoPUNoFakesDRp5_;
    float trackIsoNoPUNoFakesDRp3_;
    float trackIsoNoPUNoFakesDRp2_;
    float trackIsoNoPUNoFakesDRp1_;

    float dEdx_pixel_;
    float dEdxError_pixel_;
    int dEdx_numberOfSaturatedMeasurements_pixel_;
    unsigned int dEdx_numberOfMeasurements_pixel_;

    float dEdx_strip_;
    float dEdxError_strip_;
    int dEdx_numberOfSaturatedMeasurements_strip_;
    unsigned int dEdx_numberOfMeasurements_strip_;

    const double getTrackIsolation (const reco::Track &track, const vector<reco::Track> &tracks, const bool noPU = true, const bool noFake = true, const double outerDeltaR, const double innerDeltaR= 1.0e-12) const {

      double sumPt = 0.0;

      for (const auto &t : tracks) {
        if (noFakes && t.normalizedChi2 () > 20.0)
          continue;
        if (noFakes && t.hitPattern ().pixelLayersWithMeasurement () < 2)
          continue;
        if (noFakes && t.hitPattern ().trackerLayersWithMeasurement () < 5)
          continue;
        if (noFakes && fabs (t.d0 () / t.d0Error ()) > 5.0)
          continue;

        if (noPU && fabs( track.dz (t.vertex ())) > 3.0 * hypot (track.dzError (), t.dzError ()))
          continue;

        double dR = deltaR (track, t);
        if (dR < outerDeltaR && dR > innerDeltaR)
          sumPt += t.pt ();
        }   

      return sumPt;
      };

    const double caloTotNoPU (double, RhoType = All, CaloType = Sum) const{
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
          throw cms::Exception("FatalError") << "Unkown or not implemented rho type requested, type:" << rhoType;
        } 
  
      double rawCaloTot = 0.0;
      int intDR = dR * 10.0;
      switch (caloType) {
        case Sum:
          if(intDR == 5) rawCaloTot = caloNewDRp5();
          else if(intDR == 3) rawCaloTot = caloNewDRp3();
          else if(intDR == 2) rawCaloTot = caloNewDRp2();
          else if(intDR == 1) rawCaloTot = caloNewDRp1();
          break;
        case EM:
          if(intDR == 5) rawCaloTot = caloNewEMDRp5();
          else if(intDR == 3) rawCaloTot = caloNewEMDRp3();
          else if(intDR == 2) rawCaloTot = caloNewEMDRp2();
          else if(intDR == 1) rawCaloTot = caloNewEMDRp1();
          break;
        case Had:
	  if(intDR == 5) rawCaloTot = caloNewHadDRp5();
	  else if(intDR == 3) rawCaloTot = caloNewHadDRp3();
	  else if(intDR == 2) rawCaloTot = caloNewHadDRp2();
	  else if(intDR == 1) rawCaloTot = caloNewHadDRp1();
	  break; 
	}	
  
      double caloCorr = rho * TMath::Pi() * dR * dR;  // Define effective area as pi*r^2, where r is radius of DeltaR cone.
      double caloTotNoPU = TMath::Max(0., rawCaloTot - caloCorr);
      return caloTotNoPU;
    };
  };

  typedef std::vector<IsolatedTrack> IsolatedTrackCollection;

}  // namespace pat

#endif
