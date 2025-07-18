#ifndef Validation_HGCalValidation_BarrelVHistoProducerAlgo_h
#define Validation_HGCalValidation_BarrelVHistoProducerAlgo_h

/* \author HGCal
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "CommonTools/RecoAlgos/interface/MultiVectorManager.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociatorBaseImpl.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociatorBaseImpl.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

struct BarrelVHistoProducerAlgoHistograms {
  dqm::reco::MonitorElement* lastLayerEB;
  dqm::reco::MonitorElement* lastLayerHB;
  //1D
  std::vector<dqm::reco::MonitorElement*> h_cluster_eta;
  std::vector<dqm::reco::MonitorElement*> h_energyclustered;

  std::unordered_map<int, dqm::reco::MonitorElement*> h_clusternum_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_energyclustered_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_score_layercl2caloparticle_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_score_caloparticle2layercl_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_energy_vs_score_caloparticle2layercl_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_energy_vs_score_layercl2caloparticle_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2layercl_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2layercl_vs_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2layercl_vs_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_layercl2caloparticle_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_layercl2caloparticle_vs_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_sharedenergy_layercl2caloparticle_vs_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_num_caloparticle_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_numDup_caloparticle_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_denom_caloparticle_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_num_caloparticle_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_numDup_caloparticle_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_denom_caloparticle_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_num_layercl_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_numMerge_layercl_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_denom_layercl_eta_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_num_layercl_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_numMerge_layercl_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_denom_layercl_phi_perlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_cellAssociation_perlayer;

  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_eta;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_eta_Zorigin;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_energy;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_selfenergy;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_energyDifference;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_pt;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_phi;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nSimClusters;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nHitsInSimClusters;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_firstlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_lastlayer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_layersnum;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nHitsInSimClusters_matchedtoRecHit;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nHits_matched_energy;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nHits_matched_energy_layer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_nHits_matched_energy_layer_1SimCl;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_sum_energy_layer;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_firstlayer_matchedtoRecHit;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_lastlayer_matchedtoRecHit;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_layersnum_matchedtoRecHit;
  std::unordered_map<int, dqm::reco::MonitorElement*> h_caloparticle_fractions, h_caloparticle_fractions_weight;

  //For SimClusters
  std::unordered_map<int, dqm::reco::MonitorElement*> h_simclusternum_perlayer;

  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_denom_layercl_in_simcl_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_denom_layercl_in_simcl_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_score_layercl2simcluster_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_layercl2simcluster_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_energy_vs_score_layercl2simcluster_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_num_layercl_in_simcl_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_num_layercl_in_simcl_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_numMerge_layercl_in_simcl_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_numMerge_layercl_in_simcl_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_layercl2simcluster_vs_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_layercl2simcluster_vs_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_denom_simcluster_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_denom_simcluster_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_score_simcluster2layercl_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_simcluster2layercl_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_energy_vs_score_simcluster2layercl_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_num_simcluster_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_num_simcluster_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_numDup_simcluster_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_numDup_simcluster_phi_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_simcluster2layercl_vs_eta_perlayer;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_sharedenergy_simcluster2layercl_vs_phi_perlayer;

  // For Tracksters
  constexpr static int numberOfValidationTypes_ = 4;

  std::vector<dqm::reco::MonitorElement*> h_score_trackster2caloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_score_trackster2bestCaloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_score_trackster2bestCaloparticle2[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_score_caloparticle2trackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_scorePur_caloparticle2trackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_scoreDupl_caloparticle2trackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_trackster2caloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_trackster2bestCaloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_trackster2bestCaloparticle2[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_caloparticle2trackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_caloparticle2bestTrackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_energy_vs_score_caloparticle2bestTrackster2[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_trackster_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_trackster_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_trackster_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_trackster_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numMerge_trackster_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numMerge_trackster_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numMerge_trackster_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numMerge_trackster_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_trackster2caloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_trackster2bestCaloparticle[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_trackster2bestCaloparticle2[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2trackster[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2trackster_assoc[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2trackster_assoc2[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_trackster2bestCaloparticle_vs_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_trackster2bestCaloparticle_vs_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2trackster_assoc_vs_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_sharedenergy_caloparticle2trackster_assoc_vs_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_trackster_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_trackster_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_trackster_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_trackster_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numEff_caloparticle_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numEff_caloparticle_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numEff_caloparticle_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numEff_caloparticle_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_caloparticle_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_caloparticle_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_caloparticle_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_num_caloparticle_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numDup_trackster_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numDup_trackster_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numDup_trackster_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_numDup_trackster_pt[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_caloparticle_eta[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_caloparticle_phi[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_caloparticle_en[numberOfValidationTypes_];
  std::vector<dqm::reco::MonitorElement*> h_denom_caloparticle_pt[numberOfValidationTypes_];
  // Generic histograms
  std::vector<dqm::reco::MonitorElement*> h_tracksternum;
  std::vector<dqm::reco::MonitorElement*> h_conttracksternum;
  std::vector<dqm::reco::MonitorElement*> h_nonconttracksternum;
  std::vector<dqm::reco::MonitorElement*> h_clusternum_in_trackster;
  std::vector<std::unordered_map<int, dqm::reco::MonitorElement*>> h_clusternum_in_trackster_perlayer;
  std::vector<dqm::reco::MonitorElement*> h_multiplicityOfLCinTST;
  std::vector<dqm::reco::MonitorElement*> h_multiplicity_numberOfEventsHistogram;
  std::vector<dqm::reco::MonitorElement*> h_multiplicity_zminus_numberOfEventsHistogram;
  std::vector<dqm::reco::MonitorElement*> h_multiplicity_zplus_numberOfEventsHistogram;
  std::vector<dqm::reco::MonitorElement*> h_multiplicityOfLCinTST_vs_layercluster;
  std::vector<dqm::reco::MonitorElement*> h_multiplicityOfLCinTST_vs_layerclusterenergy;
  std::vector<dqm::reco::MonitorElement*> h_clusternum_in_trackster_vs_layer;
  std::vector<dqm::reco::MonitorElement*> h_trackster_pt;
  std::vector<dqm::reco::MonitorElement*> h_trackster_eta;
  std::vector<dqm::reco::MonitorElement*> h_trackster_phi;
  std::vector<dqm::reco::MonitorElement*> h_trackster_energy;
  std::vector<dqm::reco::MonitorElement*> h_trackster_x;
  std::vector<dqm::reco::MonitorElement*> h_trackster_y;
  std::vector<dqm::reco::MonitorElement*> h_trackster_z;
  std::vector<dqm::reco::MonitorElement*> h_trackster_firstlayer;
  std::vector<dqm::reco::MonitorElement*> h_trackster_lastlayer;
  std::vector<dqm::reco::MonitorElement*> h_trackster_layersnum;
};

using Density = hgcal_clustering::Density;

class BarrelVHistoProducerAlgo {
public:
  typedef dqm::legacy::DQMStore DQMStore;
  typedef dqm::legacy::MonitorElement MonitorElement;
  using TracksterToTracksterMap =
      ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>;
  using SimClusterToCaloParticleMap =
      ticl::AssociationMap<ticl::oneToOneMapWithFraction, std::vector<SimCluster>, std::vector<CaloParticle>>;
  enum validationType { byHits_CP = 0, byLCs, byLCs_CP, byHits };

  BarrelVHistoProducerAlgo(const edm::ParameterSet& pset);
  ~BarrelVHistoProducerAlgo();

  using Histograms = BarrelVHistoProducerAlgoHistograms;

  void bookInfo(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookCaloParticleHistos(DQMStore::IBooker& ibook, Histograms& histograms, int pdgid, unsigned int layers);

  void bookSimClusterHistos(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookSimClusterAssociationHistos(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookClusterHistos_ClusterLevel(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookClusterHistos_LCtoCP_association(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookClusterHistos_CellLevel(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookTracksterHistos(DQMStore::IBooker& ibook, Histograms& histograms, unsigned int layers);

  void bookTracksterSTSHistos(DQMStore::IBooker& ibook, Histograms& histograms, const validationType valType);

  void layerClusters_to_CaloParticles(const Histograms& histograms,
                                      edm::Handle<reco::CaloClusterCollection> clusterHandle,
                                      const reco::CaloClusterCollection& clusters,
                                      edm::Handle<std::vector<CaloParticle>> caloParticleHandle,
                                      std::vector<CaloParticle> const& cP,
                                      std::vector<size_t> const& cPIndices,
                                      std::vector<size_t> const& cPSelectedIndices,
                                      std::unordered_map<DetId, const unsigned int> const&,
                                      unsigned int layers,
                                      const ticl::RecoToSimCollection& recSimColl,
                                      const ticl::SimToRecoCollection& simRecColl,
                                      MultiVectorManager<reco::PFRecHit> const& barrelHits) const;
  void layerClusters_to_SimClusters(const Histograms& histograms,
                                    const int count,
                                    edm::Handle<reco::CaloClusterCollection> clusterHandle,
                                    const reco::CaloClusterCollection& clusters,
                                    edm::Handle<std::vector<SimCluster>> simClusterHandle,
                                    std::vector<SimCluster> const& simClusters,
                                    std::vector<size_t> const& sCIndices,
                                    const std::vector<float>& mask,
                                    std::unordered_map<DetId, const unsigned int> const&,
                                    unsigned int layers,
                                    const ticl::RecoToSimCollectionWithSimClusters& recSimColl,
                                    const ticl::SimToRecoCollectionWithSimClusters& simRecColl,
                                    MultiVectorManager<reco::PFRecHit> const& barrelHits) const;

  void tracksters_to_SimTracksters_fp(const Histograms& histograms,
                                      const int count,
                                      const TracksterToTracksterMap& trackstersToSimTrackstersMap,
                                      const TracksterToTracksterMap& simTrackstersToTrackstersMap,
                                      const validationType valType,
                                      const SimClusterToCaloParticleMap& scToCpMap,
                                      const std::vector<size_t>& cPIndices,
                                      const std::vector<size_t>& cPSelectedIndices,
                                      const edm::ProductID& cPHandle_id) const;

  void fill_info_histos(const Histograms& histograms, unsigned int layers) const;
  void fill_caloparticle_histos(const Histograms& histograms,
                                int pdgid,
                                const CaloParticle& caloparticle,
                                std::vector<SimVertex> const& simVertices,
                                unsigned int layers,
                                std::unordered_map<DetId, const unsigned int> const&,
                                MultiVectorManager<reco::PFRecHit> const& barrelHits) const;
  void fill_generic_cluster_histos(const Histograms& histograms,
                                   const int count,
                                   edm::Handle<reco::CaloClusterCollection> clusterHandle,
                                   const reco::CaloClusterCollection& clusters,
                                   edm::Handle<std::vector<CaloParticle>> caloParticleHandle,
                                   std::vector<CaloParticle> const& cP,
                                   std::vector<size_t> const& cPIndices,
                                   std::vector<size_t> const& cPSelectedIndices,
                                   std::unordered_map<DetId, const unsigned int> const&,
                                   unsigned int layers,
                                   const ticl::RecoToSimCollection& recSimColl,
                                   const ticl::SimToRecoCollection& simRecColl,
                                   MultiVectorManager<reco::PFRecHit> const& barrelHits) const;
  void fill_simCluster_histos(const Histograms& histograms,
                              std::vector<SimCluster> const& simClusters,
                              unsigned int layers) const;
  void fill_simClusterAssociation_histos(const Histograms& histograms,
                                         const int count,
                                         edm::Handle<reco::CaloClusterCollection> clusterHandle,
                                         const reco::CaloClusterCollection& clusters,
                                         edm::Handle<std::vector<SimCluster>> simClusterHandle,
                                         std::vector<SimCluster> const& simClusters,
                                         std::vector<size_t> const& sCIndices,
                                         const std::vector<float>& mask,
                                         std::unordered_map<DetId, const unsigned int> const& barrelHitMap,
                                         unsigned int layers,
                                         const ticl::RecoToSimCollectionWithSimClusters& recSimColl,
                                         const ticl::SimToRecoCollectionWithSimClusters& simRecColl,
                                         MultiVectorManager<reco::PFRecHit> const& barrelHits) const;
  void fill_cluster_histos(const Histograms& histograms, const int count, const reco::CaloCluster& cluster) const;

  double distance2(const double x1, const double y1, const double x2, const double y2) const;
  double distance(const double x1, const double y1, const double x2, const double y2) const;

  void setRecHitTools(std::shared_ptr<hgcal::RecHitTools> recHitTools);

  DetId findmaxhit(const reco::CaloCluster& cluster,
                   std::unordered_map<DetId, const unsigned int> const&,
                   MultiVectorManager<reco::PFRecHit> const& hits) const;

  struct detIdInfoInCluster {
    bool operator==(const detIdInfoInCluster& o) const { return clusterId == o.clusterId; };
    long unsigned int clusterId;
    float fraction;
  };

  struct detIdInfoInTrackster {
    bool operator==(const detIdInfoInTrackster& o) const { return tracksterId == o.tracksterId; };
    unsigned int tracksterId;
    long unsigned int clusterId;
    float fraction;
  };

  struct caloParticleOnLayer {
    unsigned int caloParticleId;
    float energy = 0;
    std::vector<std::pair<DetId, float>> hits_and_fractions;
    std::unordered_map<unsigned int, std::pair<float, float>> layerClusterIdToEnergyAndScore;
  };

private:
  double getEta(double eta) const;

  std::shared_ptr<hgcal::RecHitTools> recHitTools_;
  constexpr static int numberOfValidationTypes_ = 4;
  std::array<std::string, numberOfValidationTypes_> ref_ = {
      {"SimTrackster_fromCP_byHits", "SimTrackster_byLCs", "SimTrackster_fromCP_byLCs", "SimTrackster_byHits"}};
  std::array<std::string, numberOfValidationTypes_> refText_ = {{"SimTrackster from CP Associated by Hits",
                                                                 "SimTrackster Associated by LCs",
                                                                 "SimTrackster from CP Associated by LCs",
                                                                 "SimTrackster Associated by Hits"}};
  // Must be in sync with labels in PostProcessorHGCAL_cfi.py
  std::array<std::string, numberOfValidationTypes_> valSuffix_ = {{"_byHits_CP", "_byLCs", "_byLCs_CP", "_byHits"}};

  int barrelLayersOffset_ = 5;

  //private data members
  double minEta_, maxEta_;
  int nintEta_;
  bool useFabsEta_;
  double minEne_, maxEne_;
  int nintEne_;
  double minPt_, maxPt_;
  int nintPt_;
  double minPhi_, maxPhi_;
  int nintPhi_;
  double minMixedHitsSimCluster_, maxMixedHitsSimCluster_;
  int nintMixedHitsSimCluster_;
  double minMixedHitsCluster_, maxMixedHitsCluster_;
  int nintMixedHitsCluster_;
  double minEneCl_, maxEneCl_;
  int nintEneCl_;
  double minLongDepBary_, maxLongDepBary_;
  int nintLongDepBary_;
  double minZpos_, maxZpos_;
  int nintZpos_;
  double minTotNsimClsperlay_, maxTotNsimClsperlay_;
  int nintTotNsimClsperlay_;
  double minTotNClsperlay_, maxTotNClsperlay_;
  int nintTotNClsperlay_;
  double minEneClperlay_, maxEneClperlay_;
  int nintEneClperlay_;
  double minScore_, maxScore_;
  int nintScore_;
  double minSharedEneFrac_, maxSharedEneFrac_;
  int nintSharedEneFrac_;
  double minTSTSharedEneFracEfficiency_;
  double minTSTSharedEneFrac_, maxTSTSharedEneFrac_;
  int nintTSTSharedEneFrac_;
  double minTotNTSTs_, maxTotNTSTs_;
  int nintTotNTSTs_;
  double minTotNClsinTSTs_, maxTotNClsinTSTs_;
  int nintTotNClsinTSTs_;
  double minTotNClsinTSTsperlayer_, maxTotNClsinTSTsperlayer_;
  int nintTotNClsinTSTsperlayer_;
  double minMplofLCs_, maxMplofLCs_;
  int nintMplofLCs_;
  double minSizeCLsinTSTs_, maxSizeCLsinTSTs_;
  int nintSizeCLsinTSTs_;
  double minClEnepermultiplicity_, maxClEnepermultiplicity_;
  int nintClEnepermultiplicity_;
  double minX_, maxX_;
  int nintX_;
  double minY_, maxY_;
  int nintY_;
  double minZ_, maxZ_;
  int nintZ_;
};

#endif
