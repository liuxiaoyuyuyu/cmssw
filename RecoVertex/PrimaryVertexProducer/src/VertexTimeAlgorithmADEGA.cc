#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "vdt/vdtMath.h"

#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmADEGA.h"

#ifdef PVTX_DEBUG
#define LOG edm::LogPrint("VertexTimeAlgorithmADEGA")
#else
#define LOG LogDebug("VertexTimeAlgorithmADEGA")
#endif

VertexTimeAlgorithmADEGA::VertexTimeAlgorithmADEGA(edm::ParameterSet const& iConfig,
                                                                   edm::ConsumesCollector& iCC)
    : VertexTimeAlgorithmBase(iConfig, iCC),
      trackMTDTimeToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeVMapTag"))),
      trackMTDTimeErrorToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeErrorVMapTag"))),
      trackMTDTimeQualityToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
      trackMTDTofPiToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofPiVMapTag"))),
      trackMTDTofKToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofKVMapTag"))),
      trackMTDTofPToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTofPVMapTag"))),
      trackMTDSigmaTofPiToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofPiVMapTag"))),
      trackMTDSigmaTofKToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofKVMapTag"))),
      trackMTDSigmaTofPToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDSigmaTofPVMapTag"))),
      minTrackVtxWeight_(iConfig.getParameter<double>("minTrackVtxWeight")),
      minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")),
      probPion_(iConfig.getParameter<double>("probPion")),
      probKaon_(iConfig.getParameter<double>("probKaon")),
      probProton_(iConfig.getParameter<double>("probProton")),
      Tstart_(iConfig.getParameter<double>("Tstart")),
      coolingFactor_(iConfig.getParameter<double>("coolingFactor")),
      populationSize_(15) {}
      

void VertexTimeAlgorithmADEGA::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
  VertexTimeAlgorithmBase::fillPSetDescription(iDesc);

  iDesc.add<edm::InputTag>("trackMTDTimeVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))
      ->setComment("Input ValueMap for track time at MTD");
  iDesc.add<edm::InputTag>("trackMTDTimeErrorVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))
      ->setComment("Input ValueMap for track time uncertainty at MTD");
  iDesc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))
      ->setComment("Input ValueMap for track MVA quality value");
  iDesc.add<edm::InputTag>("trackMTDTofPiVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"))
      ->setComment("Input ValueMap for track tof as pion");
  iDesc.add<edm::InputTag>("trackMTDTofKVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"))
      ->setComment("Input ValueMap for track tof as kaon");
  iDesc.add<edm::InputTag>("trackMTDTofPVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"))
      ->setComment("Input ValueMap for track tof as proton");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofPiVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"))
      ->setComment("Input ValueMap for track tof uncertainty as pion");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofKVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"))
      ->setComment("Input ValueMap for track tof uncertainty as kaon");
  iDesc.add<edm::InputTag>("trackMTDSigmaTofPVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"))
      ->setComment("Input ValueMap for track tof uncertainty as proton");

  iDesc.add<double>("minTrackVtxWeight", 0.5)->setComment("Minimum track weight");
  iDesc.add<double>("minTrackTimeQuality", 0.8)->setComment("Minimum MVA Quality selection on tracks");

  iDesc.add<double>("probPion", 0.7)->setComment("A priori probability pions");
  iDesc.add<double>("probKaon", 0.2)->setComment("A priori probability kaons");
  iDesc.add<double>("probProton", 0.1)->setComment("A priori probability protons");

  iDesc.add<double>("Tstart", 256.)->setComment("DA initial temperature T");
  iDesc.add<double>("coolingFactor", 0.5)->setComment("DA cooling factor");
}

void VertexTimeAlgorithmADEGA::setEvent(edm::Event& iEvent, edm::EventSetup const&) {
  // additional collections required for vertex-time calculation
  trackMTDTimes_ = iEvent.get(trackMTDTimeToken_);
  trackMTDTimeErrors_ = iEvent.get(trackMTDTimeErrorToken_);
  trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken_);
  trackMTDTofPi_ = iEvent.get(trackMTDTofPiToken_);
  trackMTDTofK_ = iEvent.get(trackMTDTofKToken_);
  trackMTDTofP_ = iEvent.get(trackMTDTofPToken_);
  trackMTDSigmaTofPi_ = iEvent.get(trackMTDSigmaTofPiToken_);
  trackMTDSigmaTofK_ = iEvent.get(trackMTDSigmaTofKToken_);
  trackMTDSigmaTofP_ = iEvent.get(trackMTDSigmaTofPToken_);
}

bool VertexTimeAlgorithmADEGA::vertexTime(float& vtxTime,
                                                  float& vtxTimeError,
                                                  const TransientVertex& vtx) const {
  if (vtx.originalTracks().empty()) {
    return false;
  }

  std::vector<TrackInfo> v_trackInfo;
  v_trackInfo.reserve(vtx.originalTracks().size());

  //Loop over the tracks associated with the vertex
  for (const auto& trk : vtx.originalTracks()) {
    auto const trkWeight = vtx.trackWeight(trk);
    if (trkWeight > minTrackVtxWeight_) {
      auto const trkTimeQuality = trackMTDTimeQualities_[trk.trackBaseRef()];

      if (trkTimeQuality >= minTrackTimeQuality_) {
        auto const trkTime = trackMTDTimes_[trk.trackBaseRef()];
        auto const trkTimeError = trackMTDTimeErrors_[trk.trackBaseRef()];

        v_trackInfo.emplace_back();//insert a new element at the end of the vector, 
        //effectively increases the container size by one.      
        auto& trkInfo = v_trackInfo.back();//Returns a reference to the last element in the vector

        //trkInfo.trkWeight = trkWeight;
        trkInfo.trkTimeHyp[0] = trkTime - trackMTDTofPi_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[1] = trkTime - trackMTDTofK_[trk.trackBaseRef()];
        trkInfo.trkTimeHyp[2] = trkTime - trackMTDTofP_[trk.trackBaseRef()];
        
        trkInfo.trkTimeErrorHyp[0] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofPi_[trk.trackBaseRef()] * trackMTDSigmaTofPi_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[1] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofK_[trk.trackBaseRef()] * trackMTDSigmaTofK_[trk.trackBaseRef()]);
        trkInfo.trkTimeErrorHyp[2] =
            std::sqrt(trkTimeError * trkTimeError +
                      trackMTDSigmaTofP_[trk.trackBaseRef()] * trackMTDSigmaTofP_[trk.trackBaseRef()]);

        trkInfo.weight[0]=trkWeight / (trkInfo.trkTimeErrorHyp[0] * trkInfo.trkTimeErrorHyp[0]);
        trkInfo.weight[1]=trkWeight / (trkInfo.trkTimeErrorHyp[1] * trkInfo.trkTimeErrorHyp[1]);
        trkInfo.weight[2]=trkWeight / (trkInfo.trkTimeErrorHyp[2] * trkInfo.trkTimeErrorHyp[2]);

        LOG << "vertexTimeFromTracks:     track"
            << " pt=" << trk.track().pt() << " eta=" << trk.track().eta() << " phi=" << trk.track().phi()
            << " vtxWeight=" << trkWeight << " time=" << trkTime << " timeError=" << trkTimeError
            << " timeQuality=" << trkTimeQuality << " timeHyp[pion]=" << trkInfo.trkTimeHyp[0]
            << " timeHyp[kaon]=" << trkInfo.trkTimeHyp[1] << " timeHyp[proton]=" << trkInfo.trkTimeHyp[2];
      }
    }
  }//End looping over tracks
  
  int Ntrks=vtx.originalTracks().size();

  //initiate a population of mass vectors 
  std::vector<int> population[populationSize_];
  for(int i=0;i<populationSize_;i++){
    population[i].reserve(Ntrks);
    for(int j=0;j<Ntrks;j++){
      int  //random int in [0,2]
      population[i].push_back();
    }
  }
  
  //check if the population meets all the requirement, if not regenerate


  //choose three mass vectors
  std::vector<int> parent=population[1];

  //make v_trackInfo with PID based on the parent vector 
  std::vector<std::vector<double>> v_trackInfowPID;
  for(int i=0;i<Ntrks;i++){
    v_trackInfowPID.emplace_back();
    auto&trkwpid = v_trackInfowPID.back();
    int tmp_pid=parent[N];
    double tmp_w=v_trackInfo[i].weight(tmp_pid);
    double tmp_t0=v_trackInfo[i].trkTimeHyp(tmp_pid);
    double tmp_sigma=v_trackInfo[i].trkTimeErrorHyp(tmp_pid);
    trkwpid={tmp_w,tmp_t0,tmp_sigma};     
  }

  double chi2_parent=calcChi2(v_trackInfowPID);

  std::vector<double> populationChi2;
  populationChi2.reserve(populationSize_);
  //create offspring generation 

  //natural (Darwinian) selection 

  return false;
}

double VertexTimeAlgorithmADEGA::solveVertexTime(std::vector<std::vector<double>>& trackvector){
  double sum=0.0;
  double wsum=0.0;
  for(const auto& trk: trackvector){
    double wgt=trk[0];
    double t0=trk[1];
    double sigmat0=trk[2];
    sum+=wgt*t0;
    wsum+=wgt;
  }
  double tv=sum/wsum;
  return tv;
}

double VertexTimeAlgorithmADEGA::calcChi2(std::vector<std::vector<double>>& trackvector,double tv){
  double chi2=0.0;
  for(const auto& trk: trackvector){
    double wgt=trk[0];
    double t0=trk[1];
    double sigmat0=trk[2];
    chi2+=wgt*(t0-tv)*(t0-tv);
  }
  
  return chi2;
}
