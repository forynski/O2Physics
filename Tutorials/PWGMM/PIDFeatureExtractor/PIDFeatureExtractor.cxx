// PIDFeatureExtractor.cxx
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PIDFeatureExtractor {
  std::unique_ptr<TFile> outputFile;
  std::unique_ptr<TTree> featureTree;
  std::ofstream         csvFile;

  int   event_id;  // Use int for simplicity if using counter
  int   track_id;
  float px, py, pz, pt, p, eta, phi, theta;
  int   charge, track_type;
  float tpc_signal;
  float tpc_nsigma_pi, tpc_nsigma_ka, tpc_nsigma_pr, tpc_nsigma_el;
  int   tpc_nclusters;
  float tpc_chi2;
  float tof_beta, tof_mass;
  float tof_nsigma_pi, tof_nsigma_ka, tof_nsigma_pr, tof_nsigma_el;
  float bayes_prob_pi, bayes_prob_ka, bayes_prob_pr, bayes_prob_el;
  int   mc_pdg;
  float mc_px, mc_py, mc_pz;
  bool  has_tpc, has_tof;
  float dca_xy, dca_z;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> outputPath{"outputPath", "pid_features", "Output file base"};
  Configurable<bool>        exportCSV{"exportCSV", true, "Export CSV"};
  Configurable<bool>        exportROOT{"exportROOT", true, "Export ROOT"};
  Configurable<float>       etaMin{"etaMin", -1.5f, "Minimum eta"};
  Configurable<float>       etaMax{"etaMax",  1.5f, "Maximum eta"};
  Configurable<float>       ptMin{"ptMin",   0.1f, "Minimum pT"};
  Configurable<float>       ptMax{"ptMax",  20.0f, "Maximum pT"};

  int eventCounter = 0;  // For event ID fallback

  void init(InitContext const&) {
    std::string base = outputPath.value;  // Use .value member, not function
    if (exportROOT) {
      outputFile = std::make_unique<TFile>((base + ".root").c_str(), "RECREATE");
      featureTree = std::make_unique<TTree>("pid_features","PID features");

      featureTree->Branch("event_id",&event_id);
      featureTree->Branch("track_id",&track_id);
      featureTree->Branch("px",&px);
      featureTree->Branch("py",&py);
      featureTree->Branch("pz",&pz);
      featureTree->Branch("pt",&pt);
      featureTree->Branch("p",&p);
      featureTree->Branch("eta",&eta);
      featureTree->Branch("phi",&phi);
      featureTree->Branch("theta",&theta);
      featureTree->Branch("charge",&charge);
      featureTree->Branch("track_type",&track_type);

      featureTree->Branch("tpc_signal",&tpc_signal);
      featureTree->Branch("tpc_nsigma_pi",&tpc_nsigma_pi);
      featureTree->Branch("tpc_nsigma_ka",&tpc_nsigma_ka);
      featureTree->Branch("tpc_nsigma_pr",&tpc_nsigma_pr);
      featureTree->Branch("tpc_nsigma_el",&tpc_nsigma_el);
      featureTree->Branch("tpc_nclusters",&tpc_nclusters);
      featureTree->Branch("tpc_chi2",&tpc_chi2);

      featureTree->Branch("tof_beta",&tof_beta);
      featureTree->Branch("tof_mass",&tof_mass);
      featureTree->Branch("tof_nsigma_pi",&tof_nsigma_pi);
      featureTree->Branch("tof_nsigma_ka",&tof_nsigma_ka);
      featureTree->Branch("tof_nsigma_pr",&tof_nsigma_pr);
      featureTree->Branch("tof_nsigma_el",&tof_nsigma_el);

      featureTree->Branch("bayes_prob_pi",&bayes_prob_pi);
      featureTree->Branch("bayes_prob_ka",&bayes_prob_ka);
      featureTree->Branch("bayes_prob_pr",&bayes_prob_pr);
      featureTree->Branch("bayes_prob_el",&bayes_prob_el);

      featureTree->Branch("mc_pdg",&mc_pdg);
      featureTree->Branch("mc_px",&mc_px);
      featureTree->Branch("mc_py",&mc_py);
      featureTree->Branch("mc_pz",&mc_pz);

      featureTree->Branch("has_tpc",&has_tpc);
      featureTree->Branch("has_tof",&has_tof);
      featureTree->Branch("dca_xy",&dca_xy);
      featureTree->Branch("dca_z",&dca_z);
    }
    if (exportCSV) {
      csvFile.open((base + ".csv").c_str());
      csvFile <<
        "event_id,track_id,px,py,pz,pt,p,eta,phi,theta,charge,track_type,"
        "tpc_signal,tpc_nsigma_pi,tpc_nsigma_ka,tpc_nsigma_pr,tpc_nsigma_el,"
        "tpc_nclusters,tpc_chi2,"
        "tof_beta,tof_mass,tof_nsigma_pi,tof_nsigma_ka,tof_nsigma_pr,tof_nsigma_el,"
        "bayes_prob_pi,bayes_prob_ka,bayes_prob_pr,bayes_prob_el,"
        "mc_pdg,mc_px,mc_py,mc_pz,has_tpc,has_tof,dca_xy,dca_z\n";
    }

    const AxisSpec axisPt{200, 0,10, "pT"};
    const AxisSpec axisEta{ 60,-1.5,1.5,"eta"};
    const AxisSpec axisdEdx{300,0,300,"dE/dx"};
    const AxisSpec axisBeta{120,0,1.2,"beta"};
    const AxisSpec axisMass{100,-0.2,2.0,"mass"};
    histos.add("QC/nTracks","Tracks",kTH1F,{{10000,0,100000}});
    histos.add("QC/pt","pT",kTH1F,{axisPt});
    histos.add("QC/eta","eta",kTH1F,{axisEta});
    histos.add("QC/tpc_dEdx_vs_pt","dE/dx vs pT",kTH2F,{axisPt,axisdEdx});
    histos.add("QC/tof_beta_vs_p","beta vs p",kTH2F,{axisPt,axisBeta});
    histos.add("QC/mass_vs_p","mass vs p",kTH2F,{axisPt,axisMass});
  }

  void computeBayesianPID(float nsTPC[4], float nsTOF[4], float pri[4], float out[4]) {
    float sum=0;
    for(int i=0;i<4;i++){
      float l=std::exp(-0.5f*(nsTPC[i]*nsTPC[i] + (std::isfinite(nsTOF[i])?nsTOF[i]*nsTOF[i]:0.f)));
      out[i]=l*pri[i];
      sum+=out[i];
    }
    for(int i=0;i<4;i++) out[i]=sum>0?out[i]/sum:0.f;
  }

  void process(
    aod::Collision const& collision,
    soa::Join<
      aod::Tracks,
      aod::TracksExtra,
      aod::TracksDCA,
      aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl,
      aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFEl,
      aod::pidTOFmass, aod::pidTOFbeta,
      aod::McTrackLabels
    > const& tracks,
    aod::McParticles const&)
  {
    // Use a simple event counter since collision timestamp access isn't available
    static int eventCounter=0;
    event_id = eventCounter++;
    int idx=0;
    for(auto& t:tracks) {
      if(t.pt()<ptMin||t.pt()>ptMax) continue;
      if(t.eta()<etaMin||t.eta()>etaMax) continue;
      track_id=idx++;
      px = t.px(); py = t.py(); pz = t.pz();
      pt = t.pt(); p = t.p();
      eta = t.eta(); phi = t.phi();
      theta = 2.f*atanf(expf(-eta));
      charge = t.sign();
      track_type = t.trackType();

      has_tpc = t.hasTPC();
      if(has_tpc) {
        tpc_signal = t.tpcSignal();
        tpc_nsigma_pi = t.tpcNSigmaPi();
        tpc_nsigma_ka = t.tpcNSigmaKa();
        tpc_nsigma_pr = t.tpcNSigmaPr();
        tpc_nsigma_el = t.tpcNSigmaEl();
        tpc_nclusters = t.tpcNClsFound();
        tpc_chi2 = t.tpcChi2NCl();
      } else {
        tpc_signal = tpc_nsigma_pi = tpc_nsigma_ka = tpc_nsigma_pr = tpc_nsigma_el = -999;
        tpc_nclusters = 0; tpc_chi2 = -999;
      }

      has_tof = t.hasTOF();
      if(has_tof) {
        tof_beta = t.beta();
        tof_mass = t.mass();
        tof_nsigma_pi = t.tofNSigmaPi();
        tof_nsigma_ka = t.tofNSigmaKa();
        tof_nsigma_pr = t.tofNSigmaPr();
        tof_nsigma_el = t.tofNSigmaEl();
      } else {
        tof_beta = tof_mass = -999;
        tof_nsigma_pi = tof_nsigma_ka = tof_nsigma_pr = tof_nsigma_el = -999;
      }

      dca_xy = t.dcaXY();
      dca_z = t.dcaZ();

      float arrTPC[4] = {tpc_nsigma_pi, tpc_nsigma_ka, tpc_nsigma_pr, tpc_nsigma_el};
      float arrTOF[4] = {tof_nsigma_pi, tof_nsigma_ka, tof_nsigma_pr, tof_nsigma_el};
      float priors[4] = {1.f, 0.2f, 0.1f, 0.05f};
      float probs[4];
      computeBayesianPID(arrTPC, arrTOF, priors, probs);
      bayes_prob_pi = probs[0];
      bayes_prob_ka = probs[1];
      bayes_prob_pr = probs[2];
      bayes_prob_el = probs[3];

      auto mc = t.mcParticle();
      mc_pdg = mc.pdgCode();
      mc_px = mc.px();
      mc_py = mc.py();
      mc_pz = mc.pz();

      if(exportROOT) featureTree->Fill();
      if(exportCSV) {
        csvFile << event_id << "," << track_id << ","
                << px << "," << py << "," << pz << ","
                << pt << "," << p << ","
                << eta << "," << phi << "," << theta << ","
                << charge << "," << track_type << ","
                << tpc_signal << "," << tpc_nsigma_pi << "," << tpc_nsigma_ka << "," << tpc_nsigma_pr << "," << tpc_nsigma_el << ","
                << tpc_nclusters << "," << tpc_chi2 << ","
                << tof_beta << "," << tof_mass << "," << tof_nsigma_pi << "," << tof_nsigma_ka << "," << tof_nsigma_pr << "," << tof_nsigma_el << ","
                << bayes_prob_pi << "," << bayes_prob_ka << "," << bayes_prob_pr << "," << bayes_prob_el << ","
                << mc_pdg << "," << mc_px << "," << mc_py << "," << mc_pz << ","
                << has_tpc << "," << has_tof << ","
                << dca_xy << "," << dca_z << "\n";
      }

      histos.fill(HIST("QC/nTracks"), 1);
      histos.fill(HIST("QC/pt"), pt);
      histos.fill(HIST("QC/eta"), eta);
      if(has_tpc) histos.fill(HIST("QC/tpc_dEdx_vs_pt"), pt, tpc_signal);
      if(has_tof) {
        histos.fill(HIST("QC/tof_beta_vs_p"), p, tof_beta);
        histos.fill(HIST("QC/mass_vs_p"), p, tof_mass);
      }
    }
  }

  void finalize(){
    if(exportROOT){
      outputFile->cd();
      featureTree->Write();
      outputFile->Close();
    }
    if(exportCSV){
      csvFile.close();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{ adaptAnalysisTask<PIDFeatureExtractor>(cfgc) };
}
