// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTaskPid {
  // Histogram registry
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable parameters (with explicit template arguments)
  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histo"};
  Configurable<int> nBinsdEdx{"nBinsdEdx", 500, "N bins in dE/dx histo"};
  Configurable<float> etaMin{"etaMin", -1.5f, "Minimum #eta acceptance"};
  Configurable<float> etaMax{"etaMax", 1.5f, "Maximum #eta acceptance"};
  Configurable<float> pidSigmaCutTPC{"pidSigmaCutTPC", 3.0f, "n#sigma cut for TPC PID"};
  Configurable<float> pidSigmaCutTOF{"pidSigmaCutTOF", 3.0f, "n#sigma cut for TOF PID"};

  void init(InitContext const&) {
    // Define axes
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{30, etaMin, etaMax, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisdEdx{nBinsdEdx, 0, 300, "TPC signal"};
    const AxisSpec axisNSigma{200, -10, 10, "n#sigma"};
    const AxisSpec axisP{200, 0.1, 10, "p (GeV/c)"};
    const AxisSpec axisBeta{100, 0, 1.2, "#beta"};
    const AxisSpec axisMass{100, -0.2, 2.0, "Mass (GeV/c^{2})"};
    const AxisSpec axisCharge{3, -1.5, 1.5, "Charge"};
    const AxisSpec axisTheta{100, 0, TMath::Pi(), "#theta (rad)"};

    // Histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt});
    histos.add("betaHistogram", "Velocity #beta distribution", kTH1F, {axisBeta});
    histos.add("massHistogram", "Mass distribution (TOF)", kTH1F, {axisMass});
    histos.add("chargeHistogram", "Charge distribution", kTH1F, {axisCharge});
    histos.add("thetaHistogram", "Polar angle #theta", kTH1F, {axisTheta});
    histos.add("ptHistogramAllTracks", "p_{T} of all tracks", kTH1F, {axisPt});
    histos.add("etaHistogramAllTracks", "#eta of all tracks", kTH1F, {axisEta});
    histos.add("chargeHistogramAllTracks", "Charge of all tracks", kTH1F, {axisCharge});
    histos.add("nTracksPerCollision", "Total number of tracks per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nPionsPerCollision", "Number of identified pions per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nKaonsPerCollision", "Number of identified kaons per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nProtonsPerCollision", "Number of identified protons per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nElectronsPerCollision", "Number of identified electrons per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("tpcDedxVsPt", "TPC dE/dx vs p_{T};p_{T} (GeV/c);dE/dx (a.u.)", kTH2F, {axisPt, axisdEdx});
    histos.add("nSigmaPiVsP", "n#sigma(#pi) vs p;p (GeV/c);n#sigma(#pi)", kTH2F, {axisP, axisNSigma});
    histos.add("nSigmaPiTPC", "n#sigma_{TPC}(#pi) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaKaTPC", "n#sigma_{TPC}(K) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaPrTPC", "n#sigma_{TPC}(p) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaElTPC", "n#sigma_{TPC}(e) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaPiTOF", "n#sigma_{TOF}(#pi) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaKaTOF", "n#sigma_{TOF}(K) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaPrTOF", "n#sigma_{TOF}(p) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaElTOF", "n#sigma_{TOF}(e) vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("nSigmaPiTPCvsTOF", "n#sigma_{TPC} vs n#sigma_{TOF} (#pi)", kTH2F, {axisNSigma, axisNSigma});
    histos.add("nSigmaKaTPCvsTOF", "n#sigma_{TPC} vs n#sigma_{TOF} (K)", kTH2F, {axisNSigma, axisNSigma});
    histos.add("nSigmaPrTPCvsTOF", "n#sigma_{TPC} vs n#sigma_{TOF} (p)", kTH2F, {axisNSigma, axisNSigma});
    histos.add("nSigmaElTPCvsTOF", "n#sigma_{TPC} vs n#sigma_{TOF} (e)", kTH2F, {axisNSigma, axisNSigma});
    histos.add("ptSpectrumPion", "p_{T} spectrum for pions;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaon", "p_{T} spectrum for kaons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProton", "p_{T} spectrum for protons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectron", "p_{T} spectrum for electrons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumPionCombinedPID", "Combined PID (TPC+TOF) #pi;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaonCombinedPID", "Combined PID (TPC+TOF) K;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProtonCombinedPID", "Combined PID (TPC+TOF) p;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectronCombinedPID", "Combined PID (TPC+TOF) e;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumPionBayes", "Bayesian PID #pi;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaonBayes", "Bayesian PID K;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProtonBayes", "Bayesian PID p;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectronBayes", "Bayesian PID e;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
  }

  // Compute Bayesian PID probabilities
  void computeBayesianPID(float nsTPC[4], float nsTOF[4], float priors[4], float probs[4]) {
    float sum = 0.f;
    for (int i = 0; i < 4; ++i) {
      float nsTOF2 = std::isfinite(nsTOF[i]) ? nsTOF[i]*nsTOF[i] : 0.f;
      float likelihood = std::exp(-0.5f * (nsTPC[i]*nsTPC[i] + nsTOF2));
      probs[i] = likelihood * priors[i];
      sum += probs[i];
    }
    for (int i = 0; i < 4; ++i) {
      probs[i] = (sum > 0) ? probs[i]/sum : 0.f;
    }
  }

  // Main process function for reconstructed tracks
  void process(soa::Join<
    aod::Tracks,
    aod::TracksExtra,
    aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl,
    aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFEl,
    aod::pidTOFmass, aod::pidTOFbeta> const& tracks)
  {
    for (auto& track : tracks) {
      float pt = track.pt();
      float eta = track.eta();
      int charge = track.sign();

      // Fill histograms for all tracks before cuts
      histos.fill(HIST("ptHistogramAllTracks"), pt);
      histos.fill(HIST("etaHistogramAllTracks"), eta);
      histos.fill(HIST("chargeHistogramAllTracks"), charge);

      // Apply eta cut
      if (eta < etaMin || eta > etaMax) continue;

      float beta = track.beta();
      float mass = track.mass();
      float tpcSignal = track.tpcSignal();

      float nSigmaPi = track.tpcNSigmaPi();
      float nSigmaKa = track.tpcNSigmaKa();
      float nSigmaPr = track.tpcNSigmaPr();
      float nSigmaEl = track.tpcNSigmaEl();

      float nsTPC[4] = {nSigmaPi, nSigmaKa, nSigmaPr, nSigmaEl};
      float nsTOF[4] = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(), track.tofNSigmaEl()};
      float priors[4] = {1.0f, 0.2f, 0.1f, 0.05f};
      float probs[4] = {0.f};

      float p = track.p();
      float theta = 2.0f * atanf(expf(-eta));

      histos.fill(HIST("thetaHistogram"), theta);
      histos.fill(HIST("betaHistogram"), beta);
      histos.fill(HIST("chargeHistogram"), charge);

      if (beta > 0.3f && mass > 0.0f) {
        histos.fill(HIST("massHistogram"), mass);
      }

      computeBayesianPID(nsTPC, nsTOF, priors, probs);

      // Find species with maximum probability
      int maxSpecies = 0;
      float maxProb = probs[0];
      for (int i = 1; i < 4; ++i) {
        if (probs[i] > maxProb) {
          maxProb = probs[i];
          maxSpecies = i;
        }
      }
      if (maxProb > 0.5f) {
        switch (maxSpecies) {
          case 0: histos.fill(HIST("ptSpectrumPionBayes"), pt); break;
          case 1: histos.fill(HIST("ptSpectrumKaonBayes"), pt); break;
          case 2: histos.fill(HIST("ptSpectrumProtonBayes"), pt); break;
          case 3: histos.fill(HIST("ptSpectrumElectronBayes"), pt); break;
        }
      }

      // Fill TPC dE/dx vs pT
      histos.fill(HIST("tpcDedxVsPt"), pt, tpcSignal);
      histos.fill(HIST("nSigmaPiVsP"), p, nSigmaPi);

      // Fill histograms for TPC and TOF nSigma values
      histos.fill(HIST("nSigmaPiTPC"), pt, nSigmaPi);
      histos.fill(HIST("nSigmaKaTPC"), pt, nSigmaKa);
      histos.fill(HIST("nSigmaPrTPC"), pt, nSigmaPr);
      histos.fill(HIST("nSigmaElTPC"), pt, nSigmaEl);
      histos.fill(HIST("nSigmaPiTOF"), pt, nsTOF[0]);
      histos.fill(HIST("nSigmaKaTOF"), pt, nsTOF[1]);
      histos.fill(HIST("nSigmaPrTOF"), pt, nsTOF[2]);
      histos.fill(HIST("nSigmaElTOF"), pt, nsTOF[3]);
      histos.fill(HIST("nSigmaPiTPCvsTOF"), nSigmaPi, nsTOF[0]);
      histos.fill(HIST("nSigmaKaTPCvsTOF"), nSigmaKa, nsTOF[1]);
      histos.fill(HIST("nSigmaPrTPCvsTOF"), nSigmaPr, nsTOF[2]);
      histos.fill(HIST("nSigmaElTPCvsTOF"), nSigmaEl, nsTOF[3]);

      // PID-selected spectra using TPC only
      if (std::abs(nSigmaPi) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumPion"), pt);
      if (std::abs(nSigmaKa) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumKaon"), pt);
      if (std::abs(nSigmaPr) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumProton"), pt);
      if (std::abs(nSigmaEl) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumElectron"), pt);

      // Combined PID (TPC+TOF)
      if (std::abs(nSigmaPi) < pidSigmaCutTPC && std::abs(nsTOF[0]) < pidSigmaCutTOF)
        histos.fill(HIST("ptSpectrumPionCombinedPID"), pt);
      if (std::abs(nSigmaKa) < pidSigmaCutTPC && std::abs(nsTOF[1]) < pidSigmaCutTOF)
        histos.fill(HIST("ptSpectrumKaonCombinedPID"), pt);
      if (std::abs(nSigmaPr) < pidSigmaCutTPC && std::abs(nsTOF[2]) < pidSigmaCutTOF)
        histos.fill(HIST("ptSpectrumProtonCombinedPID"), pt);
      if (std::abs(nSigmaEl) < pidSigmaCutTPC && std::abs(nsTOF[3]) < pidSigmaCutTOF)
        histos.fill(HIST("ptSpectrumElectronCombinedPID"), pt);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTaskPid>(cfgc)
  };
}
