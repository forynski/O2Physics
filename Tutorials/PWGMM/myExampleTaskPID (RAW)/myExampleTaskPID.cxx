// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Enhanced PID analysis task using TPC and TOF PID information.
/// \author Robert Forynski, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTaskPid {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histo"};
  Configurable<int> nBinsdEdx{"nBinsdEdx", 500, "N bins in dE/dx histo"};
  Configurable<float> etaMin{"etaMin", -1.5f, "Minimum #eta acceptance"};
  Configurable<float> etaMax{"etaMax", 1.5f, "Maximum #eta acceptance"};
  Configurable<float> pidSigmaCutTPC{"pidSigmaCutTPC", 3.0f, "n#sigma cut for TPC PID"};
  Configurable<float> pidSigmaCutTOF{"pidSigmaCutTOF", 3.0f, "n#sigma cut for TOF PID"};

  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;

  // This is an example of a convenient declaration of "using"
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{30, etaMin, etaMax, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisdEdx{nBinsdEdx, 0, 300, "TPC signal"};
    const AxisSpec axisNSigma{200, -10, 10, "n#sigma"};
    const AxisSpec axisP{200, 0.1, 10, "p (GeV/c)"}; // Logarithmic scale would be better but linear for now

    const AxisSpec axisBeta{100, 0, 1.2, "#beta"};
    const AxisSpec axisMass{100, -0.2, 2.0, "Mass (GeV/c^{2})"};
    const AxisSpec axisCharge{3, -1.5, 1.5, "Charge"};
    const AxisSpec axisTheta{100, 0, TMath::Pi(), "#theta (rad)"};

    histos.add("betaHistogram", "Velocity #beta distribution", kTH1F, {axisBeta});
    histos.add("massHistogram", "Mass distribution (TOF)", kTH1F, {axisMass});
    histos.add("chargeHistogram", "Charge distribution", kTH1F, {axisCharge});
    histos.add("thetaHistogram", "Polar angle #theta", kTH1F, {axisTheta});

    // All-track histograms (before selection)
    histos.add("ptHistogramAllTracks", "p_{T} of all tracks", kTH1F, {axisPt});
    histos.add("etaHistogramAllTracks", "#eta of all tracks", kTH1F, {axisEta});
    histos.add("chargeHistogramAllTracks", "Charge of all tracks", kTH1F, {axisCharge});

    histos.add("nTracksPerCollision", "Total number of tracks per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nPionsPerCollision", "Number of identified pions per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nKaonsPerCollision", "Number of identified kaons per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nProtonsPerCollision", "Number of identified protons per collision", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("nElectronsPerCollision", "Number of identified electrons per collision", kTH1F, {{100, -0.5f, 99.5f}});


    // TODO: Bethe-Bloch initialisation to add particle identification lines to the dE/dx plot
    // constexpr float mPion = 0.13957f; // GeV/c²
    // auto hBB = histos.get<TH2>(HIST("tpcBBPi"));
    // for (double p = 0.1; p < 10; p += 0.05) {
    //   double bg = p/mPion; // Beta-gamma
    //   double bb = 70 * TMath::BetheBlochAleph(bg); 
    //   hBB->Fill(p, bb);
    // }

    // Create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt});

    // New histograms
    histos.add("tpcDedxVsPt", "TPC dE/dx vs p_{T};p_{T} (GeV/c);dE/dx (a.u.)", kTH2F, {axisPt, axisdEdx});
    histos.add("nSigmaPiVsP", "n#sigma(#pi) vs p;p (GeV/c);n#sigma(#pi)", kTH2F, {axisP, axisNSigma});  // New histogram
    
    // Histograms for TPC and TOF nSigma
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

    // TODO: Apply Bayesian approaches or detector combinations (TPC+TOF) for clean PID
    histos.add("ptSpectrumPion", "p_{T} spectrum for pions;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaon", "p_{T} spectrum for kaons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProton", "p_{T} spectrum for protons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectron", "p_{T} spectrum for electrons;p_{T} (GeV/c);Counts", kTH1F, {axisPt});

    // Combined PID histograms (NEW)
    histos.add("ptSpectrumPionCombinedPID", "Combined PID (TPC+TOF) #pi;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaonCombinedPID", "Combined PID (TPC+TOF) K;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProtonCombinedPID", "Combined PID (TPC+TOF) p;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectronCombinedPID", "Combined PID (TPC+TOF) e;p_{T} (GeV/c);Counts", kTH1F, {axisPt});


    histos.add("ptHistogramPion", "ptHistogramPion", kTH1F, {axisPt});
    histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1F, {axisPt});
    histos.add("ptHistogramProton", "ptHistogramProton", kTH1F, {axisPt});
    histos.add("ptHistogramElectron", "ptHistogramElectron", kTH1F, {axisPt});
    histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt});
    histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt});
    histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});
    histos.add("ptGeneratedElectron", "ptGeneratedElectron", kTH1F, {axisPt});

    histos.add("numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{10, -0.5f, 9.5f}});
    histos.add("multiplicityCorrelation", "multiplicityCorrelations", kTH2F, {{100, -0.5f, 99.5f}, {100, -0.5f, 99.5f}});

    // New histograms for dE/dx for reconstructed particles
    histos.add("dEdxPion", "dE/dx for Pions", kTH2F, {axisPt, axisdEdx});
    histos.add("dEdxKaon", "dE/dx for Kaons", kTH2F, {axisPt, axisdEdx});
    histos.add("dEdxProton", "dE/dx for Protons", kTH2F, {axisPt, axisdEdx});
    histos.add("dEdxElectron", "dE/dx for Electrons", kTH2F, {axisPt, axisdEdx});

    histos.add("ptSpectrumPionBayes", "Bayesian PID #pi;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumKaonBayes", "Bayesian PID K;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumProtonBayes", "Bayesian PID p;p_{T} (GeV/c);Counts", kTH1F, {axisPt});
    histos.add("ptSpectrumElectronBayes", "Bayesian PID e;p_{T} (GeV/c);Counts", kTH1F, {axisPt});


  }

  // Compute Bayesian probabilities for each species
  void computeBayesianPID(float nsTPC[4], float nsTOF[4], float priors[4], float probs[4]) {
    float sum = 0.f;
    for (int i = 0; i < 4; ++i) {
      // For tracks without TOF, set nsTOF[i]=0 or ignore that term
      float nsTOF2 = std::isfinite(nsTOF[i]) ? nsTOF[i]*nsTOF[i] : 0.f;
      float likelihood = std::exp(-0.5f * (nsTPC[i]*nsTPC[i] + nsTOF2));
      probs[i] = likelihood * priors[i];
      sum += probs[i];
    }
    // Normalize
    for (int i = 0; i < 4; ++i) {
      probs[i] = (sum > 0) ? probs[i]/sum : 0.f;
    }
  }


  // Process reconstructed tracks with TPC PID information
  // Fills histograms for nSigma distributions, TPC dE/dx vs pT, and PID-selected pT spectra for pions, kaons, protons, and electrons
  // Applies a 3σ cut for particle identification
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

    // Fill histograms for *all* tracks before applying cuts
    histos.fill(HIST("ptHistogramAllTracks"), pt);
    histos.fill(HIST("etaHistogramAllTracks"), eta);
    histos.fill(HIST("chargeHistogramAllTracks"), charge);

    // Now apply eta cut
    if (eta < etaMin || eta > etaMax) continue;

    float beta = track.beta();  // from pidTOFbeta
    float mass = track.mass();  // from pidTOFmass

    // Get TPC dE/dx signal and nSigma values for each particle type
    float tpcSignal = track.tpcSignal();

    float nSigmaPi = track.tpcNSigmaPi();  // Pion NSigma
    float nSigmaKa = track.tpcNSigmaKa();  // Kaon NSigma
    float nSigmaPr = track.tpcNSigmaPr();  // Proton NSigma
    float nSigmaEl = track.tpcNSigmaEl();  // Electron NSigma

    // Get TPC dE/dx signal and nSigma values for each particle type
    float nsTPCPi = track.tpcNSigmaPi();  // TPC nSigma for Pion
    float nsTPCKa = track.tpcNSigmaKa();  // TPC nSigma for Kaon
    float nsTPCPr = track.tpcNSigmaPr();  // TPC nSigma for Proton 
    float nsTPCEl = track.tpcNSigmaEl();  // TPC nSigma for Electron

    float nsTOFPi = track.tofNSigmaPi();  // TOF nSigma for Pion
    float nsTOFKa = track.tofNSigmaKa();  // TOF nSigma for Kaon
    float nsTOFPr = track.tofNSigmaPr();  // TOF nSigma for Proton
    float nsTOFEl = track.tofNSigmaEl();  // TOF nSigma for Electron

    float p = track.p(); // Changed from pt() to p()

    float theta = 2.0f * atanf(expf(-eta));  // You missed theta declaration
    histos.fill(HIST("thetaHistogram"), theta);
    histos.fill(HIST("betaHistogram"), beta);
    histos.fill(HIST("chargeHistogram"), charge);
    if (beta > 0.3f && mass > 0.0f) {
      histos.fill(HIST("massHistogram"), mass);
    }

    float nsTPC[4] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(), track.tpcNSigmaEl()};
    float nsTOF[4] = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(), track.tofNSigmaEl()};

    // Example flat priors (π, K, p, e)
    float priors[4] = {1.0f, 0.2f, 0.1f, 0.05f};
    float probs[4] = {0.f};

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

    // Example: fill histogram if max probability above threshold (e.g., >0.5)
    // float pt = track.pt();
    if (maxProb > 0.5f) {
      switch (maxSpecies) {
        case 0: histos.fill(HIST("ptSpectrumPionBayes"), pt); break;
        case 1: histos.fill(HIST("ptSpectrumKaonBayes"), pt); break;
        case 2: histos.fill(HIST("ptSpectrumProtonBayes"), pt); break;
        case 3: histos.fill(HIST("ptSpectrumElectronBayes"), pt); break;
      }
    }

    // Fill histograms for TPC dE/dx vs pT
    histos.fill(HIST("tpcDedxVsPt"), pt, tpcSignal);

    // Fill the new histogram
    histos.fill(HIST("nSigmaPiVsP"), p, nSigmaPi);

    // Fill histograms for TPC and TOF nSigma values
    histos.fill(HIST("nSigmaPiTPC"), pt, nsTPCPi);
    histos.fill(HIST("nSigmaKaTPC"), pt, nsTPCKa);
    histos.fill(HIST("nSigmaPrTPC"), pt, nsTPCPr);
    histos.fill(HIST("nSigmaElTPC"), pt, nsTPCEl);

    histos.fill(HIST("nSigmaPiTOF"), pt, nsTOFPi);
    histos.fill(HIST("nSigmaKaTOF"), pt, nsTOFKa);
    histos.fill(HIST("nSigmaPrTOF"), pt, nsTOFPr);
    histos.fill(HIST("nSigmaElTOF"), pt, nsTOFEl);

    histos.fill(HIST("nSigmaPiTPCvsTOF"), nsTPCPi, nsTOFPi);
    histos.fill(HIST("nSigmaKaTPCvsTOF"), nsTPCKa, nsTOFKa);
    histos.fill(HIST("nSigmaPrTPCvsTOF"), nsTPCPr, nsTOFPr);
    histos.fill(HIST("nSigmaElTPCvsTOF"), nsTPCEl, nsTOFEl);

    if (std::abs(nSigmaPi) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumPion"), pt);
    if (std::abs(nSigmaKa) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumKaon"), pt);
    if (std::abs(nSigmaPr) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumProton"), pt);
    if (std::abs(nSigmaEl) < pidSigmaCutTPC) histos.fill(HIST("ptSpectrumElectron"), pt);

    if (std::abs(nsTPCPi) < pidSigmaCutTPC && std::abs(nsTOFPi) < pidSigmaCutTOF) histos.fill(HIST("ptSpectrumPionCombinedPID"), pt);
    if (std::abs(nsTPCKa) < pidSigmaCutTPC && std::abs(nsTOFKa) < pidSigmaCutTOF) histos.fill(HIST("ptSpectrumKaonCombinedPID"), pt);
    if (std::abs(nsTPCPr) < pidSigmaCutTPC && std::abs(nsTOFPr) < pidSigmaCutTOF) histos.fill(HIST("ptSpectrumProtonCombinedPID"), pt);
    if (std::abs(nsTPCEl) < pidSigmaCutTPC && std::abs(nsTOFEl) < pidSigmaCutTOF) histos.fill(HIST("ptSpectrumElectronCombinedPID"), pt);

  
      // // Fill pT spectra for pions, kaons, protons and electrons using a 3-sigma PID cut
      // if (std::abs(nSigmaPi) < 3) {
      //   histos.fill(HIST("ptSpectrumPion"), pt);
      // }
      // if (std::abs(nSigmaKa) < 3) {
      //   histos.fill(HIST("ptSpectrumKaon"), pt);
      // }
      // if (std::abs(nSigmaPr) < 3) {
      //   histos.fill(HIST("ptSpectrumProton"), pt);
      // }
      // if (std::abs(nSigmaEl) < 3) {
      //   histos.fill(HIST("ptSpectrumElectron"), pt);
      // }
      // // Combined PID cut: |nSigma_TPC| < 3 AND |nSigma_TOF| < 3
      // if (std::abs(nsTPCPi) < 3 && std::abs(nsTOFPi) < 3) {
      //   histos.fill(HIST("ptSpectrumPionCombinedPID"), pt);
      // }
      // if (std::abs(nsTPCKa) < 3 && std::abs(nsTOFKa) < 3) {
      //   histos.fill(HIST("ptSpectrumKaonCombinedPID"), pt);
      // }
      // if (std::abs(nsTPCPr) < 3 && std::abs(nsTOFPr) < 3) {
      //   histos.fill(HIST("ptSpectrumProtonCombinedPID"), pt);
      // }
      // if (std::abs(nsTPCEl) < 3 && std::abs(nsTOFEl) < 3) {
      //   histos.fill(HIST("ptSpectrumElectronCombinedPID"), pt);
      // }
    }
  }
  PROCESS_SWITCH(myExampleTaskPid, process, "Process PID information", true);

  // Process reconstructed tracks matched to Monte Carlo truth
  // Fills histograms for reconstructed pT, eta, resolution (ΔpT), and dE/dx for specific particles (pions, kaons, protons, and electrons)
  // Uses MC truth to identify physical primaries within |y| < 0.5
  void processReco(aod::Collision const& collision, myFilteredTracks const& tracks, aod::McParticles const&)
{
  histos.fill(HIST("eventCounter"), 0.5);

  // Initialize per-collision counters
  int nTracks = 0;
  int nPions = 0;
  int nKaons = 0;
  int nProtons = 0;
  int nElectrons = 0;

  for (const auto& track : tracks) {
    if (track.tpcNClsCrossedRows() < 70) continue; // badly tracked
    ++nTracks;

    histos.fill(HIST("etaHistogram"), track.eta());
    histos.fill(HIST("ptHistogram"), track.pt());

    // Check if the track has a corresponding MC truth particle
    if (track.has_mcParticle()) { 
      auto mcParticle = track.mcParticle(); 
      histos.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());

      // Only process if the MC truth particle is a physical primary and has a |y| < 0.5
      if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) < 0.5) { 
        int pdg = abs(mcParticle.pdgCode());

        // Pion
        if (pdg == 211) {
          ++nPions;
          histos.fill(HIST("ptHistogramPion"), mcParticle.pt());
          histos.fill(HIST("dEdxPion"), track.pt(), track.tpcSignal());
        }
        // Kaon
        else if (pdg == 321) {
          ++nKaons;
          histos.fill(HIST("ptHistogramKaon"), mcParticle.pt());
          histos.fill(HIST("dEdxKaon"), track.pt(), track.tpcSignal());
        }
        // Proton
        else if (pdg == 2212) {
          ++nProtons;
          histos.fill(HIST("ptHistogramProton"), mcParticle.pt());
          histos.fill(HIST("dEdxProton"), track.pt(), track.tpcSignal());
        }
        // Electron
        else if (pdg == 11) {
          ++nElectrons;
          histos.fill(HIST("ptHistogramElectron"), mcParticle.pt());
          histos.fill(HIST("dEdxElectron"), track.pt(), track.tpcSignal());
        }
      }
    }
  }

  // Fill per-collision particle counts
  histos.fill(HIST("nTracksPerCollision"), nTracks);
  histos.fill(HIST("nPionsPerCollision"), nPions);
  histos.fill(HIST("nKaonsPerCollision"), nKaons);
  histos.fill(HIST("nProtonsPerCollision"), nProtons);
  histos.fill(HIST("nElectronsPerCollision"), nElectrons);
}

  PROCESS_SWITCH(myExampleTaskPid, processReco, "process reconstructed information", true);

  // Process Monte Carlo truth particles
  // Fills histograms for generated pT spectra of pions, kaons, protons, and electrons
  // Includes physical primary particles within |y| < 0.5 from the event generator
  void processSim(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, aod::McParticles const& mcParticles, myFilteredTracks const& tracks)
  {
    histos.fill(HIST("numberOfRecoCollisions"), collisions.size());

    std::vector<int> numberOfTracks;
    for (auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      numberOfTracks.emplace_back(groupedTracks.size());
    }
    if (collisions.size() == 2) histos.fill(HIST("multiplicityCorrelation"), numberOfTracks[0], numberOfTracks[1]);

    // Process MC truth particles
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) < 0.5) { 
        // Pion
        if (abs(mcParticle.pdgCode()) == 211) histos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
        // Kaon
        if (abs(mcParticle.pdgCode()) == 321) histos.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
        // Proton
        if (abs(mcParticle.pdgCode()) == 2212) histos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
        // Electron
        if (abs(mcParticle.pdgCode()) == 11) histos.fill(HIST("ptGeneratedElectron"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(myExampleTaskPid, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTaskPid>(cfgc)};
}
