# Data Analysis Histograms Overview

This document provides a detailed explanation of the histograms used in the data analysis of **pp LHC23zzh pass4 (LHC23g3) MC anchored data**, including their classification into **Reconstructed Histograms** and **Simulated (MC Truth) Histograms**.

## Dataset Information
Analysis performed on:
- **AO2D file** from 2024 MC production
- Anchored to pp LHC23zzh pass4 (LHC23g3) dataset
- Uses both reconstructed tracks and MC truth information

## Reconstructed Histograms:
These histograms correspond to the **reconstructed tracks** from the detector, filled with track-level information:

- **eventCounter**: Number of processed events
- **etaHistogram**: η distribution of selected tracks
- **ptHistogram**: Transverse momentum ($p_T$) of selected tracks
- **ptResolution**: $p_T$ resolution (reconstructed vs generated)
- **tpcDedxVsPt**: TPC $dE/dx$ vs $p_T$ correlation
- **nSigmaPiVsPt**, **nSigmaKaVsPt**, **nSigmaPrVsPt**, **nSigmaElVsPt**: $n\sigma$ distributions for pions, kaons, protons, and electrons vs $p_T$
- **nSigmaPiVsP**: $n\sigma(\pi)$ vs $p$ distribution
- **ptSpectrumPion**, **ptSpectrumKaon**, **ptSpectrumProton**, **ptSpectrumElectron**: $p_T$ spectra using configurable $n\sigma$ PID cuts
- **Combined PID histograms** (TPC + TOF):
  - **ptSpectrumPionCombinedPID**
  - **ptSpectrumKaonCombinedPID**
  - **ptSpectrumProtonCombinedPID**
  - **ptSpectrumElectronCombinedPID**
- Particle-specific histograms:
  - **ptHistogramPion/Kaon/Proton/Electron**: Reconstructed $p_T$ for each species
  - **dEdxPion/Kaon/Proton/Electron**: $dE/dx$ vs $p_T$ for identified particles

### New Kinematic Histograms (Run 3 PID Upgrade):
- **thetaHistogram**: Polar angle ($\theta$) computed from pseudorapidity
- **betaHistogram**: Velocity $\beta$ from TOF
- **massHistogram**: Mass estimation using TOF beta and momentum
- **chargeHistogram**: Distribution of track charge

### All Tracks vs Selected Tracks Comparison:
These histograms capture the full population of tracks (before any selection cuts) and allow comparison with selected-track distributions:
- **ptHistogramAllTracks**: $p_T$ of all reconstructed tracks
- **etaHistogramAllTracks**: η of all reconstructed tracks
- **chargeHistogramAllTracks**: Charge distribution before cuts

### Per-Collision Particle Statistics:
New histograms added to inspect particle production per collision:
- **nTracksPerCollision**: Total number of selected tracks per collision
- **nPionsPerCollision**, **nKaonsPerCollision**, **nProtonsPerCollision**, **nElectronsPerCollision**: Number of matched physical primaries (|y| < 0.5) per reconstructed collision, for each of the four most common species

## Simulated (MC Truth) Histograms:
These histograms use **MC particle-level** information:

- **ptGeneratedPion/Kaon/Proton/Electron**: Generated $p_T$ of true MC particles
- **numberOfRecoCollisions**: Number of reconstructed collisions per MC collision
- **multiplicityCorrelation**: Track multiplicity correlation between reconstructed collisions

## Analysis Task Features:
The associated analysis task ([myExampleTaskPID.cxx](myExampleTaskPID.cxx)) includes:
1. Track selection with DCA cut (< 0.2 cm) and configurable η cut
2. PID validation using:
   - TPC signal ($dE/dx$) vs momentum
   - $n\sigma$ separation for pions, kaons, protons, electrons (TPC and TOF)
   - Configurable PID cuts for flexibility (via JSON)
3. Extended physics observables:
   - $\beta$ and mass via TOF-based PID
   - $\theta$ from $\eta$, charge histograms
4. MC truth matching for:
   - $p_T$ resolution studies
   - Physical primary identification ($|y| < 0.5$)
5. Collision-track association and multiplicity studies
6. All-track histograms for comparison with selected populations
7. Per-collision particle yields for physics performance validation

## Summary of Histogram Classification:

| Reconstruction Histograms               | MC Truth Histograms         |
|-----------------------------------------|-----------------------------|
| `eventCounter`                          | `ptGenerated*`              |
| `etaHistogram`, `thetaHistogram`        | `numberOfRecoCollisions`    |
| `ptHistogram*`, `ptHistogramAllTracks`  | `multiplicityCorrelation`   |
| `dEdx*`, `tpcDedxVsPt`                  |                             |
| `nSigma*VsPt`, `nSigma*VsP`             |                             |
| `betaHistogram`, `massHistogram`        |                             |
| `chargeHistogram`, `chargeHistogramAllTracks` |                       |
| `nTracksPerCollision`, `nPionsPerCollision`   |                       |
| `nKaonsPerCollision`, `nProtonsPerCollision`  |                       |
| `nElectronsPerCollision` |                                            |

## Key Notes:
- The analysis is structured using O2Physics PID helpers (`pidTPC`, `pidTOF`, `pidTOFmass`, `pidTOFbeta`)
- Supports both single-detector and combined (TPC + TOF) PID
- Mass estimation relies on TOF-based $\beta$ with momentum using $m = p\sqrt{1/\beta^2 - 1}$
- Histogram binning and selection cuts are configurable using `Configurable<>` settings in JSON
- Full-track vs selected-track histograms aid in assessing the effect of selection criteria
- Per-collision particle statistics allow monitoring particle production per event for performance studies