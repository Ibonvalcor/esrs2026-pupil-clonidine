# Methods, in one page

**Design.** 4 healthy adults; after an adaptation nap, three experimental 90-min naps (7.1 ± 1.8 days apart) with
placebo, 0.075 mg or 0.150 mg clonidine, double-blind, in randomised order.

**Recording.** Polysomnography (EEG Fz, Cz, Pz; EOG; chin EMG; mastoid reference), ECG, and infrared video of
one eye (taped open). Sleep stages were scored in 30-s epochs.

**Pupil.** DeepLabCut tracked pupil and iris markers in the video, and ellipses fitted to them gave the
**pupil / iris diameter ratio**, which does not depend on camera distance. Invalid frames, blinks and values
outside 0.14–0.78 were removed. Values were clipped to the 1st–95th percentile and gaps were interpolated.

**Pupil infraslow rhythm (ISF) and its phase.** The pupil ratio was low-pass filtered at 0.05 Hz (Butterworth,
order 2, zero-phase), and its 200-s moving mean was subtracted. The Hilbert transform then gave the phase:
0° = pupil peak, ±180° = trough.

**Sigma power.** Central EEG (Cz) was band-passed at 13–15 Hz (Butterworth, order 2, zero-phase). Its squared
Hilbert envelope was low-pass filtered at 0.03 Hz and expressed as Δ% of the nap's mean NREM sigma power.

**Spindles.** Spindles were detected with A7 (Lacourse et al., 2019; 11–16 Hz) and curated. "Spindles (%)" is the
share of NREM samples in a phase bin that lie inside a spindle.

**Coupling (modulation depth).** NREM samples were sorted into 60 bins of pupil-ISF phase (6° each), separately
for each nap. Each nap's curve was smoothed around the circle (Gaussian, 10 bins). The modulation depth is the
curve's maximum minus its minimum (`data/phase_curves.csv` holds the unsmoothed curves).

**Rhythms.** Infraslow cycles of the pupil and of sigma power were detected trough to trough. For each nap we
report their mean amplitude (%) and duration (s), and the spindle density.

**Aperiodic exponent.** FOOOF 1.1 (fixed mode, 1–45 Hz) on 30-s windows of the frontal EEG (Fz). Figure 1 pairs it with the
pupil ratio of the same window (median of the valid samples, ≥ 50 % valid).

**Sleep.** Total sleep time, sleep efficiency, sleep-onset latency (to the first N2 epoch), stage % of time in bed,
slow-oscillation density and amplitude, spindle amplitude, and the share of spindles in trains.

**Statistics.** Each dose was compared with placebo within participants (d = dose − placebo):

- **Effect size:** Hedges' g<sub>z</sub> = J · mean(d) / sd(d), with the exact small-sample factor J.
- **Confidence interval:** 95 % BCa bootstrap of mean(d), 10,000 resamples, seed 42.
- **Direction:** "k of n" counts the participants who changed in the same direction.
- **Correlation:** the repeated-measures correlation r<sub>rm</sub> (Bakdash & Marusich, 2017) is the correlation
  that remains after each nap's mean is removed.
