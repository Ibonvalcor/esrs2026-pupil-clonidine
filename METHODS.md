# The nerdy details 🤓

Everything about how the numbers were made, in one page.

**The study.** 4 healthy adults. First an adaptation nap to get used to the lab, then three 90-min naps, 7.1 ± 1.8
days apart. Before each one they took placebo, 0.075 mg or 0.150 mg clonidine, in random order. Double-blind:
neither they nor we knew which was which.

**What we recorded.**
- **Brain and body:** EEG at Fz, Cz and Pz, plus eye movements (EOG), chin muscle (EMG) and heart (ECG).
- **The eye:** infrared video of one eye, taped open.
- **Sleep stages:** scored in 30-s epochs.

**Pupil size.**
- **Tracking:** DeepLabCut found the pupil and the iris in every video frame, and we fitted ellipses to them.
- **The measure:** the pupil / iris diameter ratio. It doesn't care how close the camera is.
- **Cleaning:** we threw out blinks, invalid frames and values outside 0.14–0.78. Then we clipped the rest to the
  1st–95th percentile and filled the gaps by interpolation.

**The pupil's slow rhythm (ISF) and its phase.**
1. Keep only the slow stuff: a low-pass filter at 0.05 Hz (Butterworth, order 2, run forwards and backwards).
2. Remove the slow drift: subtract a 200-s moving average.
3. The Hilbert transform tells us where the pupil is in its cycle at every moment (0° = peak, ±180° = trough).

**Sigma power.**
1. Take the EEG at Cz and keep 13–15 Hz (Butterworth, order 2, zero-phase).
2. Its power (the squared Hilbert envelope), low-pass filtered at 0.03 Hz.
3. Expressed as Δ% of the nap's mean NREM sigma power.

**Spindles.** Detected with A7 (Lacourse et al., 2019; 11–16 Hz), then curated. "Spindles (%)" is the share
of NREM moments in a phase bin that fall inside a spindle.

**Coupling (modulation depth).**
1. Sort the NREM moments of each nap into 60 bins of pupil phase (6° each).
2. Average sigma power (or the spindle %) in each bin. The raw curves are in `data/phase_curves.csv`.
3. Smooth each nap's curve around the circle (Gaussian, 10 bins).
4. The modulation depth is the highest point minus the lowest point.

**Rhythms.** We found the slow cycles of the pupil and of sigma power (trough to trough). For each nap we give
their mean size (amplitude, %) and length (period, s), plus the spindle density.

**Aperiodic exponent.** FOOOF 1.1 (fixed mode, 1–45 Hz) on 30-s windows of the frontal EEG (Fz). Figure 1 pairs it
with the pupil ratio of the same window (the median of its valid samples, at least 50 % valid).

**Sleep.**
- **Architecture:** total sleep time, sleep efficiency, time until the first N2 epoch, and % of time in bed spent in
  each stage.
- **Slow oscillations:** density and amplitude.
- **Spindles:** amplitude, and how many come in trains.

**Stats.** Each dose is compared with placebo within the same person (d = dose − placebo).
- **Hedges' g<sub>z</sub>:** mean(d) / sd(d), times a small-sample correction J.
- **95 % CI:** BCa bootstrap of mean(d), 10,000 resamples, seed 42.
- **"Lower in 4/4":** how many people changed in the same direction.
- **r<sub>rm</sub>:** the repeated-measures correlation (Bakdash & Marusich, 2017), i.e. the correlation that is
  left after removing each nap's own mean.
