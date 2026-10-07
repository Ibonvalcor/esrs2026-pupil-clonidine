# Created by: Ibon Vales Cortina #
# Version: 2.1
# Validated by:
# Date: 2026-10-07
# Affiliation: Brain-Body-Regulation Lab, D-HEST, ETH Zürich
# Copyright (c) 2026, Ibon Vales Cortina. All rights reserved.
# Contact: ibon.valescortina@hest.ethz.ch

"""
Every number and figure of the poster, straight from the CSV files in data/.

    python code/make_figures.py

Takes a few seconds. You get 4 figures and 2 tables in results/, and the key numbers printed on screen.

The maths lives in the first four small functions; everything below them just draws.

v2.1: simpler maths (same numbers), header added; functions, figures and panels named as on the poster.
"""

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

DATA = Path(__file__).resolve().parents[1] / 'data'
RESULTS = DATA.parent / 'results'
DOSES = ['placebo', '0.075 mg', '0.150 mg']
COLORS = {'placebo': '#6F6F6F', '0.075 mg': '#215CAF', '0.150 mg': '#B7352D'}
STAGE_COLORS = {'W': '#B37C24', 'N1': '#60946E', 'N2': '#1B979B', 'N3': '#8483BE'}
plt.rcParams.update({'font.family': ['Arial', 'DejaVu Sans'], 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False})

effects = []            # every "dose vs placebo" comparison lands here, and then in results/effects.csv


# ======================================================================== the maths
def hedges_gz(d):
    """How big is a change? d = dose - placebo, one value per participant.
    Hedges' g_z = mean change / its standard deviation, times J, a correction for tiny samples (n = 4 here)."""
    df = len(d) - 1
    J = math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2))
    sd = d.std(ddof=1)
    low, high = bootstrap_ci(d)
    return J * d.mean() / sd, J * low / sd, J * high / sd


def bootstrap_ci(d, n_boot=10000, seed=42):
    """95 % confidence interval of the mean change (BCa bootstrap).
    Resample the participants 10,000 times and look at how the mean moves around. BCa then nudges the
    2.5 / 97.5 percentiles a little to correct for bias (z0) and skew (a) of those means."""
    rng = np.random.default_rng(seed)
    boot = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(axis=1)
    z0 = norm.ppf(np.mean(boot < d.mean()))
    leave_one_out = np.array([np.delete(d, i).mean() for i in range(len(d))])
    dev = leave_one_out.mean() - leave_one_out
    a = np.sum(dev ** 3) / (6 * np.sum(dev ** 2) ** 1.5)
    z = norm.ppf([0.025, 0.975])
    percentiles = norm.cdf(z0 + (z0 + z) / (1 - a * (z0 + z)))
    return np.percentile(boot, 100 * percentiles, method='hazen')


def smooth_circle(y):
    """Gaussian smoothing over 10 phase bins (sd = 2 bins) that wraps around: the pupil cycle is a circle,
    so the last bin sits right next to the first one."""
    shifts = np.arange(-5, 5)
    weights = np.exp(-0.5 * (shifts / 2) ** 2)
    weights /= weights.sum()
    return sum(w * np.roll(y, -k) for k, w in zip(shifts, weights))


def modulation_depth(curve):
    """How much does sigma (or the spindle chance) go up and down with the pupil? Highest minus lowest point."""
    smooth = smooth_circle(curve)
    return smooth.max() - smooth.min()


# ======================================================================== small helpers
def per_dose(df, column):
    """One row per participant, one column per dose."""
    return df.pivot(index='participant', columns='dose', values=column)[DOSES]


def compare_doses(part, measure, table):
    """Both doses vs placebo, for one measure."""
    for dose in DOSES[1:]:
        d = (table[dose] - table['placebo']).to_numpy()
        gz, low, high = hedges_gz(d)
        same_way = max(sum(d > 0), sum(d < 0))
        effects.append([part, measure, dose, gz, low, high, f'{same_way}/{len(d)}'])


def forest(ax, measures):
    """g_z with its 95 % CI for both doses, one line per measure."""
    for k, dose in enumerate(DOSES[1:]):
        for row, measure in enumerate(measures):
            gz, low, high = next(e[3:6] for e in effects if e[1] == measure and e[2] == dose)
            y = row + (k - 0.5) * 0.32
            ax.plot([low, high], [y, y], color=COLORS[dose])
            ax.plot(gz, y, 'o', color=COLORS[dose], label=f'{dose} vs placebo' if row == 0 else None)
    ax.axvline(0, color='grey', ls='--', lw=0.8)
    ax.set_yticks(range(len(measures)), measures)
    ax.set_ylim(len(measures) - 0.4, -0.6)
    ax.set_xlim(-2.2, 2.2)
    ax.set_xlabel("Hedges' $g_z$ (95 % CI)")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.0))


def save(fig, name):
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / f'{name}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


# ======================================================================== poster 1: pupil size tracks the brain
def pupil_and_exponent():
    """Every 30-s window of the 12 naps: pupil size against the aperiodic exponent of the EEG."""
    w = pd.read_csv(DATA / 'pupil_vs_exponent.csv')
    # repeated-measures correlation: remove each nap's own mean, then correlate what is left
    nap = w['participant'] + w['dose']
    x = w['aperiodic_exponent'] - w.groupby(nap)['aperiodic_exponent'].transform('mean')
    y = w['pupil_iris'] - w.groupby(nap)['pupil_iris'].transform('mean')
    r_rm = np.corrcoef(x, y)[0, 1]
    print(f'Pupil vs aperiodic exponent: r_rm = {r_rm:.2f}  ({len(w)} windows, 12 naps)')

    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    for stage, color in STAGE_COLORS.items():
        s = w[w['stage'] == stage]
        ax.scatter(s['aperiodic_exponent'], s['pupil_iris'], s=5, color=color, alpha=0.3, lw=0)
        # big dot: mean +- s.e.m. of the participants' medians (participants with at least 3 windows)
        counts = s.groupby('participant').size()
        med = s.groupby('participant')[['aperiodic_exponent', 'pupil_iris']].median()[counts >= 3]
        ax.errorbar(med['aperiodic_exponent'].mean(), med['pupil_iris'].mean(), xerr=med['aperiodic_exponent'].sem(),
                    yerr=med['pupil_iris'].sem(), fmt='o', ms=9, color=color, mec='white', label=stage)
    ax.set(xlabel='Aperiodic exponent (steeper = deeper sleep)', ylabel='Pupil / iris diameter',
           title=f'Pupil size tracks cortical activity ($r_{{rm}}$ = {r_rm:.2f})')
    ax.legend(frameon=False, ncol=4, fontsize=8)
    save(fig, '1_pupil_and_exponent')


# ======================================================================== poster 3a, 3b: same rhythms
RHYTHMS = {'pupil_isf_amplitude_pct': 'Pupil-ISF amplitude', 'pupil_isf_period_s': 'Pupil-ISF period',
           'sigma_isf_amplitude_pct': 'σ-ISF amplitude', 'sigma_isf_period_s': 'σ-ISF period',
           'spindle_density_per_h': 'Spindle density'}


def same_rhythms():
    """a: pupil size per sleep stage. b: the slow rhythms of pupil and sigma, and spindle density, vs placebo."""
    rhythms = pd.read_csv(DATA / 'rhythms.csv')
    for column, name in RHYTHMS.items():
        compare_doses('b rhythms', name, per_dose(rhythms, column))

    pupil = pd.read_csv(DATA / 'pupil_by_stage.csv')
    stages = ['W', 'N2', 'N3', 'REM']
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={'width_ratios': [1, 1.1], 'wspace': 0.55})
    for i, dose in enumerate(DOSES):
        table = pupil[pupil['dose'] == dose].pivot(index='participant', columns='stage', values='pupil_iris')
        table = table.reindex(columns=stages)
        x = np.arange(4) + (i - 1) * 0.13
        a.plot(x, table.T.to_numpy(), color=COLORS[dose], alpha=0.25)                  # each participant
        a.errorbar(x, table.mean(), yerr=table.sem(), fmt='-o', lw=2, capsize=3, color=COLORS[dose], label=dose)
    n = pupil.groupby(['stage', 'dose']).size()
    a.set_xticks(range(4), [f'{st}\nn = ' + '/'.join(str(n[st, dose]) for dose in DOSES) for st in stages])
    a.set(ylabel='Pupil / iris diameter', title='a  Pupil size per stage (faint: participants)')
    a.legend(frameon=False, fontsize=8)
    forest(b, list(RHYTHMS.values()))
    b.set_title('b  Rhythms vs placebo', pad=24)
    save(fig, '3ab_same_rhythms')


# ======================================================================== poster 3c, 3d, 3e: lost coupling
def lost_coupling():
    """c: one real pupil cycle, and sigma power and spindles along the pupil cycle (group curves).
    d, e: how much sigma power and spindles move with the pupil (modulation depth), per participant."""
    curves = pd.read_csv(DATA / 'phase_curves.csv').sort_values(['participant', 'dose', 'phase_deg'])
    cycle = pd.read_csv(DATA / 'example_cycle.csv')
    phase = np.sort(curves['phase_deg'].unique())
    fig = plt.figure(figsize=(8.5, 8))
    grid = fig.add_gridspec(3, 2, height_ratios=[0.75, 1, 1], width_ratios=[1.6, 1], hspace=0.55, wspace=0.35)

    c = fig.add_subplot(grid[0, :])
    c.plot(cycle['t_s'], cycle['pupil_iris'], '.', ms=2, color='#B5BCC6')
    c.plot(cycle['t_s'], cycle['pupil_iris_lowpass'], color='#2F3640', lw=2)
    for _, row in cycle.dropna(subset=['eye_frame']).iterrows():
        c.annotate(row['eye_frame'], (row['t_s'], row['pupil_iris_lowpass']), xytext=(0, 9),
                   textcoords='offset points', ha='center', fontsize=8)
    c.set(xlabel='Time (s)', ylabel='Pupil / iris', title='c  One real pupil cycle (P1, placebo) and the group curves')

    depths = []
    for row, (column, name, unit, panel, title) in enumerate([
            ('sigma_dpct', 'σ power', 'Δ% of mean', 'd', 'Sigma coupling per participant'),
            ('spindle_pct', 'Spindles', '%', 'e', 'Spindle coupling per participant')], start=1):
        left, right = fig.add_subplot(grid[row, 0]), fig.add_subplot(grid[row, 1])
        depth = pd.DataFrame(index=sorted(curves['participant'].unique()), columns=DOSES, dtype=float)
        for dose in DOSES:
            smooth = []
            for person, g in curves[curves['dose'] == dose].groupby('participant'):
                smooth.append(smooth_circle(g[column].to_numpy()))
                depth.loc[person, dose] = modulation_depth(g[column].to_numpy())
            smooth = np.array(smooth)
            mean, sem = smooth.mean(axis=0), smooth.std(axis=0, ddof=1) / np.sqrt(len(smooth))
            left.fill_between(phase, mean - sem, mean + sem, color=COLORS[dose], alpha=0.15, lw=0)
            left.plot(phase, mean, color=COLORS[dose], lw=2, label=dose)
        left.set_xticks([-180, -90, 0, 90, 180], ['trough', 'rise', 'peak', 'fall', 'trough'])
        left.set(ylabel=f'{name} ({unit})', title=f'{name} along the pupil cycle')

        compare_doses('c-e coupling', f'{name} coupling (modulation depth)', depth)
        gz, low, high, same_way = effects[-1][3:]
        change = 100 * (depth['0.150 mg'].mean() / depth['placebo'].mean() - 1)
        print(f'{name} coupling: {depth["placebo"].mean():.1f} -> {depth["0.150 mg"].mean():.1f} at 0.150 mg '
              f'({change:+.0f} %, lower in {same_way}), g_z = {gz:.2f} [{low:.2f}, {high:.2f}]')
        depths.append(depth.assign(measure=name))

        right.plot(range(3), depth.T.to_numpy(), color='#C9CED6', zorder=1)               # each participant
        for i, dose in enumerate(DOSES):
            right.scatter([i] * len(depth), depth[dose], s=18, color=COLORS[dose], zorder=3)
            right.hlines(depth[dose].mean(), i - 0.28, i + 0.28, color=COLORS[dose], lw=3)   # the mean
        right.set_xticks(range(3), DOSES)
        right.set(ylabel='Modulation depth', title=f'{panel}  {title}\n{change:+.0f} % at 0.150 mg (lower in {same_way})')

    fig.axes[1].legend(frameon=False, fontsize=8)
    fig.axes[-2].set_xlabel('Pupil-ISF phase')
    save(fig, '3cde_lost_coupling')
    pd.concat(depths).rename_axis('participant').reset_index().round(2).to_csv(RESULTS / 'modulation_depth.csv',
                                                                                index=False)


# ======================================================================== poster 4: sleep stayed largely comparable
SLEEP = {'total_sleep_time_min': 'Total sleep time', 'sleep_efficiency_pct': 'Sleep efficiency',
         'sleep_onset_latency_min': 'Sleep-onset latency', 'wake_pct': 'Wake (%)', 'n1_pct': 'N1 (%)',
         'n2_pct': 'N2 (%)', 'n3_pct': 'N3 (%)', 'rem_pct': 'REM (%)', 'so_density_per_h': 'Slow-oscillation density',
         'so_peak_to_peak_uv': 'Slow-oscillation amplitude', 'spindle_amplitude_uv': 'Spindle amplitude',
         'spindles_in_trains_pct': 'Spindles in trains'}


def sleep():
    """Sleep architecture, slow oscillations and spindles vs placebo."""
    naps = pd.read_csv(DATA / 'sleep.csv')
    for column, name in SLEEP.items():
        compare_doses('4 sleep', name, per_dose(naps, column))
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    forest(ax, list(SLEEP.values()))
    ax.set_title('Sleep vs placebo', pad=24)
    save(fig, '4_sleep')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')            # so σ prints fine on any console
    pupil_and_exponent()
    same_rhythms()
    lost_coupling()
    sleep()
    table = pd.DataFrame(effects, columns=['part', 'measure', 'dose', 'g_z', 'ci_low', 'ci_high', 'same_direction'])
    table.round(3).to_csv(RESULTS / 'effects.csv', index=False)
    print(f'Done! {len(table)} effect sizes in results/effects.csv, figures in results/')
