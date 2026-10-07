% Created by: Ibon Vales Cortina %
% Version: 2.1
% Validated by:
% Date: 2026-10-07
% Affiliation: Brain-Body-Regulation Lab, D-HEST, ETH Zürich
% Copyright (c) 2026, Ibon Vales Cortina. All rights reserved.
% Contact: ibon.valescortina@hest.ethz.ch

"""
make_figures.py - the numbers and figures of the poster, from the CSV files in data/.

    python code/make_figures.py

Writes four figures, results/effects.csv (every effect size) and results/modulation_depth.csv, and prints the key
numbers of the poster. One function per part of the poster; read it from top to bottom.

Ibon Vales Cortina, Brain-Body Regulation Lab, ETH Zurich (2026). MIT licence.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import erf, erfinv, gammaln

ROOT = Path(__file__).resolve().parents[1]
DATA, RESULTS = ROOT / 'data', ROOT / 'results'
DOSES = ['placebo', '0.075 mg', '0.150 mg']
DOSE_COLOR = {'placebo': '#6F6F6F', '0.075 mg': '#215CAF', '0.150 mg': '#B7352D'}
STAGE_COLOR = {'W': '#B37C24', 'N1': '#60946E', 'N2': '#1B979B', 'N3': '#8483BE', 'REM': '#B0779C'}
plt.rcParams.update({'font.family': ['Arial', 'DejaVu Sans'], 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False})
EFFECTS = []                                       # one row per effect: part, measure, dose, g_z, CI, same direction


# ---------------------------------------------------------------- statistics
def hedges_gz(placebo, dose, n_boot=10000, seed=42):
    """Paired effect of a dose against placebo. d = dose - placebo per participant; Hedges' g_z = J * mean(d) / sd(d)
    with the exact small-sample correction J, its 95 % CI from a BCa bootstrap of mean(d), and how many
    participants changed in the same direction."""
    d = np.asarray(dose, float) - np.asarray(placebo, float)
    d = d[np.isfinite(d)]
    n, sd = d.size, d.std(ddof=1)
    j = np.exp(gammaln((n - 1) / 2) - 0.5 * np.log((n - 1) / 2) - gammaln((n - 2) / 2))
    lo, hi = bca_ci(d, n_boot, np.random.default_rng(seed))
    same = max(np.sum(d > 0), np.sum(d < 0))
    return j * d.mean() / sd, j * lo / sd, j * hi / sd, f'{same}/{n}'


def bca_ci(x, n_boot, rng, alpha=0.05):
    """Bias-corrected and accelerated (BCa) bootstrap 95 % CI of the mean of x."""
    boot = x[rng.integers(0, x.size, size=(n_boot, x.size))].mean(axis=1)
    pct = lambda q: np.percentile(boot, 100 * q, method='hazen')                    # noqa: E731
    p0 = np.mean(boot < x.mean())
    if p0 in (0, 1):
        return pct(alpha / 2), pct(1 - alpha / 2)
    z = lambda p: np.sqrt(2) * erfinv(2 * p - 1)                                      # noqa: E731
    phi = lambda v: 0.5 * (1 + erf(v / np.sqrt(2)))                                   # noqa: E731
    jack = np.array([np.delete(x, i).mean() for i in range(x.size)])
    a = np.sum((jack.mean() - jack) ** 3) / (6 * np.sum((jack.mean() - jack) ** 2) ** 1.5)
    z0 = z(p0)
    q = [phi(z0 + (z0 + z(p)) / (1 - a * (z0 + z(p)))) for p in (alpha / 2, 1 - alpha / 2)]
    return pct(q[0]), pct(q[1])


def add_effects(part, measure, table):
    """table: one row per participant, one column per dose."""
    for dose in DOSES[1:]:
        EFFECTS.append((part, measure, dose) + hedges_gz(table['placebo'], table[dose]))


def by_dose(df, column):
    return df.pivot(index='participant', columns='dose', values=column)[DOSES]


def save(fig, name):
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / f'{name}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


def forest(ax, rows):
    """g_z and 95 % CI of 0.075 mg and 0.150 mg against placebo, one line per measure."""
    e = pd.DataFrame([r for r in EFFECTS if r[1] in rows], columns=['part', 'measure', 'dose', 'g', 'lo', 'hi', 'k'])
    for k, dose in enumerate(DOSES[1:]):
        m = e[e['dose'] == dose]
        y = [rows.index(r) + (k - 0.5) * 0.32 for r in m['measure']]
        ax.errorbar(m['g'], y, xerr=[m['g'] - m['lo'], m['hi'] - m['g']], fmt='o', color=DOSE_COLOR[dose],
                    label=f'{dose} vs placebo')
    ax.axvline(0, color='#8C8C8C', ls='--', lw=0.8)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(-2.2, 2.2)
    ax.set_xlabel("Hedges' $g_z$ (95 % CI)")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.0))


# ---------------------------------------------------------------- 1  introduction: pupil size tracks cortical activity
def pupil_and_exponent():
    """Every 30-s window of the 12 naps: pupil / iris diameter against the aperiodic exponent of the EEG. The
    repeated-measures correlation (Bakdash & Marusich 2017) is the correlation left after removing each nap's mean.
    Big dots: per stage the mean +- s.e.m. of the participants' medians (participants with >= 3 windows)."""
    w = pd.read_csv(DATA / 'pupil_vs_exponent.csv')
    nap = w['participant'] + ' ' + w['dose']
    x = w['aperiodic_exponent'] - w.groupby(nap)['aperiodic_exponent'].transform('mean')
    y = w['pupil_iris'] - w.groupby(nap)['pupil_iris'].transform('mean')
    r_rm = np.sum(x * y) / np.sqrt(np.sum(x ** 2) * np.sum(y ** 2))
    print(f'1  pupil vs aperiodic exponent: {len(w)} windows of {nap.nunique()} naps, r_rm = {r_rm:.2f}')
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    for stage, color in STAGE_COLOR.items():
        s = w[w['stage'] == stage]
        if s.empty:
            continue
        ax.scatter(s['aperiodic_exponent'], s['pupil_iris'], s=5, color=color, alpha=0.3, lw=0)
        g = s.groupby('participant')
        med = g[['aperiodic_exponent', 'pupil_iris']].median()[g.size() >= 3]
        ax.errorbar(med.iloc[:, 0].mean(), med.iloc[:, 1].mean(), xerr=med.iloc[:, 0].sem(),
                    yerr=med.iloc[:, 1].sem(), fmt='o', ms=9, color=color, mec='white', label=stage)
    ax.set(xlabel='Aperiodic exponent (steeper = deeper sleep)', ylabel='Pupil / iris diameter',
           title=f'Pupil size tracks cortical activity ($r_{{rm}}$ = {r_rm:.2f})')
    ax.legend(frameon=False, ncol=4, fontsize=8)
    save(fig, 'fig1_pupil_and_exponent')


# ---------------------------------------------------------------- 2  same rhythms: panels a and b
RHYTHMS = {'pupil_isf_amplitude_pct': 'Pupil-ISF amplitude', 'pupil_isf_period_s': 'Pupil-ISF period',
           'sigma_isf_amplitude_pct': 'σ-ISF amplitude', 'sigma_isf_period_s': 'σ-ISF period',
           'spindle_density_per_h': 'Spindle density'}


def same_rhythms():
    """a: pupil size per sleep stage and dose (each participant faint, mean +- s.e.m. bold). b: the infraslow
    rhythms of pupil and sigma power, and spindle density, against placebo."""
    p = pd.read_csv(DATA / 'pupil_by_stage.csv')
    r = pd.read_csv(DATA / 'rhythms.csv')
    for column, label in RHYTHMS.items():
        add_effects('b  rhythms', label, by_dose(r, column))
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={'width_ratios': [1, 1.1], 'wspace': 0.55})
    stages = ['W', 'N2', 'N3', 'REM']
    for i, dose in enumerate(DOSES):
        t = p[p['dose'] == dose].pivot(index='participant', columns='stage', values='pupil_iris').reindex(
            columns=stages)
        x = np.arange(4) + (i - 1) * 0.13
        for _, row in t.iterrows():
            a.plot(x, row.to_numpy(), color=DOSE_COLOR[dose], alpha=0.25)
        a.errorbar(x, t.mean(), yerr=t.sem(), fmt='-o', lw=2, capsize=3, color=DOSE_COLOR[dose], label=dose)
    n = p.groupby(['stage', 'dose']).size().unstack()[DOSES]
    a.set_xticks(range(4), [f'{st}\nn = ' + '/'.join(map(str, n.loc[st])) for st in stages])
    a.set(ylabel='Pupil / iris diameter', title='a  Pupil size per stage (faint: participants)')
    a.legend(frameon=False, fontsize=8)
    forest(b, list(RHYTHMS.values()))
    b.set_title('b  Rhythms vs placebo', pad=24)
    save(fig, 'fig2_same_rhythms')


# ---------------------------------------------------------------- 3  lost coupling: panels c, d and e
def circular_smooth(y, win=10):
    """Gaussian moving average of a curve over the circular phase axis: window of 10 of the 60 bins, sd = window / 5
    (as MATLAB's smoothdata(y, 'gaussian', 10) on the curve wrapped round)."""
    offsets = np.arange(-(win // 2), win - win // 2)
    kernel = np.exp(-0.5 * (offsets / (win / 5)) ** 2)
    kernel /= kernel.sum()
    return np.array([np.sum(y[(i + offsets) % y.size] * kernel) for i in range(y.size)])


def lost_coupling():
    """c: one real pupil cycle. d, e: sigma power and spindle probability across the 60 bins of the pupil-ISF phase
    (each curve smoothed). The modulation depth is the maximum minus the minimum of a participant's curve: how much
    sigma power (in % of the nap's mean) or the chance of a spindle (in percentage points) rises and falls with the
    pupil cycle."""
    curves = pd.read_csv(DATA / 'phase_curves.csv').sort_values(['participant', 'dose', 'phase_deg'])
    cycle = pd.read_csv(DATA / 'example_cycle.csv')
    fig, ax = plt.subplots(3, 2, figsize=(8.5, 8), gridspec_kw={'height_ratios': [0.75, 1, 1], 'width_ratios': [1.6, 1],
                                                                 'hspace': 0.55, 'wspace': 0.35})
    ax[0, 0].remove()
    ax[0, 1].remove()
    c = fig.add_subplot(ax[1, 0].get_gridspec()[0, :])
    c.plot(cycle['t_s'], cycle['pupil_iris'], '.', ms=2, color='#B5BCC6')
    c.plot(cycle['t_s'], cycle['pupil_iris_lowpass'], color='#2F3640', lw=2)
    for _, row in cycle.dropna(subset=['eye_frame']).iterrows():
        c.annotate(row['eye_frame'], (row['t_s'], row['pupil_iris_lowpass']), xytext=(0, 9),
                   textcoords='offset points', ha='center', fontsize=8)
    c.set(xlabel='Time (s)', ylabel='Pupil / iris', title='c  One real pupil cycle (P1, placebo)')
    depths = []
    for row, (column, label) in enumerate([('sigma_dpct', 'σ power (Δ% of mean)'), ('spindle_pct', 'Spindles (%)')], 1):
        panel = 'de'[row - 1]
        phase = np.sort(curves['phase_deg'].unique())
        depth = {}
        for dose in DOSES:
            smooth = {who: circular_smooth(g[column].to_numpy())
                      for who, g in curves[curves['dose'] == dose].groupby('participant')}
            depth[dose] = pd.Series({who: y.max() - y.min() for who, y in smooth.items()})
            y = np.array(list(smooth.values()))
            m, e = y.mean(axis=0), y.std(axis=0, ddof=1) / np.sqrt(len(y))
            ax[row, 0].fill_between(phase, m - e, m + e, color=DOSE_COLOR[dose], alpha=0.15, lw=0)
            ax[row, 0].plot(phase, m, color=DOSE_COLOR[dose], lw=2, label=dose)
        depth = pd.DataFrame(depth)
        name = label.split(' (')[0]
        depths.append(depth.assign(measure=name))
        add_effects('c-e  coupling', f'Coupling of {name} (modulation depth)', depth)
        g, lo, hi, same = EFFECTS[-1][3:]
        change = 100 * (depth['0.150 mg'].mean() / depth['placebo'].mean() - 1)
        print(f'3  modulation depth of {name}: means ' + ' / '.join(f'{depth[d].mean():.1f}' for d in DOSES)
              + f' -> {change:+.0f} % at 0.150 mg, lower in {same} participants, g_z = {g:.2f} [{lo:.2f}, {hi:.2f}]')
        ax[row, 0].set_xticks([-180, -90, 0, 90, 180], ['trough', 'rise', 'peak', 'fall', 'trough'])
        ax[row, 0].set(ylabel=label, title=f'{panel}  {name} across the pupil cycle')
        b = ax[row, 1]
        for _, values in depth.iterrows():
            b.plot(range(3), values.to_numpy(), color='#C9CED6', zorder=1)
        for i, dose in enumerate(DOSES):
            b.scatter([i] * len(depth), depth[dose], s=18, color=DOSE_COLOR[dose], zorder=3)
            b.hlines(depth[dose].mean(), i - 0.28, i + 0.28, color=DOSE_COLOR[dose], lw=3, zorder=2)
        b.set_xticks(range(3), DOSES)
        b.set(ylabel='Modulation depth', title=f'{change:+.0f} % at 0.150 mg (lower in {same})')
    ax[1, 0].legend(frameon=False, fontsize=8)
    ax[2, 0].set_xlabel('Pupil-ISF phase')
    save(fig, 'fig3_lost_coupling')
    pd.concat(depths).rename_axis('participant').reset_index().round(2).to_csv(RESULTS / 'modulation_depth.csv',
                                                                                index=False)


# ---------------------------------------------------------------- 4  sleep stayed comparable
SLEEP = {'total_sleep_time_min': 'Total sleep time', 'sleep_efficiency_pct': 'Sleep efficiency',
         'sleep_onset_latency_min': 'Sleep-onset latency', 'wake_pct': 'Wake (%)', 'n1_pct': 'N1 (%)',
         'n2_pct': 'N2 (%)', 'n3_pct': 'N3 (%)', 'rem_pct': 'REM (%)', 'so_density_per_h': 'Slow-oscillation density',
         'so_peak_to_peak_uv': 'Slow-oscillation amplitude', 'spindle_amplitude_uv': 'Spindle amplitude',
         'spindles_in_trains_pct': 'Spindles in trains'}


def sleep():
    """Sleep architecture, slow oscillations and spindles against placebo."""
    s = pd.read_csv(DATA / 'sleep.csv')
    for column, label in SLEEP.items():
        add_effects('4  sleep', label, by_dose(s, column))
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    forest(ax, list(SLEEP.values()))
    ax.set_title('Sleep vs placebo', pad=24)
    save(fig, 'fig4_sleep')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')                    # σ in the printout on any console
    pupil_and_exponent()
    same_rhythms()
    lost_coupling()
    sleep()
    table = pd.DataFrame(EFFECTS, columns=['part', 'measure', 'dose', 'g_z', 'ci_low', 'ci_high', 'same_direction'])
    table.round(3).to_csv(RESULTS / 'effects.csv', index=False)
    print(f'{len(table)} effect sizes -> results/effects.csv; figures -> results/')
