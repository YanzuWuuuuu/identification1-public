"""Measured-data adaptation of paper-plot/scripts/line_aime.py.

Figure-polish pass: muted palette, logarithmic displacement axis, no fake
breakpoints or confidence bands, measured noise reference instead of constant.
Source: finalized distinct_N32, heldout frequency lines and fresh time records.
"""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'outputs_structured' / 'distinct_N32' / 'summary.json'
data = json.loads(SOURCE.read_text(encoding='utf-8'))['metrics']
f = np.asarray(data['frequency_hz'])
mask = (f >= 20) & (f <= 100)
cal = data['calibration']
series = {
    'static_rigid_nm': np.asarray(cal['static']['general_rigid_curve_m'])[mask]*1e9,
    'joint_rigid_nm': np.asarray(cal['joint']['general_rigid_curve_m'])[mask]*1e9,
    'joint_noise_nm': np.asarray(cal['joint']['general_noise_curve_m'])[mask]*1e9,
    'digital_oracle_rigid_nm': np.asarray(cal['oracle']['general_rigid_curve_m'])[mask]*1e9,
}
with (HERE/'lowband_source.csv').open('w', encoding='utf-8', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(['frequency_hz', *series])
    writer.writerows(zip(f[mask], *series.values()))

style_source = Path.home()/'.codex/skills/deepscientist-figure-polish/assets/deepscientist-academic.mplstyle'
# Quote hex colors in a local copy: unquoted '#' is a comment in mplstyle.
style_text = '\n'.join(line.split(': #')[0] + ': "#' + line.split(': #')[1] + '"'
                       if ': #' in line else line
                       for line in style_source.read_text(encoding='utf-8').splitlines())
style_local = HERE/'academic_local.mplstyle'
style_local.write_text(style_text, encoding='utf-8')
plt.style.use(style_local)
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
                     'text.usetex': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(5.7, 3.5))
ax.plot(f[mask], series['static_rigid_nm'], color='#998A7A', lw=1.5,
        ls='--', label='Static only: rigid leakage')
ax.plot(f[mask], series['joint_noise_nm'], color='#74787B', lw=1.4,
        ls=':', label='Joint: output noise')
ax.plot(f[mask], series['joint_rigid_nm'], color='#526E85', lw=2.0,
        marker='o', ms=3, markevery=3, label='Joint: rigid leakage')
ax.set(xlabel='Held-out frequency (Hz)', ylabel='RMS per excited line (nm)',
       xlim=(20,100), yscale='log', ylim=(0.01,30))
ax.set_xticks([20,40,60,80,100])
ax.grid(axis='y', which='major', alpha=.18, lw=.5)
ax.grid(axis='x', visible=False)
ax.legend(loc='upper right', frameon=False, fontsize=9)
fig.tight_layout(pad=.8)
for extension in ('png', 'svg'):
    fig.savefig(HERE/f'closed_loop_lowband.{extension}', dpi=300, facecolor='white')
plt.close(fig)
print('Saved measured curves: closed_loop_lowband.png / .svg')
