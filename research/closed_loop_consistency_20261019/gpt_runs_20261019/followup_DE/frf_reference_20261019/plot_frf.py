"""Plot saved FRFs and model endpoints only; no estimation or frequency selection."""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'
DEST = HERE / 'figures'
METHODS = ('M1', 'M2', 'M3', 'M3w', 'M5', 'M5n')
MAIN = ('M1', 'M3', 'M5n')
STYLES = {
    'M1': dict(color='#82459A', ls='-', lw=.95),
    'M2': dict(color='#AC87BB', ls='--', lw=.90),
    'M3': dict(color='#1765AE', ls='-', lw=1.10),
    'M3w': dict(color='#15979E', ls='--', lw=.95),
    'M5': dict(color='#3E3E3E', ls=':', lw=1.10),
    'M5n': dict(color='#387E58', ls='-.', lw=.95),
}
LABELS = {'M1': 'M1 direct, all bins', 'M2': 'M2 direct, lines',
          'M3': 'M3 reference IV', 'M3w': 'M3w weighted IV',
          'M5': 'M5 periodic SML', 'M5n': 'M5n local-noise SML'}
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 6.3, 'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6.3,
    'axes.linewidth': .55, 'lines.linewidth': .9,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'axes.spines.top': True, 'axes.spines.right': True,
    'savefig.facecolor': 'white',
})


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_npz(path):
    with np.load(path) as a:
        return {k: a[k] for k in a.files}


def phase_display(g, anchor):
    # Separate unwrap for every input/output curve, followed by one whole-turn shift.
    p = np.rad2deg(np.unwrap(np.angle(g), axis=0))
    shifts = 360 * np.round((anchor - p[0]) / 360)
    return p + shifts[None], shifts


def model_files(z, split, summary):
    models = {}
    statuses = {}
    sources = []
    for method in METHODS:
        row = next((x for x in summary['fits'] if x['split'] == split and x['method'] == method), {})
        path = RESULTS / f'z{z}_{split}_{method}_model.npz'
        statuses[method] = dict(available=path.exists(), converged=bool(row.get('converged', False)),
                                iterations=row.get('iterations'))
        if path.exists():
            models[method] = read_npz(path)
            sources.append(dict(name=path.name, sha256=digest(path)))
    return models, statuses, sources


def source_csv(z, fit_split, ref_split, ref, models, statuses, reference='lpm'):
    prefix = 'lrm_' if reference == 'lrm' else ''
    path = DEST / f'{prefix}z{z}_{fit_split}_vs_{ref_split}_source.csv'
    f = ref['f']
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['frequency_hz', 'output_1based', 'input_1based', 'series', 'real', 'imag',
                         'reference_se_complex', 'reference_se_db', 'reference_se_phase_deg',
                         'asymmetric_window', 'converged'])
        for k, fk in enumerate(f):
            for y in range(4):
                for u in range(4):
                    g = ref['G'][k, y, u]
                    writer.writerow([fk, y+1, u+1, reference.upper(), g.real, g.imag,
                                     ref['se_complex'][k,y,u], ref['se_db'][k,y,u], ref['se_phase_deg'][k,y,u],
                                     bool(ref['edge'][k]), 'not applicable'])
                    for method, a in models.items():
                        g = a['G'][k,y,u]
                        writer.writerow([fk, y+1, u+1, method, g.real, g.imag, '', '', '',
                                         bool(ref['edge'][k]), statuses[method]['converged']])
    return path


def edge_shading(ax, f, edge):
    middle = np.flatnonzero(~edge)
    if middle.size:
        if middle[0]:
            ax.axvspan(f[0], .5*(f[middle[0]-1]+f[middle[0]]), color='0.97', zorder=0)
        if middle[-1] < len(f)-1:
            ax.axvspan(.5*(f[middle[-1]]+f[middle[-1]+1]), f[-1], color='0.97', zorder=0)


def draw(z, fit_split, ref_split, ref, models, statuses, selected, panel_set, quantity, reference='lpm'):
    f = ref['f']
    present = [m for m in selected if m in models]
    for a in models.values():
        if not np.array_equal(a['f'], f):
            raise ValueError('Saved model and reference frequency grids differ; no interpolation authorized.')
    gref = ref['G']
    anchor = np.angle(gref[0]) * 180 / np.pi
    pref, shiftref = phase_display(gref, anchor)
    phases = {m: phase_display(models[m]['G'], anchor) for m in present}
    fig, axes = plt.subplots(4, 4, figsize=(183/25.4, 178/25.4), sharex=True)
    fig.subplots_adjust(left=.085, right=.988, bottom=.106, top=.853, hspace=.39, wspace=.32)
    for y in range(4):
        for u in range(4):
            ax = axes[y,u]
            edge_shading(ax, f, ref['edge'])
            ax.set_title(f'$u_{u+1} \\rightarrow y_{y+1}$', pad=3.3)
            if quantity == 'magnitude':
                center = 20*np.log10(np.abs(gref[:,y,u]))
                width = 2*ref['se_db'][:,y,u]
                ax.fill_between(f, center-width, center+width, facecolor='0.78', alpha=.38, linewidth=0, zorder=1)
                ax.plot(f, center, color='0.50', lw=1.15, zorder=2)
                for method in present:
                    ax.plot(f, 20*np.log10(np.abs(models[method]['G'][:,y,u])), **STYLES[method], zorder=3)
            elif quantity == 'phase':
                center = pref[:,y,u]
                width = 2*ref['se_phase_deg'][:,y,u]
                ax.fill_between(f, center-width, center+width, facecolor='0.78', alpha=.38, linewidth=0, zorder=1)
                ax.plot(f, center, color='0.50', lw=1.15, zorder=2)
                for method in present:
                    ax.plot(f, phases[method][0][:,y,u], **STYLES[method], zorder=3)
            elif quantity == 'complex_difference':
                denominator = np.abs(gref[:,y,u])
                with np.errstate(divide='ignore', invalid='ignore'):
                    random = 200*ref['se_complex'][:,y,u]/denominator
                    ax.plot(f, random, color='0.55', lw=.85, ls=':', zorder=2)
                    for method in present:
                        value = 100*np.abs(models[method]['G'][:,y,u]-gref[:,y,u])/denominator
                        ax.plot(f, value, **STYLES[method], zorder=3)
                ax.set_yscale('log')
                ax.axhline(100, color='.82', lw=.5, ls='--', zorder=0)
            ax.set_xscale('log')
            ax.set_xlim(float(f[0]), float(f[-1]))
            ax.set_xticks([3,5,10,15,25], labels=['3','5','10','15','25'])
            ax.tick_params(axis='both', which='major', direction='out', length=2, width=.5, pad=1.8)
            ax.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)
            ax.grid(which='major', color='.88', lw=.45)
            ax.grid(which='minor', color='.93', lw=.35, ls=':')
            ax.margins(x=0)
    ylabel = {'magnitude': 'Magnitude (dB; recorded units)',
              'phase': 'Independently unwrapped phase (deg)',
              'complex_difference': f'Complex relative difference to {reference.upper()} (%)'}[quantity]
    fig.text(.012,.485,ylabel,rotation=90,va='center',fontsize=8)
    fig.text(.54,.066,'Frequency (Hz)',ha='center',fontsize=8)
    split_label = 'Same-record check: 8-period models / 8-period reference' if fit_split==ref_split=='full' else f'{fit_split.capitalize()} 4-period models / {ref_split} 4-period reference'
    qualifier = ' (diagnostic)' if reference == 'lrm' else ''
    title = {'magnitude': f'Nonparametric {reference.upper()} reference{qualifier} and parametric models',
             'phase': f'Phase comparison with the {reference.upper()} reference{qualifier}',
             'complex_difference': f'Complex difference to the {reference.upper()} reference{qualifier}'}[quantity]
    fig.suptitle(f'{title}\n$z={z}$ mm | {split_label}', y=.982, fontsize=9, linespacing=1.5)
    if quantity == 'complex_difference':
        handles = [Line2D([],[],color='.55',ls=':',lw=.85,label=f'2 x random SE / |{reference.upper()}|')]
    else:
        handles = [Line2D([],[],color='.50',lw=1.15,label=f'{reference.upper()} reference'), Patch(facecolor='.78',alpha=.38,label='+/- 2 random SE')]
    for method in present:
        suffix = '' if statuses[method]['converged'] else ' [failed endpoint]'
        handles.append(Line2D([],[],**STYLES[method],label=LABELS[method]+suffix))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.924),ncol=3,
               columnspacing=1.3,handlelength=2.5,handletextpad=.6,frameon=False)
    footer = f'{reference.upper()} is an estimate, not ground truth; smoothing and inversion bias are not included in random SE.'
    second = 'All frequency bins retained; pale end regions use asymmetric local windows.'
    if quantity == 'phase':
        second = 'Each curve unwrapped separately; only a constant multiple of 360 deg aligns its first point.'
    if fit_split!=ref_split:
        second += ' Prior fixed delays are shared.'
    fig.text(.54,.033,footer,ha='center',fontsize=6.1)
    fig.text(.54,.017,second,ha='center',fontsize=5.8)
    prefix = 'lrm_' if reference == 'lrm' else ''
    stem = f'{prefix}z{z}_{fit_split}_vs_{ref_split}_{panel_set}_{quantity}'
    files = []
    for ext in ('png','pdf','svg'):
        path = DEST / (stem+'.'+ext)
        fig.savefig(path,dpi=300)
        files.append(dict(name=path.name,bytes=path.stat().st_size))
    plt.close(fig)
    return dict(stem=stem, quantity=quantity, n_axes=16, n_frequency=len(f),
                methods=present, missing=[m for m in selected if m not in models],
                failed=[m for m in present if not statuses[m]['converged']], files=files,
                reference_phase_shift_deg=shiftref.tolist(),
                model_phase_shift_deg={m: phases[m][1].tolist() for m in present})


def run(z, panels='main', reverse=False, reference='lpm'):
    DEST.mkdir(parents=True,exist_ok=True)
    summary_path=RESULTS/f'z{z}_summary.json'
    summary=json.loads(summary_path.read_text(encoding='utf-8'))
    reference_name = 'LRM q=1, window=31, diagnostic, data-only front-half CV selection' if reference == 'lrm' else 'LPM q=2, window=21'
    manifest=dict(z_mm=z, backend='matplotlib', reference=reference_name,
                  reference_is_truth=False, filtered_frequencies=False, sources=[], figures=[])
    pairs=[('full','full'),('front','back')]
    if reverse:
        pairs.append(('back','front'))
    for fit_split,ref_split in pairs:
        path=(HERE/'lrm_results'/f'z{z}_{ref_split}_lrm.npz') if reference == 'lrm' else RESULTS/f'z{z}_{ref_split}_q2_w21.npz'
        ref=read_npz(path)
        models,statuses,sources=model_files(z,fit_split,summary)
        manifest['sources'].append(dict(name=path.name,sha256=digest(path)))
        manifest['sources'].extend(sources)
        source=source_csv(z,fit_split,ref_split,ref,models,statuses,reference)
        manifest.setdefault('source_tables',[]).append(source.name)
        selections=[('main',MAIN)]
        if panels=='all':
            selections.append(('all_methods',METHODS))
        for label,selected in selections:
            for quantity in ('magnitude','phase','complex_difference'):
                manifest['figures'].append(draw(z,fit_split,ref_split,ref,models,statuses,selected,label,quantity,reference))
                print(manifest['figures'][-1]['stem'],flush=True)
    manifest['qa']={'same_frequency_grid':True,'all_frequencies_retained':True,
                    'outputs_nonempty':all(x['bytes']>0 for r in manifest['figures'] for x in r['files']),
                    'no_interpolation':True,'no_refit':True,'uncertainty_includes_smoothing_bias':False}
    prefix = 'lrm_' if reference == 'lrm' else ''
    (DEST/f'{prefix}z{z}_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--z',default='0')
    parser.add_argument('--panels',choices=('main','all'),default='main')
    parser.add_argument('--reverse',action='store_true')
    parser.add_argument('--reference',choices=('lpm','lrm'),default='lpm')
    args=parser.parse_args()
    for z in (('-0.5','0','0.5','1') if args.z=='all' else [args.z]):
        run(z,args.panels,args.reverse,args.reference)
