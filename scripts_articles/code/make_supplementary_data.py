"""Build the Supplementary Data workbook (S1_Data.xlsx).

Sheets follow the table numbering used in the manuscript and S1 Text:
default_weights, then T1 to T7.

Values are written at the precision stored in the deposited outputs; the
display rounding is controlled by cell number formats, not by rounding the
values themselves.

Requires, relative to this script:
  - ../data/case_results.json      (produced by run_cases.R)
  - ../data/benchmark_results.json (produced by run_benchmark.R)
  - ../data/suppl_sensitivity.json (produced by run_suppl.R)
  - ../data/suppl_loco.json        (produced by run_suppl.R)

Writes ../S1_Data.xlsx
"""
import json
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

with open('../data/case_results.json') as f:
    cases = json.load(f)
with open('../data/benchmark_results.json') as f:
    bench = json.load(f)
with open('../data/suppl_sensitivity.json') as f:
    sens = json.load(f)
with open('../data/suppl_loco.json') as f:
    loco = json.load(f)

# Commit of the deposited code and data that this workbook is built from.
# Update this after pushing, then re-run this script.
COMMIT = 'f30b904e9241e1a7f349d08721eb0f63def2f85e'

wb = openpyxl.Workbook()

HEADER_FILL = PatternFill(start_color='1F4E79', end_color='1F4E79', fill_type='solid')
HEADER_FONT = Font(name='Arial', size=10, bold=True, color='FFFFFF')
NORMAL_FONT = Font(name='Arial', size=10)
ITALIC_FONT = Font(name='Arial', size=10, italic=True)
TITLE_FONT = Font(name='Arial', size=12, bold=True)
BORDER = Border(left=Side(style='thin', color='CCCCCC'),
                right=Side(style='thin', color='CCCCCC'),
                top=Side(style='thin', color='CCCCCC'),
                bottom=Side(style='thin', color='CCCCCC'))

F4 = '0.0000'
F3 = '0.000'
F1 = '0.0'


def style_header_row(ws, row, ncol):
    for c in range(1, ncol + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal='center', vertical='center')
        cell.border = BORDER


def put_header(ws, headers, row=3):
    for c, h in enumerate(headers, 1):
        ws.cell(row=row, column=c, value=h)
    style_header_row(ws, row, len(headers))


def autosize(ws):
    for col in ws.columns:
        try:
            col_letter = col[0].column_letter
        except Exception:
            continue
        max_len = 0
        for cell in col:
            v = cell.value
            if v is not None:
                max_len = max(max_len, len(str(v)))
        ws.column_dimensions[col_letter].width = min(max_len + 2, 45)


def num(ws, r, c, value, fmt):
    """Write a numeric cell at stored precision, or NA when unavailable."""
    if value is None:
        ws.cell(row=r, column=c, value='NA')
    else:
        cell = ws.cell(row=r, column=c, value=value)
        cell.number_format = fmt


# ======================================================================
# README
# ======================================================================
ws = wb.active
ws.title = 'README'
README_ROWS = [
    (1, 'S1 Data. Supplementary tables', None),
    (2, 'Prioritizing candidate wildlife reservoirs from observational data with a composite evidence score', None),
    (4, 'Sheet', 'Contents'),
    (5, 'default_weights', 'Default weights for the fourteen criteria.'),
    (6, 'T1_scores', 'Per-host scores and 95% site-bootstrap intervals for the three outbreak reconstructions.'),
    (7, 'T2_criteria', 'Normalized scores for each host and criterion in the reconstructions.'),
    (8, 'T3_benchmark_summary', 'Deposited summaries for eleven scenarios, with 300 replicates per scenario.'),
    (9, 'T4_benchmark_perrep', 'All 3,300 deposited benchmark records, in their original order.'),
    (10, 'T5_weight_sensitivity', 'Existing Dirichlet weight-perturbation results (200 draws per reconstruction).'),
    (11, 'T6_leave_one_out', 'Existing leave-one-criterion-out results for the reconstructions.'),
    (12, 'T7_scenario_parameters', 'Scenario modifications, canonical defaults and initial conditions transcribed from the deposited code.'),
    (14, 'Results source', 'https://github.com/rocheben/spilloverID/tree/' + COMMIT + '/scripts_articles/data'),
    (15, 'Code source', 'https://github.com/rocheben/spilloverID/tree/' + COMMIT + '/R/simulate.R'),
    (16, 'Files', 'case_results.json; benchmark_results.json; suppl_sensitivity.json; suppl_loco.json. Table layout follows scripts_articles/code/make_supplementary_data.py.'),
    (17, 'Code version', 'spilloverID v0.2.0; MIT license; commit ' + COMMIT + '.'),
    (18, 'Data handling', 'Values are copied from the deposited outputs at their stored precision. Cell formats control display rounding.'),
    (19, 'Host labels', 'For the outbreak reconstructions, “designated host” identifies the comparator encoded by true_host in the source files. The wildlife maintenance role assumed in the Rift Valley fever reconstruction is unresolved.'),
    (20, 'Unavailable values', 'NA denotes an unavailable value in the source output. It is not zero. The no-spillover scenario has no true reservoir.'),
    (21, 'Intervals', 'T1 intervals use the reported site-level bootstrap. T5 quantiles describe weight perturbations. Additional uncertainty analyses are reported in S1 Text.'),
]
for r, a, b in README_ROWS:
    if a is not None:
        ws.cell(row=r, column=1, value=a)
    if b is not None:
        ws.cell(row=r, column=2, value=b)
ws['A1'].font = Font(name='Arial', size=14, bold=True)
ws['A2'].alignment = Alignment(wrap_text=True, vertical='top')
style_header_row(ws, 4, 2)
for r, a, b in README_ROWS:
    if r >= 5:
        ws.cell(row=r, column=1).font = Font(name='Arial', size=10, bold=True)
        ws.cell(row=r, column=2).font = NORMAL_FONT
        ws.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 24
ws.column_dimensions['B'].width = 110


# ======================================================================
# default_weights
# ======================================================================
ws = wb.create_sheet('default_weights')
default_weights = [
    ('C1',  'Spatial co-distribution',                        0.04),
    ('C2',  'Reservoir-versus-accidental host shape',             0.04),
    ('C3',  'Individual exposure association',                0.06),
    ('C4',  'Dose-response (residence time)',                 0.04),
    ('C5',  'Low cryptic fraction',                           0.03),
    ('C6',  'Exposure \u00d7 local prevalence interaction',        0.09),
    ('C7',  'Catalytic FOI coherence',                        0.04),
    ('C8',  'Temporal precedence',                            0.03),
    ('C9',  'Mediation',                                      0.06),
    ('C10', 'Mechanistic Akaike weights',                     0.07),
    ('C11', 'Specificity (negative-control pathogens)',       0.05),
    ('C12', 'Inter-epidemic persistence',                     0.16),
    ('C13', 'Endemic age-prevalence shape',                   0.14),
    ('C14', 'Spatial independence from co-circulating hosts', 0.15),
]
ws['A1'] = 'Default criterion weights'
ws['A1'].font = TITLE_FONT
put_header(ws, ['Criterion', 'Name', 'Default weight'])
for i, (cn, name, w) in enumerate(default_weights):
    r = 4 + i
    ws.cell(row=r, column=1, value=cn)
    ws.cell(row=r, column=2, value=name)
    ws.cell(row=r, column=3, value=w).number_format = '0.00'
ws.cell(row=18, column=2, value='Sum').font = Font(bold=True)
c = ws.cell(row=18, column=3, value='=SUM(C4:C17)')
c.number_format = '0.00'
c.font = Font(bold=True)
autosize(ws)


# ======================================================================
# T1. Scores for the outbreak reconstructions
# ======================================================================
ws = wb.create_sheet('T1_scores')
ws['A1'] = 'T1. Scores for the outbreak reconstructions'
ws['A1'].font = TITLE_FONT
put_header(ws, ['Case', 'Host', 'Pk', '95% CI lower', '95% CI upper',
                'Rank', 'Designated host'])
ws['C3'] = '\u03a0k'
r = 4
for key in ['sn', 'nv', 'rv']:
    c = cases[key]
    for i, host in enumerate(c['hosts']):
        ws.cell(row=r, column=1, value=c['case'])
        ws.cell(row=r, column=2, value=host)
        num(ws, r, 3, c['pi'][i], F4)
        num(ws, r, 4, c['ci_low'][i] if c.get('ci_low') else None, F4)
        num(ws, r, 5, c['ci_high'][i] if c.get('ci_high') else None, F4)
        ws.cell(row=r, column=6, value=c['rank'][i])
        ws.cell(row=r, column=7,
                value='Yes' if host == c['true_host'] else 'No')
        r += 1
autosize(ws)


# ======================================================================
# T2. Normalized criterion scores
# ======================================================================
ws = wb.create_sheet('T2_criteria')
ws['A1'] = 'T2. Normalized criterion scores'
ws['A1'].font = TITLE_FONT
crit_names = ['C%d' % i for i in range(1, 15)]
put_header(ws, ['Case', 'Host', 'Designated host'] + crit_names)
r = 4
for key in ['sn', 'nv', 'rv']:
    c = cases[key]
    for i, host in enumerate(c['hosts']):
        ws.cell(row=r, column=1, value=c['case'])
        ws.cell(row=r, column=2, value=host)
        ws.cell(row=r, column=3,
                value='Yes' if host == c['true_host'] else 'No')
        crit_row = c['criteria'][i]
        for j, cn in enumerate(crit_names):
            v = crit_row.get(cn)
            num(ws, r, 4 + j, None if v is None else float(v), F4)
        r += 1
autosize(ws)


# ======================================================================
# T3. Benchmark summary by scenario
# ======================================================================
ws = wb.create_sheet('T3_benchmark_summary')
ws['A1'] = 'T3. Benchmark summary by scenario'
ws['A1'].font = TITLE_FONT
put_header(ws, ['Scenario', 'Label', 'Replicates', 'Top-1 (%)', 'Top-2 (%)',
                'Mean margin', 'SD margin', 'Margin < 0 (%)', 'Median P true',
                'Median max P', 'Median P range', 'Median observable margin'])
ws['I3'] = 'Median \u03a0 true'
ws['J3'] = 'Median max \u03a0'
ws['K3'] = 'Median \u03a0 range'
r = 4
for row in bench['summary']:
    ws.cell(row=r, column=1, value=row['scenario'])
    ws.cell(row=r, column=2, value=row['label'])
    ws.cell(row=r, column=3, value=row['n_rep'])
    num(ws, r, 4, row.get('top1_rate'), F1)
    num(ws, r, 5, row.get('top2_rate'), F1)
    num(ws, r, 6, row.get('mean_margin'), F3)
    num(ws, r, 7, row.get('sd_margin'), F3)
    num(ws, r, 8, row.get('pct_margin_negative'), F1)
    num(ws, r, 9, row.get('median_pi_true'), F3)
    num(ws, r, 10, row.get('median_max_pi'), F3)
    num(ws, r, 11, row.get('median_pi_range'), F3)
    num(ws, r, 12, row.get('median_obs_margin'), F3)
    r += 1
autosize(ws)


# ======================================================================
# T4. Per-replicate benchmark outputs
# ======================================================================
ws = wb.create_sheet('T4_benchmark_perrep')
ws['A1'] = 'T4. Per-replicate benchmark outputs (%s simulations)' % format(
    len(bench['per_rep']), ',')
ws['A1'].font = TITLE_FONT
put_header(ws, ['Scenario', 'Replicate', 'True reservoir host ID', 'P true',
                'Max P other', 'Margin', 'Observable margin', 'Rank true',
                'Top-1', 'Top-2', 'Max P (all hosts)', 'P range'])
ws['D3'] = '\u03a0 true'
ws['E3'] = 'Max \u03a0 other'
ws['K3'] = 'Max \u03a0 (all hosts)'
ws['L3'] = '\u03a0 range'
r = 4
FLOAT_COLS = {4, 5, 6, 7, 11, 12}
for row in bench['per_rep']:
    vals = [row.get('scenario'), row.get('rep'), row.get('k_star'),
            row.get('pi_true'), row.get('pi_max_other'), row.get('margin'),
            row.get('obs_margin'), row.get('rank_true'), row.get('top1'),
            row.get('top2'), row.get('max_pi'), row.get('pi_range')]
    for c, v in enumerate(vals, 1):
        if c in FLOAT_COLS:
            num(ws, r, c, v, F4)
        elif v is None:
            ws.cell(row=r, column=c, value='NA')
        else:
            ws.cell(row=r, column=c, value=v)
    r += 1
autosize(ws)
ws.freeze_panes = 'A4'


# ======================================================================
# T5. Weight-perturbation sensitivity
# ======================================================================
ws = wb.create_sheet('T5_weight_sensitivity')
ws['A1'] = 'T5. Weight-perturbation sensitivity'
ws['A1'].font = TITLE_FONT
ws['A2'] = ('200 Dirichlet perturbations per reconstruction; '
            'concentration vector = default weights \u00d7 30.')
ws['A2'].font = ITALIC_FONT
ws.merge_cells('A2:H2')
ws['A2'].alignment = Alignment(wrap_text=True)
put_header(ws, ['Case', 'Designated host', 'Baseline P', 'Rank-1 (%)',
                'Top-2 (%)', 'P median', 'P 2.5%', 'P 97.5%'])
ws['C3'] = 'Baseline \u03a0'
ws['F3'] = '\u03a0 median'
ws['G3'] = '\u03a0 2.5%'
ws['H3'] = '\u03a0 97.5%'
r = 4
for key in ['sn', 'nv', 'rv']:
    s = sens[key]
    ws.cell(row=r, column=1, value=s['case'])
    ws.cell(row=r, column=2, value=s['true_host'])
    num(ws, r, 3, s['baseline_pi'], F4)
    num(ws, r, 4, s['pct_rank1'], F1)
    num(ws, r, 5, s['pct_rank_top2'], F1)
    num(ws, r, 6, s['pi_median'], F4)
    num(ws, r, 7, s['pi_2.5'], F4)
    num(ws, r, 8, s['pi_97.5'], F4)
    r += 1
autosize(ws)


# ======================================================================
# T6. Leave-one-criterion-out results
# ======================================================================
ws = wb.create_sheet('T6_leave_one_out')
ws['A1'] = 'T6. Leave-one-criterion-out results'
ws['A1'].font = TITLE_FONT
ws['A2'] = ('Each criterion is omitted in turn. \u0394\u03a0 is the change '
            'from the designated host\u2019s baseline score.')
ws['A2'].font = ITALIC_FONT
ws.merge_cells('A2:F2')
ws['A2'].alignment = Alignment(wrap_text=True)
put_header(ws, ['Case', 'Designated host', 'Criterion removed',
                'P designated host', 'Rank designated host',
                'DP from baseline'])
ws['D3'] = '\u03a0 designated host'
ws['F3'] = '\u0394\u03a0 from baseline'
r = 4
for key in ['sn', 'nv', 'rv']:
    l = loco[key]
    for crow in l['rows']:
        ws.cell(row=r, column=1, value=l['case'])
        ws.cell(row=r, column=2, value=l['true_host'])
        ws.cell(row=r, column=3, value=crow['criterion'])
        num(ws, r, 4, crow['pi_true'], F4)
        ws.cell(row=r, column=5, value=crow['rank_true'])
        num(ws, r, 6, crow['delta_pi'], F4)
        r += 1
autosize(ws)
ws.freeze_panes = 'A4'


# ======================================================================
# T7. Scenario parameters and initial conditions
# ======================================================================
ws = wb.create_sheet('T7_scenario_parameters')
ws['A1'] = 'T7. Scenario parameters and initial conditions'
ws['A1'].font = TITLE_FONT
ws['A2'] = ('Scenario-specific entries override the canonical defaults. '
            'Values and expressions are transcribed from R/simulate.R and '
            'run_benchmark.R.')
ws['A2'].font = ITALIC_FONT
ws.merge_cells('A2:C2')
ws['A2'].alignment = Alignment(wrap_text=True)
put_header(ws, ['Scenario / scope', 'Parameter', 'Value or definition'])
T7_ROWS = [
    ('Canonical default', 'n_sites', 30),
    ('Canonical default', 'n_hosts', 4),
    ('Canonical default', 'true_reservoir', 1),
    ('Canonical default', 'bridge_host', 'NA'),
    ('Canonical default', 'n_humans_per_site', 60),
    ('Canonical default', 'n_animals_per_host_site', 40),
    ('Canonical default', 'n_years', 15),
    ('Canonical default', 'seasons_per_year', 4),
    ('Canonical default', 'animal_time_resolution', 2),
    ('Canonical default', 'anthrop_effect_on_abundance', 'c(0.5, -0.3, 0.0, 0.2)'),
    ('Canonical default', 'anthrop_effect_on_exposure', 0.4),
    ('Canonical default', 'anthrop_effect_on_prevalence', 0.3),
    ('Canonical default', 'base_abundance', 'c(100, 80, 60, 40)'),
    ('Canonical default', 'R0_reservoir', 2.5),
    ('Canonical default', 'cross_species_beta', 0.05),
    ('Canonical default', 'recovery_rate', 0.4),
    ('Canonical default', 'seasonal_amplitude', 0.4),
    ('Canonical default', 'seasonal_phase', 0),
    ('Canonical default', 'env_noise_sd', 0.15),
    ('Canonical default', 'alpha_reservoir', 0.05),
    ('Canonical default', 'alpha_bridge', 0),
    ('Canonical default', 'alpha_environment', 0.003),
    ('Canonical default', 'alpha_cryptic_bg', 0.0005),
    ('Canonical default', 'age_min', 5),
    ('Canonical default', 'age_max', 75),
    ('Canonical default', 'residence_share_local', 0.7),
    ('Canonical default', 'serology_sensitivity', 0.92),
    ('Canonical default', 'serology_specificity', 0.97),
    ('Canonical default', 'exposure_misclass_rate', 0.1),
    ('Canonical default', 'antibody_persistence', 999),
    ('S2_bridge', 'bridge_host', 2),
    ('S2_bridge', 'alpha_reservoir', 0.035),
    ('S2_bridge', 'alpha_bridge', 0.025),
    ('S3_cryptic', 'alpha_reservoir', 0.015),
    ('S3_cryptic', 'alpha_environment', 0.025),
    ('S3_cryptic', 'alpha_cryptic_bg', 0.003),
    ('S4_confounded', 'anthrop_effect_on_abundance', 'c(1.0, -0.5, 0.0, 0.3)'),
    ('S4_confounded', 'anthrop_effect_on_exposure', 0.8),
    ('S4_confounded', 'anthrop_effect_on_prevalence', 0.7),
    ('S5_codistributed', 'anthrop_effect_on_abundance', 'c(0.5, 0.5, 0.0, -0.2)'),
    ('S6_strong_seasonal', 'seasonal_amplitude', 0.7),
    ('S6_strong_seasonal', 'env_noise_sd', 0.25),
    ('S7_low_power', 'n_humans_per_site', 20),
    ('S7_low_power', 'n_animals_per_host_site', 15),
    ('S7_low_power', 'n_sites', 15),
    ('S8_measurement_error', 'exposure_misclass_rate', 0.3),
    ('S8_measurement_error', 'serology_sensitivity', 0.8),
    ('S8_measurement_error', 'serology_specificity', 0.9),
    ('S9_weak_signal', 'alpha_reservoir', 0.01),
    ('S9_weak_signal', 'alpha_cryptic_bg', 0.002),
    ('S10_short_immunity', 'antibody_persistence', 3),
    ('S11_no_spillover', 'alpha_reservoir', 0),
    ('S11_no_spillover', 'alpha_bridge', 0),
    ('S11_no_spillover', 'alpha_environment', 0),
    ('S11_no_spillover', 'alpha_cryptic_bg', 0.0008),
    ('S11_no_spillover', 'true_reservoir', 'NA'),
    ('Initial conditions', 'Reservoir prevalence, time 1', 'min(max((1 − 1/R0_reservoir) × (1 + Normal(0, 0.1)), 0.01), 0.9) per site'),
    ('Initial conditions', 'Other-host prevalence, time 1', '0'),
    ('Initial conditions', 'Time step', 'dt = 1 / seasons_per_year'),
    ('Benchmark', 'Replicates per scenario', 300),
    ('Benchmark', 'Seed', '20260524 + 17 × replicate + scenario index (1–11)'),
    ('S11_no_spillover', 'Baseline animal seroprevalence', '0.02 added by the observation model'),
]
r = 4
for scope, param, value in T7_ROWS:
    ws.cell(row=r, column=1, value=scope)
    ws.cell(row=r, column=2, value=param)
    ws.cell(row=r, column=3, value=value)
    r += 1
ws.column_dimensions['A'].width = 30
ws.column_dimensions['B'].width = 43
ws.column_dimensions['C'].width = 100


wb.save('../S1_Data.xlsx')
print('Wrote S1_Data.xlsx (%d sheets).' % len(wb.sheetnames))
