#!/usr/bin/env bash
# ============================================================================
# Full reproduction pipeline. Run from this directory (code/).
#
#   bash run_all.sh
#
# Every step writes into ../data or ../figures and is checkpointed where it is
# long, so the pipeline can be interrupted and restarted. Total runtime is a few
# hours, dominated by the 11-scenario benchmark.
#
# R    >= 4.0.0 with lme4, jsonlite and the spilloverID package installed
#      (R CMD INSTALL ../package/spilloverID)
# Python 3 with numpy, pandas, matplotlib, openpyxl
# ============================================================================
set -euo pipefail
export BENCH_N_CORE="${BENCH_N_CORE:-$(nproc)}"
export EXT_N_CORE="${EXT_N_CORE:-$BENCH_N_CORE}"
export SWEEP_N_CORE="${SWEEP_N_CORE:-$BENCH_N_CORE}"
export PG_N_CORE="${PG_N_CORE:-$BENCH_N_CORE}"

step () { echo; echo "=== $* ==="; }

step "0. Install the package"
R CMD INSTALL ../package/spilloverID

step "1. Outbreak reconstructions"
Rscript run_cases.R

step "2. Mechanistic benchmark (long)"
Rscript run_benchmark.R

step "3. Composite versus single criteria (Suppl. Fig. F)"
Rscript comp_vs_indiv.R

step "4. Interpretation-threshold sweep (Suppl. Fig. D)"
Rscript threshold_sweep.R

step "5. Weight perturbation and leave-one-out (Suppl. Figs. H and I)"
Rscript run_suppl.R

step "6. Prospective power analysis (Suppl. Fig. J)"
Rscript power_grid.R

step "7. Bootstrap convergence (Suppl. Fig. K)"
Rscript boot_converge.R

step "8. S4 reliability weighting (Suppl. Fig. B)"
( cd suppl/S4 && Rscript 10_reliability.R && Rscript 20_reweight.R \
    && Rscript 21_check_renorm.R && python3 30_figure.py )

step "9. S5 uncertainty and threshold recalibration (Suppl. Fig. C)"
( cd suppl/S5 && Rscript 02_uncertainty.R && Rscript 03_threshold.R \
    && Rscript 04_observable_margin.R && python3 05_calibration.py \
    && Rscript 07_case_threshold.R && python3 06_figure.py )

step "10. S6 reference rankings (Suppl. Fig. E)"
( cd suppl/S6 && Rscript baselines_summary.R && python3 figure_baselines.py )

step "11. S7 taxonomic resolution (Suppl. Fig. G)"
( cd suppl/S7 && Rscript validate_ext.R && Rscript run_benchmark_ext.R \
    && python3 make_figure_ext.py )

step "12. Main and supplementary figures"
for f in make_figure_1_framework.py make_figure_2_benchmark.py \
         make_figure_3_pi_scores.py make_figure_4_criteria_heatmap.py \
         make_figure_5_margins_and_C12.py make_suppl_figure_H_sensitivity.py \
         make_suppl_figure_I_loco.py make_suppl_figure_A_per_rep_scatter.py \
         make_suppl_figures_D_F.py make_suppl_figure_J_power.py \
         make_suppl_figure_K_bootconv.py; do
  echo "  $f"; python3 "$f"
done

step "13. Supplementary Data workbook"
python3 make_supplementary_data.py

echo; echo "=== pipeline complete ==="
