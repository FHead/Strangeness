#!/usr/bin/env bash

set -euo pipefail

base_dir="$(cd "$(dirname "$0")" && pwd)"
mc_input="${MC_INPUT:-../../../../Samples/merged_mc_v2.3.root}"
data_hist_input="${DATA_HIST_INPUT:-PhiSBHistogramsData.root}"

run_case() {
  local family="$1"
  local label="$2"
  local signal_model="$3"
  local fit_min="$4"
  local fit_max="$5"
  local match_angle="$6"

  local out_dir="${base_dir}/Systematics/${family}/${label}"
  mkdir -p "${out_dir}/SignalOnlyFits" "${out_dir}/DataFits"

  "${base_dir}/ExecuteMakePhiSignalOnlyHistograms" \
    --input "${mc_input}" \
    --output "${out_dir}/PhiSignalOnlyHistograms.root" \
    --mass-min 0.99 \
    --mass-max 1.06 \
    --match-angle-max "${match_angle}"

  "${base_dir}/ExecuteFitPhiSignalOnlyShapes" \
    --input "${out_dir}/PhiSignalOnlyHistograms.root" \
    --output-dir "${out_dir}/SignalOnlyFits" \
    --fit-min "${fit_min}" \
    --fit-max "${fit_max}"

  "${base_dir}/ExecuteFitPhiSBData" \
    --signal-input "${out_dir}/PhiSignalOnlyHistograms.root" \
    --data-input "${data_hist_input}" \
    --output-dir "${out_dir}/DataFits" \
    --signal-model "${signal_model}" \
    --fit-min "${fit_min}" \
    --fit-max "${fit_max}"

  {
    echo "family=${family}"
    echo "label=${label}"
    echo "signal_model=${signal_model}"
    echo "fit_min=${fit_min}"
    echo "fit_max=${fit_max}"
    echo "match_angle_max=${match_angle}"
    echo "mc_input=${mc_input}"
    echo "data_hist_input=${data_hist_input}"
  } > "${out_dir}/manifest.txt"
}

run_case "SignalFunction" "Nominal" "GaussPlusRightTailCB" "0.99" "1.06" "0.01"
run_case "SignalFunction" "TripleGaussian" "TripleGaussian" "0.99" "1.06" "0.01"

run_case "FitRange" "Nominal" "GaussPlusRightTailCB" "0.99" "1.06" "0.01"
run_case "FitRange" "Range100to105" "GaussPlusRightTailCB" "1.00" "1.05" "0.01"

run_case "MatchingAngle" "Nominal" "GaussPlusRightTailCB" "0.99" "1.06" "0.01"
run_case "MatchingAngle" "Match0025" "GaussPlusRightTailCB" "0.99" "1.06" "0.025"
