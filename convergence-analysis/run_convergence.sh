#!/usr/bin/env bash
# Usage: (from build folder)
#   bash ../convergence-analysis/run_convergence.sh <fem_bin> <config_dir> <out_csv>
# Example: (from build folder)
#   bash ../convergence-analysis/run_convergence.sh ./TomlMain ../config/convergence-analysis ./convergence_results.csv

set -uo pipefail  # (don't stop at errors)

FEM_BIN="${1:-./TomlMain}"
CFG_DIR="${2:-../config/convergence-analysis}"
OUT_CSV="${3:-../convergence-analysis/convergence_results.csv}"

to_abs() {
  if command -v readlink >/dev/null 2>&1; then readlink -f "$1"; else
    python3 - "$1" <<'PY'
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
  fi
}

FEM_BIN="$(to_abs "$FEM_BIN")"
OUT_CSV="$(to_abs "$OUT_CSV")"
CFG_DIR="$(to_abs "$CFG_DIR")"

mapfile -t CFGS < <(find "$CFG_DIR" -type f -name '*.toml' | sort -r)
if ((${#CFGS[@]} == 0)); then
  echo "No .toml files under: $CFG_DIR" >&2
  exit 1
fi

# CSV header (overwrites)
echo "config,mesh_size,relative_error" > "$OUT_CSV"

# Error extractor: captures program output to log file
extract_errors() {
  local log_file="$1"
  shift
  # capture both stdout and stderr to log file
  "$@" > "$log_file" 2>&1
  local rc=$?
  if ((rc != 0)); then
    return $rc
  fi
  
  # Extract relative error from log file
  local relative_error="NA"
  
  if [[ -f "$log_file" ]]; then
    # Look for relative error (adjust pattern based on your program's output format)
    relative_error=$(grep -i "relative.*error\|error.*relative" "$log_file" | grep -oE '[0-9]+\.?[0-9]*([eE][+-]?[0-9]+)?' | head -1)
    
    # If no specific relative error pattern found, look for generic error patterns
    if [[ -z "$relative_error" ]]; then
      relative_error=$(grep -i "error" "$log_file" | grep -oE '[0-9]+\.?[0-9]*([eE][+-]?[0-9]+)?' | head -1)
    fi
  fi
  
  [[ -z "$relative_error" ]] && relative_error="NA"
  
  printf '%s' "$relative_error"
}

for cfg in "${CFGS[@]}"; do
  cfg_dir="$(dirname "$cfg")"
  cfg_base="$(basename "$cfg")"     
  echo "[INFO] Config: $cfg"

  # Extract mesh size from filename and convert to proper decimal format
  # Patterns: h00125 -> 0.0125, h0025 -> 0.025, h005 -> 0.05, h01 -> 0.1
  mesh_size_raw=$(echo "$cfg_base" | grep -oE 'h[0-9]+' | sed 's/h//')
  if [[ -n "$mesh_size_raw" ]]; then
    case "$mesh_size_raw" in
      "00125") mesh_size="0.0125" ;;
      "0025")  mesh_size="0.025" ;;
      "005")   mesh_size="0.05" ;;
      "01")    mesh_size="0.1" ;;
      *)       mesh_size="$mesh_size_raw" ;;  # fallback for other patterns
    esac
  else
    mesh_size="unknown"
  fi

  pushd "$cfg_dir" >/dev/null

  # --- CONVERGENCE ANALYSIS ---
  echo "=== convergence | $cfg_base | mesh_size: $mesh_size ==="
  log="logs-${cfg_base//\//_}-conv.log"
  
  # Execute FEM solver and extract error information
  if rel_err="$(extract_errors "$log" "$FEM_BIN" "$cfg_base")"; then
    echo "$cfg_base,$mesh_size,$rel_err" >> "$OUT_CSV"
    echo "  -> Relative error: $rel_err"
  else
    echo "Convergence run failed: $cfg_base — see $cfg_dir/$log" >&2
    echo "$cfg_base,$mesh_size,NA" >> "$OUT_CSV"
  fi

  popd >/dev/null
done

echo "Convergence analysis completed. Results saved to: $OUT_CSV"
