#!/usr/bin/env bash
# Usage:
#   bash ../speedup_analysis/run_speedup.sh <omp_bin> <config_dir> <out_csv> [sequential_bin]
# Example (from build folder):
#   bash ../speedup_analysis/run_speedup.sh ./TomlMain ../config/speedup-analysis ../speedup_results.csv ./sequentialTomlMain
# Optional env vars:
#   THREADS="2 4 8 12 16"
#   REPEAT=1
# REPEAT=3 bash ../speedup_analysis/run_speedup.sh ./TomlMain ../config/speedup-analysis ../speedup_results.csv ./sequentialTomlMain
# REPEAT=3 THREADS="4 8 16" bash ../speedup_analysis/run_speedup.sh ./TomlMain ../config/speedup-analysis ../speedup_results.csv ./sequentialTomlMain

set -uo pipefail  # (don't stop at errors)

OMP_BIN="${1:-./TomlMain}"
CFG_DIR="${2:-../config/speedup-analysis}"
OUT_CSV="${3:-../speedup_results.csv}"
SEQ_BIN="${4:-./sequentialTomlMain}"

THREADS="${THREADS:-2 4 8 12 16}"
REPEAT="${REPEAT:-1}"

to_abs() {
  if command -v readlink >/dev/null 2>&1; then readlink -f "$1"; else
    python3 - "$1" <<'PY'
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
  fi
}

OMP_BIN="$(to_abs "$OMP_BIN")"
SEQ_BIN="$(to_abs "$SEQ_BIN")"
OUT_CSV="$(to_abs "$OUT_CSV")"
CFG_DIR="$(to_abs "$CFG_DIR")"

mapfile -t CFGS < <(find "$CFG_DIR" -type f -name '*.toml' | sort)
if ((${#CFGS[@]} == 0)); then
  echo "No .toml files under: $CFG_DIR" >&2
  exit 1
fi

# CSV header (overwrites)
echo "config,mode,threads,run,elapsed_s" > "$OUT_CSV"

# Robust timer: uses /usr/bin/time with output to temporary file (only seconds)
time_cmd() {
  local out_tmp="$(mktemp)"
  # -f '%e' = only elapsed seconds; stdout of the program discarded, stderr to run log
  /usr/bin/time -f '%e' -o "$out_tmp" "$@" 1>/dev/null
  local rc=$?
  if ((rc != 0)); then
    rm -f "$out_tmp"
    return $rc
  fi
  local t
  t="$(cat "$out_tmp")"
  rm -f "$out_tmp"
  printf '%s' "$t"
}

for cfg in "${CFGS[@]}"; do
  cfg_dir="$(dirname "$cfg")"
  cfg_base="$(basename "$cfg")"     
  echo "[INFO] Config: $cfg"

  pushd "$cfg_dir" >/dev/null

  # --- SEQUENTIAL ---
  for r in $(seq 1 "$REPEAT"); do
    echo "=== sequential | $cfg_base | run $r/$REPEAT ==="
    log="logs-${cfg_base//\//_}-seq-$r.log"
    # execute sequential
    if t="$(time_cmd "$SEQ_BIN" "$cfg_base" 2> "$log")"; then
      echo "$cfg_base,sequential,1,$r,$t" >> "$OUT_CSV"
    else
      echo "Sequential run failed: $cfg_base (run $r) — see $cfg_dir/$log" >&2
      echo "$cfg_base,sequential,1,$r,NA" >> "$OUT_CSV"
    fi
  done

  # --- OPENMP ---
  for th in $THREADS; do
    for r in $(seq 1 "$REPEAT"); do
      echo "=== openmp | $cfg_base | ${th}t | run $r/$REPEAT ==="
      log="logs-${cfg_base//\//_}-omp-${th}-r${r}.log"
      # pass ENV OMP_NUM_THREADS as argument <nThreads>
      if t="$(OMP_NUM_THREADS="$th" time_cmd "$OMP_BIN" "$cfg_base" "$th" 2> "$log")"; then
        echo "$cfg_base,openmp,$th,$r,$t" >> "$OUT_CSV"
      else
        echo "OpenMP run failed: $cfg_base with $th threads (run $r) — see $cfg_dir/$log" >&2
        echo "$cfg_base,openmp,$th,$r,NA" >> "$OUT_CSV"
      fi
    done
  done

  popd >/dev/null
done

echo "All runs completed. Results saved to: $OUT_CSV"
