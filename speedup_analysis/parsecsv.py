#!/usr/bin/env python3
# parsecsv.py

# python3 parsecsv.py -i ../speedup_result -o parsed_results.csv 

import argparse
import os
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Average multiple runs from speedup_results.csv")
    ap.add_argument("--input", "-i", required=True, help="Path to speedup_results.csv")
    ap.add_argument("--output", "-o", default=None, help="Output CSV path (default: <input> with _avg suffix)")
    ap.add_argument("--keep-path", action="store_true",
                    help="Keep full config path instead of just basename")
    args = ap.parse_args()

    in_csv = args.input
    out_csv = args.output or os.path.splitext(in_csv)[0] + "_avg.csv"

    # Leggi CSV, ignora eventuali righe rotte (es. linee isolate con numeri)
    df = pd.read_csv(in_csv, on_bad_lines="skip")

    # Normalizza colonne
    expected_cols = {"config","mode","threads","run","elapsed_s"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV is missing columns: {missing}")

    # Conversions & cleaning
    df["threads"] = pd.to_numeric(df["threads"], errors="coerce")
    df["run"]     = pd.to_numeric(df["run"], errors="coerce")
    df["elapsed_s"] = pd.to_numeric(df["elapsed_s"], errors="coerce")
    df = df.dropna(subset=["config","mode","threads","run","elapsed_s"])

    # Usa solo il nome file della config (come volevi nel CSV)
    if args.keep_path:
        df["config_key"] = df["config"]
    else:
        df["config_key"] = df["config"].apply(lambda p: os.path.basename(str(p)))

    # Aggrega
    grp_cols = ["config_key","mode","threads"]
    agg = (df.groupby(grp_cols)["elapsed_s"]
             .agg(mean="mean", std="std", count="count")
             .reset_index())

    # Calcola speedup vs sequenziale 1-thread della stessa config
    ref = (agg[(agg["mode"]=="sequential") & (agg["threads"]==1)]
             .loc[:, ["config_key","mean"]]
             .rename(columns={"mean":"seq1_mean"}))
    merged = agg.merge(ref, on="config_key", how="left")
    merged["speedup_vs_seq1"] = merged["seq1_mean"] / merged["mean"]
    merged["efficiency"] = merged["speedup_vs_seq1"] / merged["threads"]

    # Ordina per comodità
    merged = merged.sort_values(by=["config_key","mode","threads"]).reset_index(drop=True)

    # Salva
    merged.to_csv(out_csv, index=False, float_format="%.4f")
    print(f"[OK] Wrote: {out_csv}")

    # Stampa un riassunto rapido
    with pd.option_context("display.max_rows", None, "display.width", 120):
        print(merged)

if __name__ == "__main__":
    main()
