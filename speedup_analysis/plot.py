#!/usr/bin/env python3
import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

def last_int(s: str):
    """Extract last integer in a string (used to order meshes)."""
    m = re.search(r'(\d+)(?!.*\d)', s)
    return int(m.group(1)) if m else None

def main(csv_path: str = "speedup_results.csv", out_png: str = "speedup_plot.png"):
    # Load CSV and drop junk rows (e.g., trailing 'f' line)
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    needed = {"config_key","mode","threads","mean"}
    df = df[[c for c in df.columns if c in needed]].copy()
    # Coerce numerics; drop bad rows
    df["threads"] = pd.to_numeric(df["threads"], errors="coerce")
    df["mean"]    = pd.to_numeric(df["mean"],    errors="coerce")
    df = df.dropna(subset=["config_key","mode","threads","mean"])

    # Normalize config basename and mesh name
    df["config"] = df["config_key"].apply(os.path.basename)
    df["mesh"]   = df["config"].str.replace(".toml", "", regex=False)

    # Build label "sequential (1)" or "openmp (N)"
    def mk_label(row):
        m = str(row["mode"]).strip().lower()
        if m == "sequential":
            return "sequential (1)"
        return f"openmp ({int(row['threads'])})"
    df["label"] = df.apply(mk_label, axis=1)

    # Average in case there are multiple lines for the same (mesh, label)
    agg = (
        df.groupby(["mesh","label"], as_index=False)["mean"]
          .mean()
          .rename(columns={"mean":"elapsed_mean"})
    )

    # Order meshes by the last integer in the name; fallback alphabetical
    unique_meshes = agg["mesh"].unique()
    sort_map = {}
    for m in unique_meshes:
        n = last_int(m)
        sort_map[m] = (0, n) if n is not None else (1, m)
    agg["sort_key"] = agg["mesh"].map(sort_map)
    agg = agg.sort_values(["sort_key","label"]).drop(columns=["sort_key"])

    # Pivot to have one curve per label
    pivot = agg.pivot(index="mesh", columns="label", values="elapsed_mean")

    # Plot
    plt.figure(figsize=(10, 5))
    for col in pivot.columns:
        plt.plot(pivot.index, pivot[col], marker="o", label=col)

    plt.xlabel("mesh (config)")
    plt.ylabel("elapsed time [s]")
    plt.title("Timing per mesh (mean over runs)")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="mode (threads)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved plot -> {out_png}")
    try:
        plt.show()
    except Exception:
        pass

if __name__ == "__main__":
    csv = sys.argv[1] if len(sys.argv) >= 2 else "speedup_results.csv"
    png = sys.argv[2] if len(sys.argv) >= 3 else "speedup_plot.png"
    main(csv, png)
