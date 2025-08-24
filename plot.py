#!/usr/bin/env python3
import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

def last_int(s: str):
    """Estrae l'ultimo intero nel nome (per ordinare le mesh). Ritorna None se assente."""
    m = re.search(r'(\d+)(?!.*\d)', s)
    return int(m.group(1)) if m else None

def main(csv_path: str = "speedup_results.csv", out_png: str = "speedup_plot.png"):
    # Carica il CSV
    df = pd.read_csv(csv_path)

    # Normalizza: usa solo il basename della config (nel caso il CSV contenesse path interi)
    df["config"] = df["config"].apply(os.path.basename)

    # Pulisci tempi: elimina righe con 'NA' o non numeriche
    df = df[pd.to_numeric(df["elapsed_s"], errors="coerce").notna()].copy()
    df["elapsed_s"] = df["elapsed_s"].astype(float)

    # Crea label per la curva (mode + threads). Tratta il sequenziale come 1 thread.
    def mk_label(row):
        if str(row["mode"]).lower() == "sequential":
            return "sequential (1)"
        return f"openmp ({int(row['threads'])})"

    df["label"] = df.apply(mk_label, axis=1)

    # Estrai un nome mesh più corto (senza .toml)
    df["mesh"] = df["config"].str.replace(".toml", "", regex=False)

    # Media per (mesh, label)
    agg = (
        df.groupby(["mesh", "label"], as_index=False)["elapsed_s"]
          .mean()
          .rename(columns={"elapsed_s": "elapsed_mean"})
    )

    # Ordina le mesh usando l’ultimo numero nel nome, fallback alfabetico
    order_key = []
    for m in agg["mesh"].unique():
        n = last_int(m)
        order_key.append((m, (0, n) if n is not None else (1, m)))
    sort_map = {m: k for m, k in order_key}
    agg["sort_key"] = agg["mesh"].map(sort_map)
    agg = agg.sort_values(["sort_key", "label"]).drop(columns=["sort_key"])

    # Pivot per avere colonne = label (curve)
    pivot = agg.pivot(index="mesh", columns="label", values="elapsed_mean")

    # Plot
    plt.figure(figsize=(10, 5))
    for col in pivot.columns:
        plt.plot(pivot.index, pivot[col], marker="o", label=col)  # niente colori espliciti

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
        # in ambienti senza display, semplicemente non mostra
        pass

if __name__ == "__main__":
    csv = sys.argv[1] if len(sys.argv) >= 2 else "speedup_results.csv"
    png = sys.argv[2] if len(sys.argv) >= 3 else "speedup_plot.png"
    main(csv, png)
