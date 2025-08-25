#!/usr/bin/env python3
import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

def last_int(s: str):
    """Prende l'ultimo intero in una stringa (per ordinare le mesh)."""
    m = re.search(r'(\d+)(?!.*\d)', s)
    return int(m.group(1)) if m else None

def main(csv_path: str = "speedup_results.csv",
         out_png: str = "speedup_vs_threads.png"):
    # Carica CSV e tieni solo le colonne necessarie (niente medie)
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)

    cols_needed = ["config_key", "threads", "speedup_vs_seq1"]
    for c in cols_needed:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found in CSV")

    # Cast numerici
    df["threads"] = pd.to_numeric(df["threads"], errors="coerce")
    df["speedup_vs_seq1"] = pd.to_numeric(df["speedup_vs_seq1"], errors="coerce")

    # Pulisci righe non valide
    df = df.dropna(subset=["config_key", "threads", "speedup_vs_seq1"])

    # Etichetta mesh: basename senza .toml
    df["mesh"] = df["config_key"].apply(os.path.basename).str.replace(".toml", "", regex=False)

    # Se per caso ci sono duplicati (stesso mesh, stesso #thread), tieni la prima occorrenza
    df = df.drop_duplicates(subset=["mesh", "threads"], keep="first")

    # Pivot: righe=threads, colonne=mesh, valori=speedup
    pivot = df.pivot(index="threads", columns="mesh", values="speedup_vs_seq1")

    # Ordina threads e colonne mesh
    pivot = pivot.sort_index()
    mesh_order = sorted(
        pivot.columns,
        key=lambda m: (0, last_int(m)) if last_int(m) is not None else (1, m)
    )
    pivot = pivot[mesh_order]

    # Plot
    plt.figure(figsize=(10, 5))
    for col in pivot.columns:
        plt.plot(pivot.index, pivot[col], marker="o", label=col)

    plt.xlabel("threads")
    plt.ylabel("speedup vs sequential (1 thread)")
    plt.title("Speedup by threads")
    plt.xticks(pivot.index.tolist())
    plt.legend(title="mesh (config)", ncol=2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved plot -> {out_png}")
    try:
        plt.show()
    except Exception:
        pass

if __name__ == "__main__":
    csv = sys.argv[1] if len(sys.argv) >= 2 else "speedup_results.csv"
    png = sys.argv[2] if len(sys.argv) >= 3 else "speedup_vs_threads.png"
    main(csv, png)
