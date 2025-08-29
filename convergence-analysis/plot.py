#!/usr/bin/env python3
# plot.py

# Usage: python plot.py <input_csv> <output_png>

# Example: python plot.py "convergence_results.csv" "convergence_plot.png"

import sys, os, re
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

def main(csv_path: str = "convergence_results.csv", out_png: str = "convergence_plot.png"):
    # Load CSV for convergence analysis
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    needed = {"config", "mesh_size", "relative_error"}
    df = df[[c for c in df.columns if c in needed]].copy()
    
    # Coerce numerics; drop bad rows
    df["mesh_size"] = pd.to_numeric(df["mesh_size"], errors="coerce")
    df["relative_error"] = pd.to_numeric(df["relative_error"], errors="coerce")
    df = df.dropna(subset=["config", "mesh_size", "relative_error"])
    
    # Sort by mesh size (h value)
    df = df.sort_values("mesh_size")
    
    # Plot convergence: log(error) vs log(h)
    plt.figure(figsize=(10, 6))
    
    # Log-log plot
    plt.loglog(df["mesh_size"], df["relative_error"], 'o-', label="Relative Error", markersize=8, linewidth=2)
    
    # Add reference lines for common convergence rates
    h_min, h_max = df["mesh_size"].min(), df["mesh_size"].max()
    err_min, err_max = df["relative_error"].min(), df["relative_error"].max()
    
    # Reference line for O(h^2) convergence
    h_ref = np.logspace(np.log10(h_min), np.log10(h_max), 100)
    # Scale reference line to be visible
    err_ref_scale = err_max / (h_max**2)
    err_ref_h2 = err_ref_scale * h_ref**2
    plt.loglog(h_ref, err_ref_h2, '--', alpha=0.7, label="O(h²)", color='red')
    
    # Reference line for O(h^1) convergence
    err_ref_scale_h1 = err_max / h_max
    err_ref_h1 = err_ref_scale_h1 * h_ref
    plt.loglog(h_ref, err_ref_h1, '--', alpha=0.7, label="O(h¹)", color='green')
    
    plt.xlabel("Mesh size (h)")
    plt.ylabel("Relative Error")
    plt.title("Convergence Analysis: Relative Error vs Mesh Size")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved convergence plot -> {out_png}")
    
    # Calculate convergence rate
    if len(df) >= 2:
        # Fit line to log(error) vs log(h) to get convergence rate
        log_h = np.log(df["mesh_size"])
        log_err = np.log(df["relative_error"])
        coeffs = np.polyfit(log_h, log_err, 1)
        convergence_rate = coeffs[0]
        print(f"Estimated convergence rate: O(h^{convergence_rate:.2f})")

if __name__ == "__main__":
    csv = sys.argv[1] if len(sys.argv) >= 2 else "convergence_results.csv"
    png = sys.argv[2] if len(sys.argv) >= 3 else "convergence_plot.png"
    main(csv, png)
