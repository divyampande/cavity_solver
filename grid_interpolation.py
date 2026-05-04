#!/usr/bin/env python3
"""
Grid Convergence Table Generator
Reads centerline velocity data from multiple grids, performs linear interpolation,
and outputs a fully formatted LaTeX table.
"""

import numpy as np
import os
import sys

# ==========================================
# Configuration
# ==========================================
# Add or remove your grid sizes here (e.g., if you run a 201x201 grid, just add 201)
grids = [51, 76, 101, 126, 151]

# The folder where your data files are stored
data_dir = (
    "results/grid_study"  # Change this to "results" or "results/grid_study" if needed
)

# The locations where you want to sample the velocity
target_locs = [0.25, 0.50, 0.75, 0.95]


def generate_latex_table():
    u_data = {}
    v_data = {}

    # 1. Read and Interpolate Data
    for g in grids:
        # Define file paths based on your naming convention
        u_file = os.path.join(data_dir, f"u_centerline_{g}x{g}.dat")
        v_file = os.path.join(data_dir, f"v_centerline_{g}x{g}.dat")

        # Check if files exist before trying to read
        if not os.path.exists(u_file) or not os.path.exists(v_file):
            print(f"Error: Could not find files for grid {g}x{g}. Checked paths:")
            print(f"  - {u_file}")
            print(f"  - {v_file}")
            sys.exit(1)

        # Read u-velocity (y, u)
        data_u = np.loadtxt(u_file, comments="#")
        y_u, u_val = data_u[:, 0], data_u[:, 1]
        u_data[g] = np.interp(target_locs, y_u, u_val)

        # Read v-velocity (x, v)
        data_v = np.loadtxt(v_file, comments="#")
        x_v, v_val = data_v[:, 0], data_v[:, 1]
        v_data[g] = np.interp(target_locs, x_v, v_val)

    # 2. Build the LaTeX Table Header
    num_cols = len(grids) + 1
    col_format = "@{}l" + "c" * len(grids) + "@{}"

    print("%" + "=" * 50)
    print("% COPY AND PASTE THIS INTO YOUR .TEX FILE")
    print("%" + "=" * 50)
    print("\\begin{table}[!t]")
    print("\\centering")
    print("\\caption{Grid Convergence Analysis: Velocity Values at Selected Points}")
    print("\\label{tab:convergence}")
    print(f"\\begin{{tabular}}{{{col_format}}}")
    print("\\toprule")

    # Build the dynamic header row based on grids
    header_cols = ["\\textbf{Location}"] + [
        f"\\textbf{{{g}$\\times${g}}}" for g in grids
    ]
    print(" & ".join(header_cols) + " \\\\")

    # 3. Build u-velocity section
    print("\\midrule")
    print(f"\\multicolumn{{{num_cols}}}{{c}}{{\\textit{{u-velocity at x = 0.5}}}} \\\\")
    print("\\midrule")

    for i, loc in enumerate(target_locs):
        row = [f"y = {loc:.2f}"]
        for g in grids:
            row.append(f"{u_data[g][i]:8.5f}")
        print(" & ".join(row) + " \\\\")

    # 4. Build v-velocity section
    print("\\midrule")
    print(f"\\multicolumn{{{num_cols}}}{{c}}{{\\textit{{v-velocity at y = 0.5}}}} \\\\")
    print("\\midrule")

    for i, loc in enumerate(target_locs):
        row = [f"x = {loc:.2f}"]
        for g in grids:
            row.append(f"{v_data[g][i]:8.5f}")
        print(" & ".join(row) + " \\\\")

    # 5. Close table
    print("\\bottomrule")
    print("\\end{tabular}")
    print("\\end{table}")


if __name__ == "__main__":
    generate_latex_table()
