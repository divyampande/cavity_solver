#!/usr/bin/env python3
"""
Post-processing script for Lid-Driven Cavity simulation
Generates all required plots as per assignment requirements
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10

def read_field_data(filename):
    """Read 2D field data from output file"""
    data = np.loadtxt(filename, comments='#')
    x = data[:, 0]
    y = data[:, 1]
    values = data[:, 2]
    
    # Determine grid dimensions
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    # Reshape to 2D arrays
    X = x.reshape(ny, nx)
    Y = y.reshape(ny, nx)
    Z = values.reshape(ny, nx)
    
    return X, Y, Z

def read_centerline_data(filename):
    """Read 1D centerline data"""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1]

def read_ghia_data():
    """Read Ghia et al. benchmark data"""
    ghia_data = {}
    
    with open('data/ghia_data.txt', 'r') as f:
        lines = f.readlines()
    
    # Parse u-velocity data
    u_start = None
    v_start = None
    
    for i, line in enumerate(lines):
        if 'u-velocity along vertical centerline' in line:
            u_start = i + 2  # Skip header lines
        elif 'v-velocity along horizontal centerline' in line:
            v_start = i + 2
    
    # Read u-velocity
    y_u_ghia = []
    u_ghia = []
    for i in range(u_start, v_start - 2):
        line = lines[i].strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) == 2:
                y_u_ghia.append(float(parts[0]))
                u_ghia.append(float(parts[1]))
    
    # Read v-velocity
    x_v_ghia = []
    v_ghia = []
    for i in range(v_start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) == 2:
                x_v_ghia.append(float(parts[0]))
                v_ghia.append(float(parts[1]))
    
    ghia_data['y_u'] = np.array(y_u_ghia)
    ghia_data['u'] = np.array(u_ghia)
    ghia_data['x_v'] = np.array(x_v_ghia)
    ghia_data['v'] = np.array(v_ghia)
    
    return ghia_data

def plot_streamlines(output_dir='results'):
    """Plot 1: Streamlines"""
    X, Y, psi = read_field_data(f'{output_dir}/psi.dat')
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Streamlines
    levels = np.linspace(psi.min(), psi.max(), 30)
    cs = ax.contour(X, Y, psi, levels=levels, colors='black', linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=8, fmt='%.2e')
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Streamlines (Stream Function ψ)', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/streamlines.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: streamlines.png")
    plt.close()

def plot_vorticity(output_dir='results'):
    """Plot 2: Vorticity contours"""
    X, Y, omega = read_field_data(f'{output_dir}/omega.dat')
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Vorticity contours
    levels = np.linspace(omega.min(), omega.max(), 40)
    cf = ax.contourf(X, Y, omega, levels=levels, cmap='RdBu_r')
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label('Vorticity (ω)', rotation=270, labelpad=20)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Vorticity Contours', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/vorticity.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: vorticity.png")
    plt.close()

def plot_pressure(output_dir='results'):
    """Plot 3: Pressure distribution"""
    X, Y, p = read_field_data(f'{output_dir}/p.dat')
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Pressure contours
    levels = np.linspace(p.min(), p.max(), 40)
    cf = ax.contourf(X, Y, p, levels=levels, cmap='viridis')
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label('Pressure (p)', rotation=270, labelpad=20)
    
    # Add contour lines
    cs = ax.contour(X, Y, p, levels=10, colors='black', linewidths=0.5, alpha=0.5)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Pressure Distribution', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pressure.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: pressure.png")
    plt.close()

def plot_u_velocity_centerline(output_dir='results'):
    """Plot 4: u-velocity along vertical centerline (x=0.5)"""
    y, u = read_centerline_data(f'{output_dir}/u_centerline.dat')
    
    # Read Ghia data
    try:
        ghia = read_ghia_data()
        has_ghia = True
    except:
        print("Warning: Could not read Ghia et al. data")
        has_ghia = False
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot current simulation
    ax.plot(u, y, 'b-', linewidth=2, label='Current Simulation')
    
    # Plot Ghia data if available
    if has_ghia:
        ax.plot(ghia['u'], ghia['y_u'], 'ro', markersize=6, 
                label='Ghia et al. (1982)', markerfacecolor='none', markeredgewidth=1.5)
    
    ax.set_xlabel('u-velocity', fontsize=11)
    ax.set_ylabel('y', fontsize=11)
    ax.set_title('u-velocity profile at x = 0.5 (Re = 100)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)
    ax.set_xlim([-0.3, 1.1])
    ax.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/u_velocity_centerline.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: u_velocity_centerline.png")
    plt.close()

def plot_v_velocity_centerline(output_dir='results'):
    """Plot 5: v-velocity along horizontal centerline (y=0.5)"""
    x, v = read_centerline_data(f'{output_dir}/v_centerline.dat')
    
    # Read Ghia data
    try:
        ghia = read_ghia_data()
        has_ghia = True
    except:
        print("Warning: Could not read Ghia et al. data")
        has_ghia = False
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot current simulation
    ax.plot(x, v, 'b-', linewidth=2, label='Current Simulation')
    
    # Plot Ghia data if available
    if has_ghia:
        ax.plot(ghia['x_v'], ghia['v'], 'ro', markersize=6,
                label='Ghia et al. (1982)', markerfacecolor='none', markeredgewidth=1.5)
    
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('v-velocity', fontsize=11)
    ax.set_title('v-velocity profile at y = 0.5 (Re = 100)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)
    ax.set_xlim([0, 1])
    ax.set_ylim([-0.3, 0.3])
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/v_velocity_centerline.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: v_velocity_centerline.png")
    plt.close()

def plot_velocity_vectors(output_dir='results'):
    """Additional plot: Velocity vectors"""
    X, Y, U = read_field_data(f'{output_dir}/u.dat')
    _, _, V = read_field_data(f'{output_dir}/v.dat')
    
    # Subsample for clarity
    skip = 3
    X_sub = X[::skip, ::skip]
    Y_sub = Y[::skip, ::skip]
    U_sub = U[::skip, ::skip]
    V_sub = V[::skip, ::skip]
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Velocity magnitude as background
    magnitude = np.sqrt(U**2 + V**2)
    cf = ax.contourf(X, Y, magnitude, levels=30, cmap='coolwarm', alpha=0.7)
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label('Velocity Magnitude', rotation=270, labelpad=20)
    
    # Velocity vectors
    ax.quiver(X_sub, Y_sub, U_sub, V_sub, scale=5, width=0.003, headwidth=4)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Velocity Field', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/velocity_field.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: velocity_field.png")
    plt.close()

def main():
    """Main post-processing routine"""
    output_dir = 'results'
    
    print("\n" + "="*50)
    print("Lid-Driven Cavity Post-Processing")
    print("="*50 + "\n")
    
    # Check if results exist
    if not os.path.exists(f'{output_dir}/psi.dat'):
        print("Error: Results not found. Run the simulation first.")
        return
    
    print("Generating plots...")
    print()
    
    # Generate all required plots
    plot_streamlines(output_dir)
    plot_vorticity(output_dir)
    plot_pressure(output_dir)
    plot_u_velocity_centerline(output_dir)
    plot_v_velocity_centerline(output_dir)
    plot_velocity_vectors(output_dir)
    
    print()
    print("="*50)
    print("Post-processing complete!")
    print(f"All plots saved in '{output_dir}/' directory")
    print("="*50)
    print()

if __name__ == '__main__':
    main()
