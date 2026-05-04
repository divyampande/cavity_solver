#!/bin/bash
# Grid refinement study script for lid-driven cavity
# Tests different grid resolutions and compares results

echo "=========================================="
echo "Grid Refinement Study"
echo "=========================================="
echo ""

# Grid sizes to test
grids=(50 75 100 200)

# Create directory for grid study results
mkdir -p results/grid_study

for n in "${grids[@]}"; do
    echo "Running simulation with ${n}x${n} grid..."
    
    # Update input file
    sed -i "s/^nx.*/nx  $n/" cavity.in
    sed -i "s/^ny.*/ny  $n/" cavity.in
    
    # Run simulation
    ./cavity_solver > results/grid_study/log_${n}x${n}.txt
    
    # Save results
    cp results/u_centerline.dat results/grid_study/u_centerline_${n}x${n}.dat
    cp results/v_centerline.dat results/grid_study/v_centerline_${n}x${n}.dat
    
    # Generate plots
    python3 postprocess.py > /dev/null 2>&1
    
    # Save plots
    cp results/u_velocity_centerline.png results/grid_study/u_centerline_${n}x${n}.png
    cp results/v_velocity_centerline.png results/grid_study/v_centerline_${n}x${n}.png
    
    echo "  ✓ Completed ${n}x${n} grid"
    echo ""
done

echo "=========================================="
echo "Grid study complete!"
echo "Results saved in results/grid_study/"
echo "=========================================="
echo ""

# Create comparison plot
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

grids = [50, 75, 100, 200]
colors = ['blue', 'green', 'red', 'purple']

# Load Ghia data
ghia = {}
with open('data/ghia_data.txt', 'r') as f:
    lines = f.readlines()

u_start, v_start = None, None
for i, line in enumerate(lines):
    if 'u-velocity along vertical' in line:
        u_start = i + 2
    elif 'v-velocity along horizontal' in line:
        v_start = i + 2

y_u_ghia, u_ghia = [], []
for i in range(u_start, v_start - 2):
    line = lines[i].strip()
    if line and not line.startswith('#'):
        parts = line.split()
        if len(parts) == 2:
            y_u_ghia.append(float(parts[0]))
            u_ghia.append(float(parts[1]))

x_v_ghia, v_ghia = [], []
for i in range(v_start, len(lines)):
    line = lines[i].strip()
    if line and not line.startswith('#'):
        parts = line.split()
        if len(parts) == 2:
            x_v_ghia.append(float(parts[0]))
            v_ghia.append(float(parts[1]))

# Plot u-velocity comparison
for i, n in enumerate(grids):
    data = np.loadtxt(f'results/grid_study/u_centerline_{n}x{n}.dat', comments='#')
    ax1.plot(data[:,1], data[:,0], color=colors[i], linewidth=2, 
             label=f'{n}×{n} grid', alpha=0.7)

ax1.plot(u_ghia, y_u_ghia, 'ko', markersize=6, label='Ghia et al. (1982)',
         markerfacecolor='none', markeredgewidth=1.5)
ax1.set_xlabel('u-velocity', fontsize=11)
ax1.set_ylabel('y', fontsize=11)
ax1.set_title('u-velocity at x=0.5 (Grid Convergence)', fontweight='bold')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# Plot v-velocity comparison
for i, n in enumerate(grids):
    data = np.loadtxt(f'results/grid_study/v_centerline_{n}x{n}.dat', comments='#')
    ax2.plot(data[:,0], data[:,1], color=colors[i], linewidth=2,
             label=f'{n}×{n} grid', alpha=0.7)

ax2.plot(x_v_ghia, v_ghia, 'ko', markersize=6, label='Ghia et al. (1982)',
         markerfacecolor='none', markeredgewidth=1.5)
ax2.set_xlabel('x', fontsize=11)
ax2.set_ylabel('v-velocity', fontsize=11)
ax2.set_title('v-velocity at y=0.5 (Grid Convergence)', fontweight='bold')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/grid_study/grid_convergence.png', dpi=300, bbox_inches='tight')
print("Grid convergence comparison saved!")
EOF

echo "Grid convergence comparison plot created!"
echo "See: results/grid_study/grid_convergence.png"
