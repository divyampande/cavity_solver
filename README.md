# Lid-Driven Cavity Solver
## Stream Function-Vorticity Formulation

A Fortran implementation of the lid-driven cavity problem using the stream function-vorticity method with Python visualization.

---

## Problem Description

Solves the 2D incompressible Navier-Stokes equations for flow in a square cavity:
- Domain: 1×1 unit square
- Top wall (lid): moving with velocity u = 1 unit/s
- Other walls: stationary (no-slip)
- Reynolds number: Re = 100 (configurable)

### Governing Equations (Non-dimensional)

Stream function-vorticity formulation:
- Vorticity transport: ∂ω/∂t + u·∇ω = (1/Re)∇²ω
- Stream function: ∇²ψ = -ω
- Velocities: u = ∂ψ/∂y, v = -∂ψ/∂x
- Vorticity: ω = ∂v/∂x - ∂u/∂y

---

## Directory Structure

```
cavity_solver/
├── src/
│   ├── main.f90         # Main program driver
│   ├── solver.f90       # Numerical solver routines
│   ├── boundary.f90     # Boundary conditions
│   └── input.f90        # Input file reader
├── results/             # Output directory (created automatically)
├── data/
│   └── ghia_data.txt    # Benchmark data from Ghia et al. (1982)
├── cavity.in            # Input parameters
├── postprocess.py       # Visualization script
├── Makefile             # Build system
└── README.md            # This file
```

---

## Requirements

### Fortran Compilation
- **gfortran** (GNU Fortran compiler)
  - Install on Ubuntu/Debian: `sudo apt install gfortran`
  - Install on macOS: `brew install gcc`

### Python Visualization
- **Python 3.x** with:
  - numpy
  - matplotlib

Install Python dependencies:
```bash
pip install numpy matplotlib
```

---

## Building and Running

### 1. Compile the Code

```bash
make
```

This will:
- Create the `build/` directory
- Compile all Fortran modules
- Link the executable `cavity_solver`

### 2. Configure Parameters

Edit `cavity.in` to adjust simulation parameters:

```
nx  50              # Grid points in x (try 50, 75, 101)
ny  50              # Grid points in y
Re  100.0           # Reynolds number
dt  0.001           # Time step
omega_sor  1.7      # SOR relaxation parameter (1.0-2.0)
max_time_steps  100000
tol_steady  1.0e-6  # Convergence tolerance
```

### 3. Run the Simulation

```bash
./cavity_solver
```

Or use the Makefile:
```bash
make run
```

The solver will:
- Read parameters from `cavity.in`
- Iterate until steady state is reached
- Write results to `results/` directory

**Expected output files:**
- `psi.dat` - Stream function
- `omega.dat` - Vorticity
- `u.dat`, `v.dat` - Velocity components
- `p.dat` - Pressure
- `u_centerline.dat`, `v_centerline.dat` - Centerline profiles

### 4. Visualize Results

```bash
python3 postprocess.py
```

This generates plots in `results/`:
1. `streamlines.png` - Streamlines
2. `vorticity.png` - Vorticity contours
3. `pressure.png` - Pressure distribution
4. `u_velocity_centerline.png` - u(y) at x=0.5 vs. Ghia et al.
5. `v_velocity_centerline.png` - v(x) at y=0.5 vs. Ghia et al.
6. `velocity_field.png` - Velocity vectors (bonus)

---

## Grid Convergence Study

To perform grid refinement as required by the assignment:

```bash
# 51×51 grid
sed -i 's/^nx.*/nx  50/' cavity.in
sed -i 's/^ny.*/ny  50/' cavity.in
make run
python3 postprocess.py
mv results/u_velocity_centerline.png results/u_centerline_51x51.png

# 76×76 grid
sed -i 's/^nx.*/nx  75/' cavity.in
sed -i 's/^ny.*/ny  75/' cavity.in
make run
python3 postprocess.py
mv results/u_velocity_centerline.png results/u_centerline_76x76.png

# 101×101 grid
sed -i 's/^nx.*/nx  100/' cavity.in
sed -i 's/^ny.*/ny  100/' cavity.in
make run
python3 postprocess.py
mv results/u_velocity_centerline.png results/u_centerline_101x101.png
```

---

## Numerical Method Details

### Stream Function-Vorticity Method

**Advantages:**
- Automatically satisfies continuity (incompressibility)
- No pressure-velocity coupling issues
- Well-suited for 2D problems

**Algorithm:**
1. Initialize fields (ψ=0, ω=0)
2. For each time step:
   a. Solve ∇²ψ = -ω using SOR
   b. Compute velocities from ψ
   c. Advance ω using explicit time-stepping
   d. Apply boundary conditions
3. Iterate until steady state

**Key Features:**
- **Poisson solver**: Successive Over-Relaxation (SOR)
- **Vorticity advection**: Upwind differencing (for stability)
- **Vorticity diffusion**: Central differencing
- **Boundary vorticity**: Derived from no-slip condition

---

## Validation

Results are compared against benchmark data from:

> Ghia, U., Ghia, K.N. and Shin, C.T. (1982)  
> "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method"  
> *Journal of Computational Physics*, 48, 387-411.  
> DOI: [10.1016/0021-9991(82)90058-4](http://dx.doi.org/10.1016/0021-9991(82)90058-4)

---

## Troubleshooting

### Compilation Issues
- **Error: gfortran not found**
  → Install gfortran: `sudo apt install gfortran`

- **Module not found errors**
  → Clean and rebuild: `make clean && make`

### Convergence Issues
- **Not converging to steady state**
  → Reduce time step `dt` in `cavity.in`
  → Increase `max_time_steps`
  → Try different SOR parameter (1.5-1.9)

- **Poisson solver not converging**
  → Increase `max_iter_psi`
  → Reduce `tol_psi`

### Stability Issues
- **Solution blowing up**
  → Reduce `dt` (CFL condition)
  → Ensure upwind scheme is active
  → Check boundary conditions

---

## Assignment Deliverables Checklist

- [x] Grid size: 51×51 (can increase to 101×101)
- [x] Reynolds number: Re = 100
- [x] Non-dimensional governing equations
- [x] Stream function-vorticity method
- [x] Output plots:
  - [x] 1. Streamlines
  - [x] 2. Vorticity contours
  - [x] 3. Pressure distribution
  - [x] 4. u-velocity at x=0.5 vs y
  - [x] 5. v-velocity at y=0.5 vs x
  - [x] 6. Comparison with Ghia et al. (1982)

---

## Code Modification Guide

### Changing the Problem
To solve different cavity problems, modify:

**Different Reynolds number:**
- Edit `Re` in `cavity.in`

**Different lid velocity:**
- Edit `u_lid` in `cavity.in`

**Different domain size:**
- Currently normalized to 1×1
- To change: modify grid spacing in `main.f90`

### Adding Features

**Output additional quantities:**
1. Compute in `main.f90`
2. Add `write_field()` call before deallocation
3. Update `postprocess.py` to read and plot

**Different boundary conditions:**
1. Modify `apply_psi_bc()` and `apply_omega_bc()` in `boundary.f90`

---

## Performance Notes

Typical run times (Intel Core i7):
- 51×51 grid: ~10 seconds
- 76×76 grid: ~30 seconds
- 101×101 grid: ~2 minutes

For faster convergence:
- Optimize SOR parameter (1.7-1.9 typically best)
- Use coarser time steps for initial iterations
- Implement multigrid methods (advanced)

---

## References

1. Ghia, U., Ghia, K.N. and Shin, C.T. (1982). High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method. *J. Comput. Phys.*, 48, 387-411.

2. Anderson, J.D. (1995). *Computational Fluid Dynamics: The Basics with Applications*. McGraw-Hill.

3. Peyret, R. and Taylor, T.D. (1983). *Computational Methods for Fluid Flow*. Springer-Verlag.

---

## License

This code is provided for educational purposes as part of AM5630 coursework at IIT Madras.

---

## Author

Divyam  
M.Tech Aerospace Engineering  
IIT Madras  
May 2026

---

## Acknowledgments

- Course: AM5630 - Foundations of Computational Fluid Dynamics
- Institution: Indian Institute of Technology Madras
- Reference data: Ghia et al. (1982)
