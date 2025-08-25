# Multi-Dimensional Finite Element Method (FEM)
C++ implementation for solving PDEs in 1D, 2D, and 3D

**Authors**: Federico Quartieri, Alessandro Ruzza, Daniele Salvi \
**Professor**: Luca Formaggia \
**Course**: Advanced Methods for Scientific Computing (2024/25)\
**Institution**: Politecnico di Milano  \
**License**: MIT License

## How to run

To clone the repository, run the following:
```bash
$ git clone https://github.com/AMSC-24-25/17-fem1d-17-fem1d.git
$ cd 17-fem1d-17-fem1d
$ git submodule update --init --recursive
```

To compile, run the following:
```bash
$ mkdir build     # if build directory does not exist
$ cd build
$ cmake ..
$ make [-j N]      # -j flag enables parallel compilation
```

To execute, run one of the following after compiling:
```bash
# Time-dependent and Static problems with TOML configuration  
# (num_threads parameter is optional, defaults to max)
$ ./TomlMain ../config/<config_file>.toml [num_threads]

# Direct uniform 1D grid steady execution (grid spans [0-length])
$ ./17_fem1d_17_fem1d 1d <length> <numNodes> [num_threads]
# Direct mesh-based steady execution
$ ./17_fem1d_17_fem1d [2-3]d <path to .msh file> [num_threads]

# Direct uniform 1D grid time-dependent execution (grid spans [0-length])
$ ./17_fem1d_17_fem1d_td 1d <length> <numNodes> [num_threads]
# Direct mesh-based time-dependent execution
$ ./17_fem1d_17_fem1d_td [2-3]d <path to .msh file> [num_threads]

# Run tests
$ ctest
```

### Run the benchmark
From inside `build/`, run the script on a directory of `.toml` configs (example set lives in `../config/speedup-analysis`). The script times each run, sweeps several thread counts, and writes a CSV one folder up (project root).


#### Optional parameters:
```
$ export THREADS="2 4 8 12 16" # thread counts to test
$ export REPEAT=3 # repetitions per (config,threads)
```

#### Run (from build/):
```
$ [Optional parameters] bash ../speedup_analysis/run_speedup.sh ./TomlMain ../config/speedup-analysis ../speedup_analysis/speedup_results.csv ./sequentialTomlMain
```

What it does
- For each `*.toml` in the target directory:
- Runs `sequentialTomlMain <config>`.
- Runs `TomlMain <config> <nThreads>` for each value in `$THREADS`.
- Captures elapsed time.
- Writes CSV with columns: `config,mode,threads,run,elapsed_s`


#### Process results (optional averaging & derived metrics)
If you ran multiple repetitions (`REPEAT>1`), you can condense to means/std and compute speedup/efficiency:

```
$ cd speedup_analysis
$ python3 parsecsv.py -i speedup_results.csv -o speedup_results_agg.csv
```

#### Plot timings by mesh
Plots mean elapsed time per mesh, with one curve per (mode,threads). Works on the raw CSV.

```
$ python3 plot.py speedup_results_agg.csv times_by_mesh.png
```

#### Plot speedup vs threads
Plots speedup (relative to the sequential time for each mesh) as a function of threads. Uses the raw CSV (single run per point is fine).

```
$ python3 plot_speedup_vs_threads.py speedup_results_agg.csv speedup_vs_threads.png
```


### Scripts

The `scripts/` directory contains helpful scripts for building the project, running all TOML test cases, and executing examples from the `config/good_examples/` folder.
It also contains the python script we used to obtain the quadrature point values.

## What is it?

This is a comprehensive C++ implementation of the Finite Element Method (FEM) capable of solving partial differential equations in 1D, 2D, and 3D domains. The implementation supports both time-dependent and steady-state problems with flexible boundary conditions.

**Key Features:**
- **Multi-dimensional support**: 1D segments, 2D triangular meshes, 3D tetrahedral meshes
- **Time-dependent problems**: Theta method
- **Parallel computing**: OpenMP support 
- **Flexible configuration**: TOML-based problem specification (at runtime)
- **Advanced quadrature**: Multiple integration orders for flexible assembly
- **Mixed boundary conditions**: Dirichlet, Neumann, and time-dependent conditions

## Architecture Overview

The project follows a modular, template-based design that enables compile-time optimization and dimensional flexibility.

### Core Classes

#### **Point<dim>**
Template class representing points in dim-dimensional space. Supports vector operations and coordinate access with compile-time dimension checking.

#### **Cell<dim>** 
Represents finite elements (intervals, triangles, tetrahedra) with:
- Barycentric coordinate transformations
- Jacobian computation for integration
- Shape function gradient calculation
- Volume/area/length measurement

#### **Grid<dim>**
Manages the computational mesh:
- Reads GMSH format files (.msh)
- Stores connectivity information
- Provides element and boundary access
- Supports mixed element types

#### **Function<dim, returnDim>**
Template-based function representation with:
- Analytical and numerical function support
- Automatic gradient computation
- Operator overloading for function arithmetic
- Time-dependent function support via `fun_td<dim, returnDim>`

#### **BoundaryConditions\<dim\>**
Comprehensive boundary condition handling system:
- **Dirichlet conditions**: Specify exact function values on boundary
- **Neumann conditions**: Specify exact flux values on boundary 
- **Time-dependent conditions**: Boundary values that vary with both space and time
- **Mixed conditions**: Combinations of Dirichlet and Neumann on different boundaries
- **Tag-based assignment**: GMSH physical tag integration for complex geometries

### Finite Element Solvers

#### **Fem\<dim\>** (Static Problems)
Solves steady-state PDEs of the form:
```
-∇·(D∇u) + b·∇u + cu = f
```
- **Assembly**: Builds stiffness matrix and load vector
- **Quadrature**: High-order Gaussian integration
- **Solving**: Direct (SparseLU) and iterative (BiCGSTAB) solvers
- **Parallelization**: OpenMP-accelerated assembly

#### **FemTD\<dim>** (Time-Dependent Problems)  
Solves parabolic PDEs:
```
∂u/∂t - ∇·(D∇u) + b·∇u + cu = f(x,t)
```
- **Time stepping**: θ-method (implicit Euler, Crank-Nicolson)
- **Mass matrix**: Consistent mass formulation
- **Stability**: Adaptive time stepping support
- **Thread safety**: Thread-local expression pools for parallel assembly

### Advanced Features

#### **Quadrature System**
Template-based numerical integration:
- **OrderTwoQuadrature**: 3-point barycentric rules for triangles/tetrahedra
- **OrderFourQuadrature**: 6-point barycentric rules for higher accuracy
- **GaussLegendre**: 1D and 2D quadrature for boundary integrals
- **Adaptive**: Future support for higher-order quadratures.

#### **Expression Parsing**
Mathematical expression evaluation using `exprtk`:
- **ThreadExpressionPool**: Thread-safe expression evaluation
- **TOML integration**: Configuration-driven problem setup
- **Runtime flexibility**: Change coefficients without recompilation

#### **Parallel Computing**
OpenMP-based parallelization:
- **Matrix assembly**: Parallel loop over elements
- **Thread-local storage**: Eliminates race conditions
- **Load balancing**: Dynamic scheduling for irregular meshes
- **Scalability**: Tested up to 16 threads with good efficiency

## Configuration System

Problems are specified using TOML configuration files located in the `config/` directory:

### Example: Complete Time-Dependent Configuration
```toml
# Complete Time-Dependent Example - showcases all available TOML options with mixed boundary conditions on a 2D square domain
# Exact solution: u(x,y,t) = x*y*t

[problem]
dimension = 2                                           # Spatial dimension (1, 2, or 3)
mesh_file = "../mesh/mesh-square-h0.050000_gmsh22.msh"  # Path to GMSH mesh file
output_file = "output/complete_example_td"              # Output file prefix
time_dependent = true                                   # Enable time-dependent solver
# For 1D problems, you can also specify:
# 1d_start = 0.0                                # 1D grid start coordinate
# 1d_end = 1.0                                  # 1D grid end coordinate  
# 1d_size = 100                                 # Number of 1D grid points

[quadrature]
type = "order2"                                 # Quadrature rule: "order2" or "order4"

[equation]
# PDE coefficients as mathematical expressions (support x, y, z, pi, e)
# NOTE: the quotes "" are important and must always be used
diffusion_function = "1.0"                      # Diffusion coefficient D(x,y,z)
transport_function_x = "0.0"                    # x-component of transport b_x(x,y,z)
transport_function_y = "0.0"                    # y-component of transport b_y(x,y,z)
transport_function_z = "0.0"                    # z-component of transport b_z(x,y,z)
reaction_function = "x"                         # Reaction coefficient c(x,y,z)
forcing_function = "0"                          # Static forcing (ignored for TD problems)

[time_dependent]
final_time = 1.0                                # Simulation end time
time_step = 0.01                                # Time step dt
theta = 0.5                                     # Theta-method: 0=Explicit, 0.5=Crank-Nicolson, 1=Implicit
initial_condition = "0.0"                       # u0(x,y,z) = 0 at t=0
forcing_function_td = "x*y + x*x*y*t"           # Time-dependent forcing term f(x,y,z,t)

# Boundary conditions (can have multiple, mixed types)
# Exact solution: u(x,y,t) = x*y*t
[[boundary_conditions]]
type = "neumann"                              # Dirichlet: specify u = g
tag = 0                                         # Boundary tag from mesh (x=0)
function = "0.0"                                # Fallback static function
time_function = "-y*t"                          # (-du/dx)             

[[boundary_conditions]]  
type = "dirichlet"                              # Dirichlet: specify u = g
tag = 1                                         # Different boundary tag (x=1)
function = "0.0"                                # Fallback static function
time_function = "x*y*t"                         # Exact solution: u(x,y,t) = x*y*t

[[boundary_conditions]]
type = "dirichlet"                              # Another Dirichlet condition
tag = 2                                         # Another boundary tag (y=0)
function = "0.0"                                # Fallback static function
time_function = "0.0"                         # Exact solution: u(x,0,t) = 0

[[boundary_conditions]]
type = "neumann"                              # Another Dirichlet condition  
tag = 3                                         # Final boundary tag (y=1)
function = "0.0"                                # Fallback static function
time_function = "x*t"                           # (du/dy)    
```

### Available Test Cases
- **Complete example**: `complete_example_td.toml` - Comprehensive time-dependent example with all options
- **1D Heat equation**: `heat_equation_1d_td.toml`
- **Advection-diffusion**: `advection_diffusion.toml`  
- **Reaction-diffusion**: `reaction_diffusion_1d.toml`
- **Poisson problems**: `poisson.toml`
- **3D manufactured solutions**: `manufactured_3d_td.toml`
- **Variable coefficients**: `variable_coefficients.toml`
- **Oscillatory problems**: `oscillatory_problem.toml`
- **Transport dominated**: `transport_dominated.toml`
- **Ring geometry**: `ring_center_05_05.toml`
- **2D Ring time-dependent**: `ring2d_td_center_00.toml`
- **Heat conduction**: `heat_conduction.toml`
- **2D time-dependent diffusion**: `diffusion_2d_td.toml`

See `config/good-examples/` for additional validated test cases with known analytical solutions.


## Performance & Validation
### Numerical Accuracy
- **Spatial convergence**: O(h²) for linear elements
- **Temporal convergence**: O(dt) for implicit Euler, O(dt²) for Crank-Nicolson
- **Manufactured solutions**: Machine precision for polynomial exact solutions
- **Complex geometries**: Consistent accuracy on unstructured meshes

### Parallel Performance  
- **OpenMP scaling**: speedup evaluated up to 16 cores
- **Thread safety**: Lock-free assembly using thread-local storage
- **Memory efficiency**: Sparse matrix formats (COO → CSR conversion)
- **Cache optimization**: Element-wise loop blocking

### Result Validation  
- **Method of manufactured solutions**: Agreement with analytical solutions

## Dependencies
### External Dependencies

- **C++20**: Modern language features and concepts
- **Eigen3**: Linear algebra and sparse matrix operations  
- **OpenMP**: Parallel computing support

### Dependencies included in this repository
- **GoogleTest**: Unit testing framework
- **TOML11**: Configuration file parsing
- **exprtk**: Mathematical expression evaluation

### Dependency for custom mesh generation
- **GMSH**: Mesh generation and I/O
## Future Development

### Possible future Extensions
- **Higher-order elements**: Quadratic and cubic shape functions
- **Adaptive mesh refinement**: Error-driven h-refinement
- **Domain decomposition**: MPI-based parallel computing
- **Multigrid solvers**: Efficient solution of large systems
- **Nonlinear problems**: Newton-Raphson iteration
- **Fluid dynamics**: Navier-Stokes equations

## Results Gallery

### Visualizations and Data
- [Google Drive Folder](https://drive.google.com/your-folder-link) - Contains simulation results, plots, and speedup data.


