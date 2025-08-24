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
# (Num_Threads parameter is optional, defaults to max)
$ ./TomlMain ../config/<config_file>.toml [Num_Threads]

# Direct mesh-based steady execution
$ ./17_fem1d_17_fem1d [2-3]d <path to .msh file>

# Direct mesh-based time-dependend execution
$ ./17_fem1d_17_fem1d_td [2-3]d <path to .msh file>

# Run tests
$ ctest
```

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

### Example: 2D Time-Dependent Diffusion-Reaction
```toml
[problem]
dimension = 2
mesh_file = "../mesh/mesh-square-h0.050000_gmsh22.msh"
time_dependent = true

[equation]
diffusion_coefficient = 1.0
reaction_coefficient = 1.0

[time_dependent]
final_time = 1.0
time_step = 0.01
theta = 0.5  # Crank-Nicolson
initial_condition = "sin(2.0 * pi * (x + y))"
forcing_function_td = "sin(2.0 * pi * (x + y)) * 2.0 * pi * cos(2.0 * pi*t)"

[[boundary_conditions]]
type = "dirichlet"
tag = 0
time_function = "sin(2.0 * pi * (x + y)) * sin(2.0 * pi*t)"
```

### Available Test Cases
- **Heat equation**: `heat_equation_1d_td.toml`
- **Advection-diffusion**: `advection_diffusion.toml`  
- **Reaction-diffusion**: `reaction_diffusion_1d.toml`
- **Poisson problems**: `poisson.toml`
- **3D manufactured solutions**: `manufactured_3d_td.toml`


<!--
## Performance & Validation
### Numerical Accuracy
- **Spatial convergence**: O(h²) for linear elements
- **Temporal convergence**: O(dt) for implicit Euler, O(dt²) for Crank-Nicolson
- **Manufactured solutions**: Machine precision for polynomial exact solutions
- **Complex geometries**: Consistent accuracy on unstructured meshes

### Parallel Performance  
- **OpenMP scaling**: Near-linear speedup up to 8 cores
- **Thread safety**: Lock-free assembly using thread-local storage
- **Memory efficiency**: Sparse matrix formats (COO → CSR conversion)
- **Cache optimization**: Element-wise loop blocking

### Validation Results
- **Method of manufactured solutions**: Verified convergence rates
- **Benchmark problems**: Agreement with analytical solutions
- **Stability analysis**: CFL condition respected for explicit schemes
- **Energy conservation**: Verified for conservative problems
-->

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

<!--
## Results Gallery
-->

