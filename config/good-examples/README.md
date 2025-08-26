# Good Examples - Test Problems with Known Solutions

This directory contains TOML configuration files for finite element problems with known exact solutions, designed for visual verification and testing.

## Problem Categories

### 1D Problems (4 examples)

#### `polynomial_1d.toml`
- **Exact solution**: u(x) = x^3 - x^2 + x
- **PDE coefficients**: D=2, b=1, r=0.5
- **Boundary conditions**: Dirichlet at both ends
- **Expected behavior**: Smooth cubic polynomial

#### `sinusoidal_1d.toml`
- **Exact solution**: u(x) = sin(2pi*x)
- **PDE coefficients**: D=1, b=0.5, r=2
- **Boundary conditions**: Dirichlet at x=0, Neumann at x=1
- **Expected behavior**: Smooth sine wave

#### `longer_sinusoidal_1d.toml`
- **Exact solution**: Extended sinusoidal solution to [0, 5] interval
- **PDE coefficients**: Varying parameters
- **Boundary conditions**: Mixed Dirichlet/Neumann
- **Expected behavior**: Extended sinusoidal pattern

#### `mixed_bc_1d.toml`
- **Exact solution**: u(x) = x(1-x)e^x
- **PDE coefficients**: D=0.5, b=1.5, r=1
- **Boundary conditions**: Dirichlet at x=0, Neumann at x=1
- **Expected behavior**: Parabolic shape with exponential growth

### 2D Problems (5 examples)

#### `polynomial_2d.toml`
- **Exact solution**: u(x,y) = x(1-x)\*y(1-y)
- **PDE coefficients**: D=1, b=(0.5,0), r=1
- **Boundary conditions**: Three Dirichlet faces, one Neumann
- **Expected behavior**: Smooth bilinear surface, zero on boundaries

#### `trigonometric_2d.toml`
- **Exact solution**: u(x,y) = sin(pi\*x)·sin(pi\*y)
- **PDE coefficients**: D=2, b=(1,0.5), r=0.5
- **Boundary conditions**: Dirichlet on all boundaries (naturally zero)
- **Expected behavior**: Sinusoidal dome, zero on all boundaries

#### `variable_coeffs_2d.toml`
- **Exact solution**: u(x,y) = xy(1-x)(1-y)
- **PDE coefficients**: D=(1+x+y), b=(x,y), r=xy (spatially varying)
- **Boundary conditions**: Two Dirichlet, two Neumann
- **Expected behavior**: Smooth bilinear surface, zero on boundaries

#### `exponential_2d.toml`
- **Exact solution**: u(x,y) = exp(-x)\*cos(pi\*y)
- **PDE coefficients**: D=0.5, b=(2,-1), r=3
- **Boundary conditions**: Mixed Dirichlet/Neumann
- **Expected behavior**: Exponential decay in x, cosine oscillation in y

#### `gaussian_bump_2d.toml`
- **Exact solution**: u(x,y) = exp(-((x-0.5)\^2 + (y-0.5)\^2)/0.1)
- **PDE coefficients**: D=1, b=(0.1,0.1), r=0.1
- **Boundary conditions**: Mixed Dirichlet/Neumann
- **Expected behavior**: Gaussian bump centered at (0.5, 0.5)

### 3D Problems (5 examples)

#### `polynomial_3d.toml`
- **Exact solution**: u(x,y,z) = xyz(1-x)(1-y)(1-z)
- **PDE coefficients**: D=1, b=(0,0,0), r=0 (pure diffusion)
- **Boundary conditions**: Dirichlet on all faces (u=0)
- **Expected behavior**: Smooth polynomial volume, zero on all boundaries

#### `trigonometric_3d.toml`
- **Exact solution**: u(x,y,z) = sin(pi\*x)\*sin(pi\*y)\*sin(pi\*z)
- **PDE coefficients**: D=2, b=(1,0.5,0.8), r=0.5
- **Boundary conditions**: Dirichlet on all faces (zero on 3 pairs)
- **Expected behavior**: 3D sinusoidal pattern, zero on most boundaries

#### `exponential_3d.toml`
- **Exact solution**: u(x,y,z) = exp(-x)\*cos(pi\*y)\*cos(pi\*z)
- **PDE coefficients**: D=1, b=(2,1.5,1), r=2
- **Boundary conditions**: Dirichlet on all faces
- **Expected behavior**: Exponential decay in x, cosine oscillations in y,z

#### `mixed_bc_3d.toml`
- **Exact solution**: u(x,y,z) = x\*y\*z\*(1-x)\*(1-y)\*(1-z)
- **PDE coefficients**: D=0.5, b=(1,0.8,0.6), r=1.5
- **Boundary conditions**: Mixed Dirichlet/Neumann
- **Expected behavior**: Smooth trilinear surface, zero on x-boundaries

#### `gaussian_3d.toml`
- **Exact solution**: u(x,y,z) = exp(-5\*((x-0.5)\^2 + (y-0.5)\^2 + (z-0.5)\^2))
- **PDE coefficients**: D=1, b=(0.2,0.15,0.1), r=0.3
- **Boundary conditions**: Dirichlet on all faces
- **Expected behavior**: 3D Gaussian bump centered at (0.5, 0.5, 0.5)

## Usage

These problems are designed for:
1. **Code verification**: Compare numerical solution with known exact solution
2. **Visual inspection**: Solutions have clear, recognizable patterns
3. **Boundary condition testing**: Mix of Dirichlet and Neumann conditions
4. **Coefficient testing**: Non-zero diffusion, transport, and reaction terms

## Verification Process

1. Run the FEM solver on these configurations (use TomlMain)
2. Compare numerical solution with exact analytical solution 
4. Visualize results to check for expected patterns

## Expected Error Behavior

- Polynomial solutions should converge at optimal rates
- Smooth trigonometric solutions should show spectral accuracy
- Problems with variable coefficients test coefficient handling
- Mixed boundary conditions verify BC implementation

All problems include non-trivial diffusion, transport, and reaction terms to thoroughly test the FEM implementation.
