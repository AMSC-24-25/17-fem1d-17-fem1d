#include <iostream>
#include <fstream>
#include "fem.hpp"
#include "fem_td.hpp"            // la classe time-dependent che ti ho dato
#include "grid.hpp"
#include "boundary_conditions.hpp"
#include "boundary_conditions_td.hpp"

using std::cout;
using std::endl;

constexpr double PI = EIGEN_PI;

/**
 * Simple FEM solver for 1D and 2D problems
 * Usage: ./fem 1d L N  or  ./fem 2d mesh.msh
 */
int main(int argc, char *argv[])
{
    cout << "-------------17-FEM1D PROJECT-----------" << endl;
    if(argc < 3){
        cout << "Usage: " << argv[0] << " [1d L N] or [2d mesh.msh] or [3d mesh.msh]" << endl;
        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " 1d 5 20" << endl;
        cout << "  " << argv[0] << " 2d square.msh" << endl;
        cout << "  " << argv[0] << " 3d cube.msh" << endl;
        return -1;
    }


    if (argv[1][0] == '1') {
    // 1D case: fix argument parsing
        if (argc < 4) {
            cout << "1D usage: " << argv[0] << " 1d L N" << endl;
            return -1;
        }
        
        Grid1D grid(0, atof(argv[2]), atoi(argv[3]));

        Function<1,1> forcing(
            [](Point<1> p) -> double { //value
                // return sin(2*PI*p[0]);
                return (4.0 * M_PI * M_PI * 2.0 + 5.0) * sin(2.0 * M_PI * p[0]) 
                + 2.0 * M_PI * 4.0 * cos(2.0 * M_PI * p[0]);
            }
        );

        // Function<1,1> diffusion_term = OneFunction<1,1>();
        Function<1,1> diffusion_term = Function<1,1>([](Point<1> p) -> double { return 8.0; });
        Function<1,1> transport_term = Function<1,1>([](Point<1> p) -> double { return 0.5; });
        Function<1,1> reaction_term = Function<1,1>([](Point<1> p) -> double {return 2.0; });

        BoundaryConditions_td<1,1> boundary_conditions;
        boundary_conditions.addDirichlet(0, [](Point<1> p, double t) -> double { return 0; });
        boundary_conditions.addNeumann(1, [](Point<1> p, double t) -> double { return 2*M_PI*std::cos(2.0 * M_PI * p[0]) * t; });

        OrderTwoQuadrature<1> quadrature;
        FemTD<1> femtd(grid, diffusion_term, transport_term, reaction_term, boundary_conditions, quadrature);

        femtd.set_forcing([](const Point<1>& p, double t) {
            double exact = std::sin(2.0 * M_PI * p[0]) * t;
            double exactD1 = 2*M_PI*std::cos(2.0 * M_PI * p[0]) * t;
            double exactD2 = -4*M_PI*M_PI*std::sin(2.0 * M_PI * p[0]) * t;
            return 8*exactD2 + 0.5*exactD1 + 2*exact;
        });

        femtd.set_initial_condition(Function<1,1>([](const Point<1>& p){
            return 0;
        }));

        const double T = 1.0;     // tempo finale
        const double dt = 0.025;   // passo
        const double theta = 1; // 0=Esplicito, 1=Implicit Euler, 0.5=Crank–Nicolson

        // Output file prefix (facoltativi)
        femtd.run(T, dt, theta, /*vtu_prefix*/ "output/u", /*csv_prefix*/ "");
        // system("python ../scripts/plot_sol.py");
    }
    else if (argv[1][0] == '2') {
        // 2D case
        BoundaryConditions_td<2,1> boundary_conditions;
        Function<2,1> diffusion([](Point<2> p) { return 1.0; });
        Function<2,1> reaction([](Point<2> p) { return 1.0; });
        Function<2,2> transport([](Point<2> p) { return Point<2>(-0.0, -0.0); });

        // 2. Configurazione delle condizioni al contorno PRIMA del parsing

        OrderTwoQuadrature<2> quadrature;

        fun_td<2,1> exactFunc = [](Point<2> p, double t) -> double {
            return std::sin(2.0 * M_PI * p[0] * p[1]) * t;
        };

        // Configurazione con mix di Dirichlet e Neumann
        boundary_conditions.addDirichlet(0, exactFunc);
        boundary_conditions.addDirichlet(1, exactFunc);
        boundary_conditions.addDirichlet(2, exactFunc);
        boundary_conditions.addDirichlet(3, exactFunc);

        cout << "Boundary conditions:" << endl;
        cout << "  Tag 0: Dirichlet u = 0.0" << endl;
        cout << "  Tag 1: Dirichlet u = 0.0" << endl;
        cout << "  Tag 2: Dirichlet u = 0.0" << endl;
        cout << "  Tag 3: Dirichlet u = 0.0" << endl;

        Grid<2> grid;
        grid.parseFromMsh(argv[2]);

        FemTD<2> femtd(grid, diffusion, transport, reaction, boundary_conditions, quadrature);

        femtd.set_forcing([](const Point<2>& p, double t) -> double {
            double exact = std::sin(2.0 * M_PI * p[0] * p[1]) * t;
            double exactD2 = -4*M_PI*M_PI*std::sin(2.0 * M_PI * p[0] * p[1]) * t;
            return 1*exactD2 + 1*exact;
        });

        femtd.set_initial_condition(Function<2,1>([](const Point<2>&){
            return 0.0;
        }));

        // Pre-assembla M e K (tempo-indipendenti)
        femtd.assemble_time_invariant();

        const double T = 1.0;     // tempo finale
        const double dt = 1e-2;   // passo
        const double theta = 0.5; // 0=Esplicito, 1=Implicit Euler, 0.5=Crank–Nicolson

        // Output file prefix (facoltativi)
        femtd.run(T, dt, theta, /*vtu_prefix*/ "output/u", /*csv_prefix*/ "");
    }
    else if (argv[1][0] == '3'){
        // 3D time-dependent case
        cout << "=== 3D Time-Dependent Problem Configuration ===" << endl;

        // 1. Problem functions definition
        Function<3,1> diffusion([](Point<3> p) { return 1.0; });
        Function<3,1> reaction([](Point<3> p) { return 0.0; });
        Function<3,3> transport([](Point<3> p) { return Point<3>(0.0, 0.0, 0.0); });

        // 2. Configure boundary conditions BEFORE mesh parsing
        BoundaryConditions_td<3,1> boundary_conditions;

        // Manufactured solution: u(x, y, z, t) = sin(2π x y z) * t
        auto exact_sol = [](const Point<3>& p, double t) -> double {
            return std::sin(2.0 * M_PI * p[0] * p[1] * p[2]) * t;
        };

        // Dirichlet su alcune facce (ad esempio tag 0, 1, 2)
        boundary_conditions.addDirichlet(0, exact_sol);
        boundary_conditions.addDirichlet(1, exact_sol);
        boundary_conditions.addDirichlet(2, exact_sol);

        boundary_conditions.addNeumann(5, [](const Point<3>& p, double t) {
            double xyz = p[0]*p[1]*p[2];
            double pi2 = 2.0*M_PI;
            double cos_term = std::cos(pi2*xyz);
            return pi2 * p[0]*p[1] * cos_term * t; // dudx
        });
        boundary_conditions.addNeumann(3, [](const Point<3>& p, double t) {
            double xyz = p[0]*p[1]*p[2];
            double pi2 = 2.0*M_PI;
            double cos_term = std::cos(pi2*xyz);
            return pi2 * p[0]*p[2] * cos_term * t; // dudy
        });
        boundary_conditions.addNeumann(4, [](const Point<3>& p, double t) {
            double xyz = p[0]*p[1]*p[2];
            double pi2 = 2.0*M_PI;
            double cos_term = std::cos(pi2*xyz);
            return -pi2 * p[0]*p[1] * cos_term * t; // dudz
        });

        // cout << "Boundary conditions configured:" << endl;
        // cout << "  Physical tag 0 (back face): Dirichlet u = 0" << endl;
        // cout << "  Physical tag 1 (front face): Dirichlet u = 0" << endl;
        // cout << "  Physical tag 2 (left face): Dirichlet u = 0" << endl;
        // cout << "  Physical tag 3 (right face): Dirichlet u = 0" << endl;
        // cout << "  Physical tag 4 (bottom face): Dirichlet u = 0" << endl;
        // cout << "  Physical tag 5 (top face): Neumann g = 1.0" << endl;
        // cout << "NOTE: 'front' is considered as the face whose normal is (1, 0, 0)." << endl;


        cout << "Boundary conditions configured:" << endl;
        cout << "  Physical tag 0: Dirichlet u = exact" << endl;
        cout << "  Physical tag 1: Dirichlet u = exact" << endl;
        cout << "  Physical tag 2: Dirichlet u = exact" << endl;
        cout << "  Physical tag 3: Neumann du/dn = exact" << endl;
        cout << "  Physical tag 4: Neumann du/dn = exact" << endl;
        cout << "  Physical tag 5: Neumann du/dn = exact" << endl;

        OrderTwoQuadrature<3> quadrature;

        // 3. Mesh parsing
        Grid<3> grid;
        grid.parseFromMsh(argv[2]);

        cout << "Mesh parsed successfully:" << endl;
        cout << "  Number of elements: " << grid.getNumElements() << endl;
        cout << "  Number of nodes: " << grid.getNumNodes() << endl;

        // 4. Create the time-dependent FEM solver
        FemTD<3> femtd(grid, diffusion, transport, reaction, boundary_conditions, quadrature);

        // 5. Set manufactured solution forcing: u(x,y,z,t) = sin(2π x y z) * t
        femtd.set_forcing([](const Point<3>& p, double t) -> double {
            double xyz = p[0] * p[1] * p[2];
            double s = std::sin(2.0 * M_PI * xyz);
            double exact = s * t;
            double pi2 = 2.0 * M_PI;
            double pi2_2 = pi2 * pi2;
            double dsdx = pi2 * p[1] * p[2] * std::cos(2.0 * M_PI * xyz);
            double dsdy = pi2 * p[0] * p[2] * std::cos(2.0 * M_PI * xyz);
            double dsdz = pi2 * p[0] * p[1] * std::cos(2.0 * M_PI * xyz);
            double d2sdx2 = -pi2_2 * p[1] * p[1] * p[2] * p[2] * s;
            double d2sdy2 = -pi2_2 * p[0] * p[0] * p[2] * p[2] * s;
            double d2sdz2 = -pi2_2 * p[0] * p[0] * p[1] * p[1] * s;
            double laplacian = d2sdx2 + d2sdy2 + d2sdz2;
            double exactD2 = t * laplacian;
            double exactDt = s; // time derivative
            return exactDt - exactD2;
        });

        // 6. Set initial condition
        femtd.set_initial_condition(Function<3,1>([](const Point<3>&){
            return 0.0;  // u(x,y,z,0) = 0
        }));

        // 7. Pre-assemble time-invariant matrices
        femtd.assemble_time_invariant();

        // 8. Time stepping parameters
        const double T = 1.0;     // final time (smaller for 3D)
        const double dt = 1e-2;   // time step
        const double theta = 0.5; // Crank-Nicolson

        cout << "Starting 3D time-dependent simulation:" << endl;
        cout << "  Final time T = " << T << endl;
        cout << "  Time step dt = " << dt << endl;
        cout << "  Theta method = " << theta << " (0=Explicit, 0.5=Crank-Nicolson, 1=Implicit)" << endl;

        // 9. Run the simulation
        femtd.run(T, dt, theta, /*vtu_prefix*/ "output/u3d", /*csv_prefix*/ "");
        
        cout << "3D time-dependent simulation completed!" << endl;
    }
    else {
        cout << "First argument must be '1d' or '2d' or '3d'" << endl;
        return -1;
    }

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'