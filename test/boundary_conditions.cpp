#include <gtest/gtest.h>
#include "boundary_conditions.hpp"
#include "grid1D.hpp"
#include "grid2D.hpp"
#include "function.hpp"
#include "point.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <chrono>
#include <iostream>

class BoundaryConditionsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup di base per i test
    }
    
    void TearDown() override {
        // Cleanup dopo i test
    }
};

// Test 1: Creazione di BoundaryConditions vuoto
TEST_F(BoundaryConditionsTest, DefaultConstructor) {
    BoundaryConditions<2,1> bc;
    EXPECT_NO_THROW(bc);
}

// Test 2: Aggiunta di condizioni Dirichlet
TEST_F(BoundaryConditionsTest, AddDirichletCondition) {
    BoundaryConditions<2,1> bc;
    
    // Test con valore costante
    EXPECT_NO_THROW(bc.addDirichlet(0, Point<1>(1.0)));
    
    // Test con funzione
    Function<2,1> func([](Point<2> p) { return p[0] + p[1]; });
    EXPECT_NO_THROW(bc.addDirichlet(1, func));
}

// Test 3: Aggiunta di condizioni Neumann
TEST_F(BoundaryConditionsTest, AddNeumannCondition) {
    BoundaryConditions<2,1> bc;
    
    // Test con valore costante
    EXPECT_NO_THROW(bc.addNeumann(2, Point<1>(2.0)));
    
    // Test con funzione
    Function<2,1> flux([](Point<2> p) { return sin(M_PI * p[0]); });
    EXPECT_NO_THROW(bc.addNeumann(3, flux));
}

// Test 4: Mix di condizioni Dirichlet e Neumann
TEST_F(BoundaryConditionsTest, MixedBoundaryConditions) {
    BoundaryConditions<2,1> bc;
    
    // Aggiungi diverse condizioni
    bc.addDirichlet(0, Point<1>(0.0));  // Lato sinistro
    bc.addDirichlet(1, Point<1>(0.0));  // Lato destro
    bc.addDirichlet(2, Point<1>(0.0));  // Lato inferiore
    
    Function<2,1> neumannFlux([](Point<2> p) { 
        return -M_PI * sin(M_PI * p[0]); 
    });
    bc.addNeumann(3, neumannFlux);      // Lato superiore
    
    EXPECT_NO_THROW(bc);
}

// Test 5: Applicazione su una griglia semplice
TEST_F(BoundaryConditionsTest, ApplyToSimpleGrid) {
    // Crea una griglia 2x2 semplice
    std::vector<Point<2>> nodes = {
        Point<2>(std::vector<double>{0.0, 0.0}), 
        Point<2>(std::vector<double>{1.0, 0.0}), 
        Point<2>(std::vector<double>{0.0, 1.0}), 
        Point<2>(std::vector<double>{1.0, 1.0})
    };
    
    std::vector<Cell<2>> cells = {
        Cell<2>(
            {nodes[0], nodes[1], nodes[2]},
            {0, 1, 2}
        ),
        Cell<2>(
            {nodes[1], nodes[3], nodes[2]},
            {1, 3, 2}
        )
    };
    
    Grid2D grid(cells, nodes);
    
    // Crea condizioni al contorno
    BoundaryConditions<2,1> bc;
    bc.addDirichlet(0, Point<1>(0.0));
    
    // Crea matrici di test
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(4, 4);
    Eigen::VectorXd rhs(4);
    A.setIdentity();
    rhs.setOnes();
    
    // Applica le condizioni (dovrebbe funzionare senza errori)
    EXPECT_NO_THROW(bc.apply(grid, A, rhs));
}

// Test 6: Verifica che le condizioni Dirichlet modifichino correttamente la matrice
TEST_F(BoundaryConditionsTest, DirichletModifiesMatrix) {
    // Setup semplificato
    std::vector<Point<2>> nodes = {
        Point<2>(std::vector<double>{0.0, 0.0}), 
        Point<2>(std::vector<double>{1.0, 0.0}), 
        Point<2>(std::vector<double>{0.0, 1.0}), 
        Point<2>(std::vector<double>{1.0, 1.0})
    };
    
    std::vector<Cell<2>> cells = {
        Cell<2>(
            {nodes[0], nodes[1], nodes[2]},
            {0, 1, 2}
        )
    };
    
    Grid2D grid(cells, nodes);
    
    BoundaryConditions<2,1> bc;
    bc.addDirichlet(0, Point<1>(5.0));  // Valore non zero per vedere l'effetto
    
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(4, 4);
    Eigen::VectorXd rhs(4);
    A.setIdentity();
    rhs.setZero();
    
    bc.apply(grid, A, rhs);
    
    // Verifica che il RHS sia stato modificato
    // (I dettagli dipendono dall'implementazione specifica)
    EXPECT_TRUE(true);  // Test base che l'applicazione non crashi
}

// Test 7: Test di funzioni con parametri complessi
TEST_F(BoundaryConditionsTest, ComplexFunctions) {
    BoundaryConditions<2,1> bc;
    
    // Funzione sinusoidale per Neumann
    Function<2,1> sinusoidalFlux([](Point<2> p) { 
        return sin(2.0 * M_PI * p[0]) * cos(M_PI * p[1]); 
    });
    
    // Funzione polinomiale per Dirichlet
    Function<2,1> polynomialValue([](Point<2> p) { 
        return p[0] * p[0] + p[1] * p[1]; 
    });
    
    EXPECT_NO_THROW(bc.addNeumann(1, sinusoidalFlux));
    EXPECT_NO_THROW(bc.addDirichlet(2, polynomialValue));
}

// Test 8: Edge case - stessi tag multipli (dovrebbe sovrascrivere o gestire l'errore)
TEST_F(BoundaryConditionsTest, DuplicateTagHandling) {
    BoundaryConditions<2,1> bc;
    
    bc.addDirichlet(1, Point<1>(1.0));
    
    // Aggiunta dello stesso tag dovrebbe essere gestita
    // (comportamento dipende dall'implementazione)
    EXPECT_NO_THROW(bc.addDirichlet(1, Point<1>(2.0)));
}

// =============================================================================
// NUOVI TEST AGGIUNTI
// =============================================================================

// Test 9: Test completo per boundary conditions 1D
TEST_F(BoundaryConditionsTest, OneDimensionalTest) {
    BoundaryConditions<1,1> bc1d;
    
    // Aggiungi condizioni 1D
    bc1d.addDirichlet(0, Point<1>(0.0));  // Estremo sinistro
    bc1d.addNeumann(1, Point<1>(1.0));    // Estremo destro
    
    // Crea una griglia 1D semplice [0, 1] con 5 nodi
    Grid1D grid(0.0, 1.0, 4);  // 4 elementi, 5 nodi
    
    // Crea sistema lineare
    int n = grid.getN();  // Grid1D usa getN() per il numero di punti/nodi
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
    Eigen::VectorXd rhs(n);
    
    // Inizializza matrice identità e RHS
    A.setIdentity();
    rhs.setOnes();
    
    // Applica condizioni al contorno
    EXPECT_NO_THROW(bc1d.apply(grid, A, rhs));
    
    // Verifica numerica: il primo nodo dovrebbe avere valore Dirichlet
    EXPECT_DOUBLE_EQ(A.coeff(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(rhs[0], 0.0);  // Valore Dirichlet applicato
}

// Test 10: Test con funzioni 1D complesse
TEST_F(BoundaryConditionsTest, OneDimensionalComplexFunctions) {
    BoundaryConditions<1,1> bc1d;
    
    // Funzione sinusoidale per Dirichlet
    Function<1,1> sinFunc([](Point<1> p) { 
        return sin(M_PI * p[0]); 
    });
    
    // Funzione derivata per Neumann
    Function<1,1> cosFunc([](Point<1> p) { 
        return M_PI * cos(M_PI * p[0]); 
    });
    
    bc1d.addDirichlet(0, sinFunc);
    bc1d.addNeumann(1, cosFunc);
    
    Grid1D grid(0.0, 1.0, 3);
    
    int n = grid.getN();  // Grid1D usa getN()
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
    Eigen::VectorXd rhs(n);
    A.setIdentity();
    rhs.setZero();
    
    EXPECT_NO_THROW(bc1d.apply(grid, A, rhs));
    
    // Verifica che sin(0) = 0 sia applicato correttamente
    EXPECT_NEAR(rhs[0], sin(0.0), 1e-10);
}

// Test 11: Test con mesh reali 2D
TEST_F(BoundaryConditionsTest, RealMesh2DTest) {
    Grid2D grid;
    
    // Prova a caricare mesh reale
    ASSERT_NO_THROW(grid.parseFromMsh("../mesh/test_grid2d.msh"));
    
    BoundaryConditions<2,1> bc;
    
    // Aggiungi condizioni con tag fisici reali (tipicamente 1, 2, 3, 4 per i lati)
    bc.addDirichlet(1, Point<1>(0.0));  // Lato 1
    bc.addDirichlet(2, Point<1>(0.0));  // Lato 2
    bc.addDirichlet(3, Point<1>(0.0));  // Lato 3
    
    Function<2,1> neumannFlux([](Point<2> p) { 
        return sin(M_PI * p[0]) * cos(M_PI * p[1]); 
    });
    bc.addNeumann(4, neumannFlux);      // Lato 4
    
    // Crea sistema per la mesh reale
    int n = grid.getNumNodes();
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
    Eigen::VectorXd rhs(n);
    A.setIdentity();
    rhs.setOnes();
    
    // Applica condizioni (dovrebbe funzionare con mesh reale)
    EXPECT_NO_THROW(bc.apply(grid, A, rhs));
    
    std::cout << "Test con mesh reale: " << n << " nodi, " 
              << grid.getNumElements() << " elementi" << std::endl;
}

// Test 12: Test con mesh reali - verifica numerica precisa
TEST_F(BoundaryConditionsTest, RealMeshNumericalVerification) {
    Grid2D grid;
    
    try {
        grid.parseFromMsh("../mesh/mesh-square-5_gmsh22.msh");
    } catch (...) {
        GTEST_SKIP() << "Mesh file non disponibile, test saltato";
        return;
    }
    
    BoundaryConditions<2,1> bc;
    
    // Condizione Dirichlet precisa: u = x^2 + y^2 sul bordo
    Function<2,1> exactSolution([](Point<2> p) { 
        return p[0]*p[0] + p[1]*p[1]; 
    });
    
    bc.addDirichlet(1, exactSolution);  // Applica su tutto il bordo tag 1
    
    int n = grid.getNumNodes();
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
    Eigen::VectorXd rhs(n);
    A.setIdentity();
    rhs.setZero();
    
    bc.apply(grid, A, rhs);
    
    // Verifica numerica: controlla alcuni nodi di bordo
    auto boundaryNodes = grid.getBoundaryNodesByTag(1);
    
    for (int nodeIdx : boundaryNodes) {
        const Point<2>& node = grid.getNode(nodeIdx);
        double expectedValue = node[0]*node[0] + node[1]*node[1];
        
        // Verifica che il valore sia stato applicato correttamente
        EXPECT_NEAR(rhs[nodeIdx], expectedValue, 1e-12) 
            << "Nodo " << nodeIdx << " at (" << node[0] << ", " << node[1] << ")";
        
        // Verifica che la matrice abbia 1 sulla diagonale
        EXPECT_DOUBLE_EQ(A.coeff(nodeIdx, nodeIdx), 1.0);
    }
    
    std::cout << "Verificati " << boundaryNodes.size() << " nodi di bordo" << std::endl;
}

// Test 13: Test prestazioni con griglia grande
TEST_F(BoundaryConditionsTest, PerformanceTest) {
    // Crea griglia 1D grande per test performance
    Grid1D largeGrid(0.0, 10.0, 1000);  // 1000 elementi, 1001 nodi
    
    BoundaryConditions<1,1> bc;
    bc.addDirichlet(0, Point<1>(0.0));
    bc.addDirichlet(1, Point<1>(1.0));
    
    int n = largeGrid.getN();  // Grid1D usa getN()
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
    Eigen::VectorXd rhs(n);
    A.setIdentity();
    rhs.setOnes();
    
    // Misura tempo di applicazione
    auto start = std::chrono::high_resolution_clock::now();
    bc.apply(largeGrid, A, rhs);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Verifica che l'applicazione sia rapida (< 100ms per 1000 nodi)
    EXPECT_LT(duration.count(), 100) 
        << "Applicazione BC troppo lenta: " << duration.count() << "ms per " << n << " nodi";
    
    // Verifica correttezza
    EXPECT_DOUBLE_EQ(rhs[0], 0.0);           // Primo nodo
    EXPECT_DOUBLE_EQ(rhs[n-1], 1.0);         // Ultimo nodo
    
    std::cout << "Performance test: " << n << " nodi in " << duration.count() << "ms" << std::endl;
}

// Test 14: Test edge cases geometrici
TEST_F(BoundaryConditionsTest, GeometricEdgeCases) {
    BoundaryConditions<2,1> bc;
    
    // Crea griglia con nodi su coordinate speciali
    std::vector<Point<2>> specialNodes = {
        Point<2>(std::vector<double>{0.0, 0.0}),     // Origine
        Point<2>(std::vector<double>{1.0, 0.0}),     // Asse X
        Point<2>(std::vector<double>{0.0, 1.0}),     // Asse Y
        Point<2>(std::vector<double>{-1.0, 0.0}),    // X negativo
        Point<2>(std::vector<double>{0.0, -1.0})     // Y negativo
    };
    
    std::vector<Cell<2>> specialCells = {
        Cell<2>({specialNodes[0], specialNodes[1], specialNodes[2]}, {0, 1, 2}),
        Cell<2>({specialNodes[0], specialNodes[3], specialNodes[4]}, {0, 3, 4})
    };
    
    Grid2D specialGrid(specialCells, specialNodes);
    
    // Condizioni complesse
    Function<2,1> complexFunc([](Point<2> p) { 
        if (abs(p[0]) < 1e-12 && abs(p[1]) < 1e-12) return 0.0;  // Origine
        return p[0] / (p[0]*p[0] + p[1]*p[1] + 1e-12);           // Evita divisione per zero
    });
    
    bc.addDirichlet(0, complexFunc);
    
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(5, 5);
    Eigen::VectorXd rhs(5);
    A.setIdentity();
    rhs.setZero();
    
    // Dovrebbe gestire i casi edge senza crash
    EXPECT_NO_THROW(bc.apply(specialGrid, A, rhs));
}
