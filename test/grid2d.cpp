#include <gtest/gtest.h>
#include "grid2d.hpp"
#include <cmath>
#include <vector>

class Grid2DTest : public ::testing::Test {
};

TEST_F(Grid2DTest, Parse_test) {
    Grid2D grid;
    grid.parseFromMsh("../mesh/test_grid2d.msh");
    return GTEST_SUCCEED(); 
}
TEST_F(Grid2DTest, Parse_reprint_test) {
    Grid2D grid;
    grid.parseFromMsh("../mesh/test_grid2d.msh");

    std::ofstream outFile("../test/test_grid2d_out.msh");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        exit(-1);
    }
    outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    outFile << "$Nodes\n" << grid.getNumNodes() << "\n";
    for (int i = 0; i < grid.getNumElements(); ++i) {
        const Cell<2>& cell = grid.getCell(i);
        for (int j = 0; j < cell.getN(); ++j) {
            const Point<2>& node = cell[j];
            outFile << i + 1 << " " << node.coords[0] << " " << node.coords[1] << " 0.0\n";
        }
    }
    outFile << "$EndNodes\n";

    outFile << "$Elements\n" << grid.getNumElements() << "\n";
    for (int i = 0; i < grid.getNumElements(); ++i) {
        const Cell<2>& cell = grid.getCell(i);
        outFile << i + 1 << " 2 2 0 0 ---- UNFINISHED ----\n";
        // << cell.getType() << " " << cell.getNodeIndices() << "\n";
    }
    outFile << "$EndElements\n";

    outFile.close();
    std::cout << "Reprinted grid to ../test/test_grid2d_out.msh" << std::endl; 

    return GTEST_SUCCEED(); 
}
