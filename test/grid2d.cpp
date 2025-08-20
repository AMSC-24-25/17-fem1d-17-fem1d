#include <gtest/gtest.h>
#include "grid2D.hpp"
#include <cmath>
#include <vector>

class GridTest : public ::testing::Test {
};

TEST_F(GridTest, Parse_2D_test) {
    Grid2D grid;
    grid.parseFromMsh("../mesh/test_grid2d.msh");
    return GTEST_SUCCEED(); 
}
TEST_F(GridTest, Parse_3D_test) {
    Grid3D grid;
    grid.parseFromMsh("../mesh/test_grid3d.msh");
    return GTEST_SUCCEED(); 
}
TEST_F(GridTest, Parse_and_output_2d_test) {
    Grid2D grid;
    grid.parseFromMsh("../mesh/test_grid2d.msh");

    std::ofstream outFile("../test/test_grid2d_out.msh");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        exit(-1);
    }
    outFile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    outFile << "$Nodes\n" << grid.getNumNodes() << "\n";
    for (int i = 0; i < grid.getNumNodes(); ++i) {
        const Point<2>& node = grid.getNode(i);
        outFile << i + 1 << " " << node.coords[0] << " " << node.coords[1] << " 0.0\n";
    }
    outFile << "$EndNodes\n";

    outFile << "$Elements\n" << grid.getNumElements() << "\n";
    for (int i = 0; i < grid.getNumElements(); ++i) {
        const Cell<2>& cell = grid.getCell(i);
        outFile << i + 1 << " 2 2 0 0 ";
        for (int j = 0; j < cell.getN(); ++j) {
            outFile << cell.getNodeIndex(j) << " ";
        }
        outFile << "\n";
    }
    outFile << "$EndElements\n";

    outFile.close();
    std::cout << "Reprinted grid to ../test/test_grid2d_out.msh" << std::endl; 

    return GTEST_SUCCEED(); 
}
