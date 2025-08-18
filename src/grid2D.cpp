#include "grid2D.hpp"

bool compare(std::string line, std::string text) {
    line.erase(line.find_last_not_of(" \r\n") + 1);
    return line.compare(text) == 0;
}

// Assume format .msh 2.2
void Grid2D::parseFromMsh(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(-1);
    }
    
    int index, dump;
    std::string line;
    
    // Skip until line $Nodes
    while (std::getline(file, line) && !compare(line, "$Nodes"));

    // Read number of nodes
    std::getline(file, line); 
    int numNodes = atoi(line.c_str());
    unique_nodes.clear();
    unique_nodes.reserve(numNodes);

    // Read nodes
    while (std::getline(file, line) && !compare(line, "$EndNodes")) {
        double x,y,z;
        // Parse line of format index x y z into variables
        if (line.empty()) continue; // Accept empty lines
        if (sscanf(line.c_str(), "%d %lf %lf %lf", &index, &x, &y, &z) != 4) {
            std::cerr << "Error parsing node line: " << line << std::endl;
            exit(-1);
        }
        unique_nodes.emplace_back(Point<2>(x, y));
    }

    // Skip $Elements line
    std::getline(file, line);
    if (!compare(line, "$Elements")) {
        std::cerr << "Expected $Elements, found: " << line << std::endl;
        exit(-1);
    }

    // Read number of elements
    std::getline(file, line);
    int numElements = atoi(line.c_str());
    cells.clear();
    cells.reserve(numElements);

    int elemType, elemPhysTag, elemGeomTag;
    // Read elements
    while (std::getline(file, line) && !compare(line, "$EndElements")) {
        char elemLine[256];
        if (sscanf(line.c_str(), "%d %d %d %d %d %255s", &index, &elemType, &dump, &elemPhysTag, &elemGeomTag, elemLine) != 6) {
            std::cerr << "Error parsing element line: " << line << std::endl;
            exit(-1);
        }

        if(elemType == 1){ // Line element == boundary
            BoundaryCell<1>::NodeVector elementNodes;
            BoundaryCell<1>::NodeIndexes elementNodeIndices;

            std::istringstream iss(elemLine);
            int nodeIdx;
            while (iss >> nodeIdx) {
                if (nodeIdx < 1 || nodeIdx > unique_nodes.size()) {
                    std::cerr << "Node index out of bounds: " << nodeIdx << std::endl;
                    exit(-1);
                }
                elementNodes.push_back(unique_nodes[nodeIdx - 1]);
                elementNodeIndices.push_back(nodeIdx - 1);
            }

            boundary_cells.emplace_back(BoundaryCell<1>(elementNodes, elementNodeIndices, elemPhysTag));
        } 
        else{ // 2D element == internal
            Cell<2>::NodeVector elementNodes;
            Cell<2>::NodeIndexes elementNodeIndices;
            std::istringstream iss(elemLine);
            int nodeIdx;
            while (iss >> nodeIdx) {
                if (nodeIdx < 1 || nodeIdx > unique_nodes.size()) {
                    std::cerr << "Node index out of bounds: " << nodeIdx << std::endl;
                    exit(-1);
                }
                elementNodes.push_back(unique_nodes[nodeIdx - 1]);
                elementNodeIndices.push_back(nodeIdx - 1);
            }
            
            cells.emplace_back(Cell<2>(elementNodes, elementNodeIndices));
        }
    }

    file.close();
}

Grid2D::PhiFunctionVector2D Grid2D::getPhiFunctions() const {
    PhiFunctionVector2D phiFunctions;
    
    // Per ogni cella, crea le shape functions per i suoi nodi
    for (unsigned int cellIdx = 0; cellIdx < cells.size(); ++cellIdx) {
        const Cell<2>& cell = cells[cellIdx];
        for (unsigned int nodeIdx = 0; nodeIdx < cell.getN(); ++nodeIdx) {
            phiFunctions.push_back(PhiFunction2D(nodeIdx, cell));
        }
    }
    
    return phiFunctions;
}
Grid2D::PhiFunctionVector2D Grid2D::getBoundaryPhiFunctions() const {
    PhiFunctionVector2D phiFunctions;
    
    // Per ogni cella, crea le shape functions per i suoi nodi
    for (unsigned int cellIdx = 0; cellIdx < boundary_cells.size(); ++cellIdx) {
        const BoundaryCell<1>& boundary_cell = boundary_cells[cellIdx];
        for (unsigned int nodeIdx = 0; nodeIdx < boundary_cell.getN(); ++nodeIdx) {
            phiFunctions.push_back(BoundaryPhiFunction2D(nodeIdx, boundary_cell));
        }
    }
    
    return phiFunctions;
}