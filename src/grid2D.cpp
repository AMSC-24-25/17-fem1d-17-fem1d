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

    // Read elements
    while (std::getline(file, line) && !compare(line, "$EndElements")) {
        Cell<2>::NodeVector elementNodes;
        Cell<2>::NodeIndexes elementNodeIndices;
        char elemLine[256];
        if (sscanf(line.c_str(), "%d %d %d %d %d %255s", &index, &dump, &dump, &dump, &dump, elemLine) != 6) {
            std::cerr << "Error parsing element line: " << line << std::endl;
            exit(-1);
        }
        
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

    file.close();
}