#include "grid2D.hpp"

// Assume format .msh 2.2
void Grid2D::parseFromMsh(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    int dump;
    std::string line;
    
    // Skip until line $Nodes
    while (std::getline(file, line) && !line.compare("$Nodes"));
    
    NodeVector nodes;
    // Read number of nodes
    std::getline(file, line); 
    int numNodes = atoi(line.c_str());
    nodes.reserve(numNodes);
    // Read nodes
    while (std::getline(file, line) && !line.compare("$EndNodes")) {
        double x,y,z;
        // Parse line of format index x y z into variables
        if (line.empty()) continue; // Accept empty lines
        if (sscanf(line.c_str(), "%d %lf %lf %lf", &dump, &x, &y, &z) != 4) {
            std::cerr << "Error parsing node line: " << line << std::endl;
            exit(-1);
        }
        nodes.emplace_back(Point<2>(x, y));
    }

    std::getline(file, line);
    if (line != "$Elements") {
        std::cerr << "Expected $Elements, found: " << line << std::endl;
        exit(-1);
    }
    std::getline(file, line);
    int numElements = atoi(line.c_str());
    cells.clear();
    cells.reserve(numElements);
    while (std::getline(file, line) && !line.compare("$EndElements")) {
        Cell<2>::NodeVector elementNodes;
        std::string elemLine;
        if (sscanf(line.c_str(), "%d %d %d %d %d %s", &dump, &dump, &dump, &dump, &dump, &elemLine) != 6) {
            std::cerr << "Error parsing element line: " << line << std::endl;
            exit(-1);
        }
        
        std::istringstream iss(elemLine);
        int nodeIdx;
        while (iss >> nodeIdx) {
            if (nodeIdx < 1 || nodeIdx > nodes.size()) {
                std::cerr << "Node index out of bounds: " << nodeIdx << std::endl;
                exit(-1);
            }
            elementNodes.push_back(nodes[nodeIdx - 1]);
        }

        cells.emplace_back(Cell<2>(elementNodes));
    }

    file.close();
}