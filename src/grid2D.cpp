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
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        int numTags;
        iss >> index >> elemType >> numTags;
        
        // Skip all tags
        for (int i = 0; i < numTags; i++) {
            int tag;
            iss >> tag;
        }

        if(elemType == 1){ // Line element == boundary
            BoundaryCell<1>::NodeVector elementNodes;
            BoundaryCell<1>::NodeIndexes elementNodeIndices;

            std::vector<int> nodeIds(2);
            for (int i = 0; i < 2; ++i) {
                iss >> nodeIds[i];
                
                if (nodeIds[i] < 1 || nodeIds[i] > unique_nodes.size()) {
                    std::cerr << "Boundary node index out of bounds: " << nodeIds[i] 
                              << " (max: " << unique_nodes.size() << ")" << std::endl;
                    continue; // Skip this element
                }

                elementNodes.push_back(unique_nodes[nodeIds[i] - 1]);
                elementNodeIndices.push_back(nodeIds[i] - 1);
            }

            boundary_cells.emplace_back(BoundaryCell<1>(elementNodes, elementNodeIndices, elemPhysTag));
        }else if(elemType == 2) { // Triangle element == internal
            Cell<2>::NodeVector elementNodes;
            Cell<2>::NodeIndexes elementNodeIndices;

            std::vector<int> nodeIds(3);
            for (int i = 0; i < 3; ++i) {
                iss >> nodeIds[i];

                if (nodeIds[i] < 1 || nodeIds[i] > unique_nodes.size()) {
                    std::cerr << "Triangle node index out of bounds: " << nodeIds[i] 
                              << " (max: " << unique_nodes.size() << ")" << std::endl;
                    exit(-1);
                }

                elementNodes.push_back(unique_nodes[nodeIds[i] - 1]);
                elementNodeIndices.push_back(nodeIds[i] - 1);
            }

            cells.emplace_back(Cell<2>(elementNodes, elementNodeIndices));
        } else {
            // Unsupported element type, we require triangular elements
            std::cerr << "Unsupported element type: " << elemType << std::endl;
            exit(-1);
        }
    }

    file.close();
}