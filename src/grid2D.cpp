#include "grid2D.hpp"

template void Grid<2>::parseFromMsh(const std::string&);
template void Grid<3>::parseFromMsh(const std::string&);

bool compare(std::string line, std::string text) {
    line.erase(line.find_last_not_of(" \r\n") + 1);
    return line.compare(text) == 0;
}

// Assume format .msh 2.2
template<unsigned int dim>
void Grid<dim>::parseFromMsh(const std::string& filename) {
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
        std::vector<double> coords(3, 0.0);
        // Parse line of format index x y z into variables
        if (line.empty()) continue; // Accept empty lines
        if (sscanf(line.c_str(), "%d %lf %lf %lf", &index, &coords[0], &coords[1], &coords[2]) != 4) {
            std::cerr << "Error parsing node line: " << line << std::endl;
            exit(-1);
        }
        coords.resize(dim);
        unique_nodes.emplace_back(Point<dim>(coords));
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
        
        // Read physical tag (first tag) and skip the rest
        elemPhysTag = 0; // Default value
        for (int i = 0; i < numTags; i++) {
            int tag;
            iss >> tag;
            if (i == 0) elemPhysTag = tag; // Save physical tag
        }

        if (isBoundaryCell(elemType)) {
            parseBoundaryCell(iss, elemPhysTag);
        } else {
            parseInternalCell(iss);
        }
    }

    file.close();
}

template<>
bool Grid<2>::isBoundaryCell(int elemType) {
    if (elemType == 1) return true;
    if (elemType == 2) return false; // Correct volume element type
    else {
        std::cerr << "Unsupported element type: " << elemType << std::endl;
        exit(-1);
    }
}
template<>
bool Grid<3>::isBoundaryCell(int elemType) {
    if (elemType == 2) return true;
    if (elemType == 4) return false; // Correct volume element type
    else {
        std::cerr << "Unsupported element type: " << elemType << std::endl;
        exit(-1);
    }
}

template<unsigned int dim>
void Grid<dim>::parseBoundaryCell(std::istringstream& iss, int boundaryId){
    typename BoundaryCell<dim-1>::NodeVector elementNodes;
    typename BoundaryCell<dim-1>::NodeIndexes elementNodeIndices;

    unsigned int nElem = dim;
    // 1D mesh -> 0D boundary -> point -> 1 element
    // 2D mesh -> 1D boundary -> line -> 2 elements
    // 3D mesh -> 2D boundary -> triangle -> 3 elements

    for (int i = 0; i < nElem; ++i) {
        unsigned int nodeIndex;
        iss >> nodeIndex;

        if (nodeIndex < 1 || nodeIndex > unique_nodes.size()) {
            std::cerr << "Boundary node index out of bounds: " << nodeIndex
                        << " (max: " << unique_nodes.size() << ")" << std::endl;
            exit(-1); // Skip this element
        }

        elementNodes.push_back(unique_nodes[nodeIndex - 1]);
        elementNodeIndices.push_back(nodeIndex - 1);
    }

    boundary_cells.emplace_back(BoundaryCell<dim-1>(elementNodes, elementNodeIndices, boundaryId));
}

template<unsigned int dim>
void Grid<dim>::parseInternalCell(std::istringstream& iss) {
    typename Cell<dim>::NodeVector elementNodes;
    typename Cell<dim>::NodeIndexes elementNodeIndices;

    unsigned int nElem = dim+1;
    // 1D mesh -> 1D element -> line -> 2 elements
    // 2D mesh -> 2D element -> triangle -> 3 elements
    // 3D mesh -> 3D element -> tetrahedron -> 4 elements

    unsigned int nodeIndex;
    for (int i = 0; i < nElem; ++i) {
        iss >> nodeIndex;

        if (nodeIndex < 1 || nodeIndex > unique_nodes.size()) {
            std::cerr << "Triangle node index out of bounds: " << nodeIndex
                        << " (max: " << unique_nodes.size() << ")" << std::endl;
            exit(-1);
        }

        elementNodes.push_back(unique_nodes[nodeIndex - 1]);
        elementNodeIndices.push_back(nodeIndex - 1);
    }

    cells.emplace_back(Cell<dim>(elementNodes, elementNodeIndices));
}
