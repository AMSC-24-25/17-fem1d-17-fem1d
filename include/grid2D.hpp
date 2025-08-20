#ifndef GRID_2D
#define GRID_2D

#include "cell.hpp"
#include "point.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>

class Grid2D{
    using CellVector = std::vector<Cell<2>>;
    using BoundaryCellVector = std::vector<BoundaryCell<1>>;
    using NodeVector = Cell<2>::NodeVector;
    using BoundaryNodeVector = Cell<1>::NodeVector;
    using IndexVector = Cell<2>::NodeIndexes;

    CellVector cells;
    BoundaryCellVector boundary_cells;
    NodeVector unique_nodes;

public:
    Grid2D(CellVector cells, NodeVector unique_nodes) : cells(cells), unique_nodes(unique_nodes)
    {}

    Grid2D() : cells(), unique_nodes() {}

    void addCell(const Cell<2>& cell) {
        cells.push_back(cell);
    }
    unsigned int getNumElements() const {
        return cells.size();
    }
    unsigned int getNumNodes() const {
        return unique_nodes.size();
    }
    const Cell<2>& getCell(unsigned int i) const {
        if (i >= cells.size()) {
            std::cerr << "Index out of bounds in Grid2D::getCell\n";
            exit(-1);
        }
        return cells[i];
    }
    const CellVector& getCells() const {
        return cells;
    }
    const Point<2>& getNode(unsigned int i) const {
        if (i >= unique_nodes.size()) {
            std::cerr << "Index out of bounds in Grid2D::getNode\n";
            exit(-1);
        }
        return unique_nodes[i];
    }
    const NodeVector& getUniqueNodes() const {
        return unique_nodes;
    }

    void parseFromMsh(const std::string& filename);

    // Boundary cells access methods
    unsigned int getNumBoundaryCells() const {
        return boundary_cells.size();
    }

    const BoundaryCell<1>& getBoundaryCell(unsigned int i) const {
        if (i >= boundary_cells.size()) {
            std::cerr << "Index out of bounds in Grid2D::getBoundaryCell\n";
            exit(-1);
        }
        return boundary_cells[i];
    }

    const BoundaryCellVector& getBoundaryCells() const {
        return boundary_cells;
    }

    // Boundary nodes extraction methods
    IndexVector getBoundaryNodes() const {
        IndexVector boundaryNodes;
        std::set<unsigned int> uniqueNodes; // Use set to avoid duplicates
        
        for (const BoundaryCell<1>& boundaryCell : boundary_cells) {
            const IndexVector& nodeIndices = boundaryCell.getNodeIndexes();
            uniqueNodes.insert(nodeIndices.begin(), nodeIndices.end());
        }
        
        // Convert set to vector
        boundaryNodes.assign(uniqueNodes.begin(), uniqueNodes.end());
        return boundaryNodes;
    }

    IndexVector getBoundaryNodesByTag(int physicalTag) const {
        IndexVector boundaryNodes;
        std::set<unsigned int> uniqueNodes; // Use set to avoid duplicates

        for (const BoundaryCell<1>& boundaryCell : boundary_cells) {
            if (boundaryCell.getBoundaryId() == physicalTag) {
                const IndexVector& nodeIndices = boundaryCell.getNodeIndexes();
                uniqueNodes.insert(nodeIndices.begin(), nodeIndices.end());
            }
        }
        
        // Convert set to vector
        boundaryNodes.assign(uniqueNodes.begin(), uniqueNodes.end());
        return boundaryNodes;
    }

    std::vector<BoundaryCell<1>> getBoundaryEdgesByTag(int physicalTag) const {
        std::vector<BoundaryCell<1>> boundaryEdges;
        
        for (const auto& boundaryCell : boundary_cells) {
            if (boundaryCell.getPhysicalTag() == physicalTag) {
                boundaryEdges.push_back(boundaryCell);
            }
        }
        
        return boundaryEdges;
    }

    bool isBoundaryNode(int nodeIndex) const {
        for (const auto& boundaryCell : boundary_cells) {
            const auto& nodeIndices = boundaryCell.getNodeIndexes();
            for (int idx : nodeIndices) {
                if (idx == nodeIndex) {
                    return true;
                }
            }
        }
        return false;
    }

    // Physical tags utility methods
    IndexVector getPhysicalTags() const {
        std::set<int> uniqueTags;

        for (const BoundaryCell<1>& boundaryCell : boundary_cells) {
            uniqueTags.insert(boundaryCell.getBoundaryId());
        }
        
        // Convert set to vector
        IndexVector tags(uniqueTags.begin(), uniqueTags.end());
        return tags;
    }
};

#endif
