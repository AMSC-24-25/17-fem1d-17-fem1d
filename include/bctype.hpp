/**
 * @file bctype.hpp
 * @brief Boundary condition type enumeration
 */
#ifndef BCTYPE_HPP
#define BCTYPE_HPP

/**
 * @brief Types of boundary conditions for PDEs
 */
enum class BCType {
    DIRICHLET,  // u = specified value
    NEUMANN     // du/dn = specified value
};

#endif // BCTYPE_HPP