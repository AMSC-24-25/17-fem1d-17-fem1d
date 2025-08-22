#ifndef BCTYPE_HPP
#define BCTYPE_HPP

// Enum for the type of boundary condition
enum class BCType {
    DIRICHLET,  // u = specified value
    NEUMANN     // du/dn = specified value
};

#endif // BCTYPE_HPP