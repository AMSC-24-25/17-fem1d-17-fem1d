# Fem1d
1D Finite Element

## What is it?

The Finite Element Method (FEM) is a numerical method used to solve partial differential equations in problems ranging from one-dimensional to multi-dimensional. 

FEM works by dividing the system into small parts, called finite elements, and discretizing the space using a mesh. This converts a boundary value problem into a system of algebraic equations. The method estimates the unknown function over the domain by assembling equations for the finite elements into a larger system that represents the entire problem. FEM then finds an approximate solution by minimizing an associated error function using the calculus of variations.


## How does it work?

We implemented a program capable of solving partial differential equations problems in one dimension.

All the functions of the program are implemented using a specif class characterized by two attributes: one for the function itself and one for its gradient (defined as std funcitons). Then there are two constructors: respectively the first need both the attributes while the second one need only the function because the gradient is 


The most important class is fem1d. It is responsible of handling all the data of the problem, building the algebraic linear system associated to the analytical PDE through finite element method. This is done through
<ul>
    <li>Constructor: takes already instantiated data for the problem. Here all state variables of the class are instantiated.</li>
    <li>Assemble method: actually computes 
</ul>



m MEFF hguorhtnisu  lacitylana EDP eht o

## What to do next?

