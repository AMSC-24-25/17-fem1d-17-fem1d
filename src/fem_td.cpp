#include "fem_td.hpp"
#include <fstream>
#include <iostream>
#include <cmath>


// Time-dependent FEM constructor: initialize matrices and parallel structures
template<unsigned int dim>
FemTD<dim>::FemTD(Grid<dim> grid,
                  Function<dim,1> diffusion,
                  Function<dim,dim> transport,
                  Function<dim,1> reaction,
                  const BoundaryConditions_td<dim,1>& bc,
                  QuadratureRule<dim> quadrature)
: mesh_(std::move(grid))
, diffusion_(std::move(diffusion))
, transport_(std::move(transport))
, reaction_(std::move(reaction))
, bc_(bc)
, quad_(std::move(quadrature))
, forcing_td_([](const Point<dim>&, double){ return 0.0; })
, u0_([](const Point<dim>&){ return 0.0; })
{
#ifdef _OPENMP
    // Ensure Eigen uses OpenMP if available
    nthreads = omp_get_max_threads();
    Eigen::setNbThreads(nthreads);
#endif

    const int N = mesh_.getNumNodes();
    M_.resize(N,N);
    K_.resize(N,N);
    A_.resize(N,N);
    rhs_.resize(N); rhs_.setZero();
    u_.resize(N);   u_.setZero();
    u_old_.resize(N); u_old_.setZero();
    f_new.resize(N); f_new.setZero();
    f_old.resize(N); f_old.setZero();

    // The following are preallocated for improved memory efficiency

    const unsigned int expPerElem = (dim + 1) * (dim + 1);
    const unsigned int numElements = mesh_.getNumElements();
    tripletM.reserve(expPerElem * numElements);
    tripletK.reserve(expPerElem * numElements);
#ifdef _OPENMP
    F_threads.resize(nthreads, VectorXd::Zero(N));
    tripletM_thr.resize(nthreads);
    tripletK_thr.resize(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        tripletM_thr[i].reserve(expPerElem * (numElements/nthreads + 1));
        tripletK_thr[i].reserve(expPerElem * (numElements/nthreads + 1));
    }
#endif
}

template<unsigned int dim>
void FemTD<dim>::set_forcing(ForcingTD f_td) { forcing_td_ = std::move(f_td); }

template<unsigned int dim>
void FemTD<dim>::set_initial_condition(Function<dim,1> u0) { u0_ = std::move(u0); }


// Assemble time-invariant matrices M (mass) and K (stiffness) with parallel assembly
template<unsigned int dim>
void FemTD<dim>::assemble_time_invariant() {
    const int numElements = mesh_.getNumElements();
    const unsigned matSize = dim + 1;
    tripletM.clear();
    tripletK.clear();

#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<Triplet>& tM = tripletM_thr[tid];
        std::vector<Triplet>& tK = tripletK_thr[tid];
        tM.clear();
        tK.clear();

        #pragma omp for schedule(dynamic, 32)
        for (int e=0; e < numElements; ++e) assemble_M_and_K_element(e, tM, tK);
    }
    for (std::vector<Triplet>& v: tripletM_thr) tripletM.insert(tripletM.end(), v.begin(), v.end());
    for (std::vector<Triplet>& v: tripletK_thr) tripletK.insert(tripletK.end(), v.begin(), v.end());
#else
    for (int e=0; e<numElements; ++e) assemble_M_and_K_element(e, tripletM, tripletK);
#endif

    M_.setFromTriplets(tripletM.begin(), tripletM.end());
    K_.setFromTriplets(tripletK.begin(), tripletK.end());
}


// Element assembly: compute local mass M and stiffness K matrices for time-dependent PDE
template<unsigned int dim>
void FemTD<dim>::assemble_M_and_K_element(int elem,
                                          std::vector<Triplet>& tripletM,
                                          std::vector<Triplet>& tripletK)
{
    const Cell<dim>& cell = mesh_.getCell(elem);
    const unsigned matSize = dim + 1;

    std::vector<Point<dim>> grad_phi;
    std::vector<Point<dim>> quadraturePoints;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    quad_.getQuadratureData(cell, grad_phi, quadraturePoints, phi, weights);

    MatrixXd M_local = MatrixXd::Zero(matSize, matSize);
    MatrixXd K_local = MatrixXd::Zero(matSize, matSize);

    for (int q=0; q<(int)quadraturePoints.size(); ++q) {
        const Point<dim>& p = quadraturePoints[q];
        const double weight = weights[q];
        const std::vector<double>& phi_q = phi[q];

        const double w_mu = weight * diffusion_.value(p);
        const Point<dim> w_b = transport_.value(p) * weight;
        const double w_r = weight * reaction_.value(p);

        for (unsigned i=0; i<matSize; ++i) {
            const double phi_i = phi_q[i];
            const Point<dim>& grad_phi_i = grad_phi[i];
            
            for (unsigned j=0; j<matSize; ++j) {
                const double phi_j = phi_q[j];
                const Point<dim>& grad_phi_j = grad_phi[j];
                
                // Pre-compute common terms
                const double phi_i_phi_j = phi_i * phi_j;
                
                // M matrix (mass)
                M_local(i,j) += weight * phi_i_phi_j;

                // K matrix components
                const double diff_c = w_mu * (grad_phi_i * grad_phi_j);
                const double adv_c  = (w_b * grad_phi_j) * phi_i;
                K_local(i,j) += (diff_c + adv_c) + w_r * phi_i_phi_j;
            }
        }
    }

    for (unsigned i=0; i<matSize; ++i) {
        const int I = cell.getNodeIndex(i);
        for (unsigned j=0; j<matSize; ++j) {
            const int J = cell.getNodeIndex(j);
            if (std::abs(M_local(i,j)) > 1e-14) tripletM.emplace_back(I,J,M_local(i,j));
            if (std::abs(K_local(i,j)) > 1e-14) tripletK.emplace_back(I,J,K_local(i,j));
        }
    }
}


// Build time-dependent forcing vector f(x,t) at given time t
template<unsigned int dim>
void FemTD<dim>::build_load(double t) {
    f_new.setZero();
    
    const unsigned matSize = dim + 1;
    const int numElements = mesh_.getNumElements();

#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        VectorXd& F_local = F_threads[tid];
        F_local.setZero();

        #pragma omp for schedule(dynamic, 32)
        for (int i = 0; i < numElements; ++i) {
            const Cell<dim>& cell = mesh_.getCell(i);

            std::vector<Point<dim>> grad_phi;
            std::vector<Point<dim>> quadraturePoints;
            std::vector<std::vector<double>> phi;
            std::vector<double> weights;
            quad_.getQuadratureData(cell, grad_phi, quadraturePoints, phi, weights);

            for (int q=0; q < quadraturePoints.size(); ++q) {
                const double forcingVal = forcing_td_(quadraturePoints[q], t);
                const double w_forcing = weights[q] * forcingVal;
                const std::vector<double>& phi_q = phi[q];
                
                for (unsigned j=0; j<matSize; ++j) {
                    const int jGlobal = cell.getNodeIndex(j);
                    F_local[jGlobal] += w_forcing * phi_q[j];
                }
            }
        }
    }
    
    // Sum all thread-local vectors
    for (int tid = 0; tid < nthreads; ++tid) {
        f_new += F_threads[tid];
    }
#else
    for (int i = 0; i < numElements; ++i) {
        const Cell<dim>& cell = mesh_.getCell(i);

        std::vector<Point<dim>> grad_phi;
        std::vector<Point<dim>> quadraturePoints;
        std::vector<std::vector<double>> phi;
        std::vector<double> weights;
        quad_.getQuadratureData(cell, grad_phi, quadraturePoints, phi, weights);

        for (int q=0; q < quadraturePoints.size(); ++q) {
            const double w_forcingVal = weights[q] * forcing_td_(quadraturePoints[q], t);
            for (unsigned j=0; j<matSize; ++j) {
                const int jGlobal = cell.getNodeIndex(j);
                f_new[jGlobal] += w_forcingVal * phi[q][j];
            }
        }
    }
#endif
}


// Apply initial condition u(x,0) = u0(x) to solution vector
template<unsigned int dim>
void FemTD<dim>::apply_initial_condition() {
    #pragma omp parallel for
    for (int i=0; i<mesh_.getNumNodes(); ++i) {
        const Point<dim>& X = mesh_.getNode(i);
        u_[i] = u0_.value(X);
    }
    u_old_ = u_;
}


// Time step: solve theta-method scheme
template<unsigned int dim>
double FemTD<dim>::step(double t_new, double dt, double theta) {
    const int N = mesh_.getNumNodes();

    f_old = f_new;
    build_load(t_new);

    rhs_ = ( (1.0/dt) * (M_ * u_old_) )
         - ( (1.0 - theta) * (K_ * u_old_) )
         + theta * f_new + (1.0 - theta) * f_old;

    bc_.apply(mesh_, A_, rhs_, t_new);

    // solve
    if (A_.nonZeros() < 4e3) {
        std::cout << "[TD] Solving small system (" << A_.nonZeros() << ") with SparseLU\n";
        Eigen::SparseLU<SparseMat> solver;
        solver.analyzePattern(A_);
        solver.factorize(A_);
        u_ = solver.solve(rhs_);
    } else {
        std::cout << "[TD] Solving large system (" << A_.nonZeros() << ") with BiCGSTAB\n";
        Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double>> solver;
        solver.setMaxIterations(10000);
        solver.setTolerance(1e-8);
        solver.compute(A_);
        u_ = solver.solve(rhs_);
    }

    const double diff = (u_ - u_old_).norm();
    u_old_ = u_;
    return diff;
}


// Time integration loop: solve time-dependent PDE from t=0 to t=T
template<unsigned int dim>
void FemTD<dim>::run(double T, double dt, double theta,
                     const std::string& vtu_prefix,
                     const std::string& csv_prefix)
{
    assemble_time_invariant();
    A_ = (1.0/dt) * M_ + theta * K_;

    apply_initial_condition();

    if (!vtu_prefix.empty()) outputVtu(vtu_prefix + "_t0.vtu");
    if (!csv_prefix.empty()) outputCsv(csv_prefix + "_t0.csv");

    double t = 0.0;
    unsigned step_id = 0;
    while (t < T - dt/2.0) {
        t += dt;
        ++step_id;
        const double change = step(t, dt, theta);

        if (!vtu_prefix.empty()) outputVtu(vtu_prefix + "_t" + std::to_string(step_id) + ".vtu");
        if (!csv_prefix.empty()) outputCsv(csv_prefix + "_t" + std::to_string(step_id) + ".csv");

        std::cout << "[TD] step " << step_id << " t=" << t
                  << "  ||u^{n+1}-u^{n}||=" << change << "\n";
    }
}

// Output solution to CSV format: coordinates and values
template<unsigned int dim>
void FemTD<dim>::outputCsv(const std::string& filename) const {
    std::ofstream csv(filename);
    if (!csv) { std::cerr << "Cannot open " << filename << "\n"; return; }
    csv << "# x " << (dim>=2?"y ":"") << (dim==3?"z ":"") << "u\n";
    for (unsigned i=0; i<mesh_.getNumNodes(); ++i) {
        const Point<dim>& X = mesh_.getNode(i);
        for (unsigned d=0; d<dim; ++d) csv << X[d] << " ";
        csv << u_[i] << "\n";
    }
}


// Output solution to VTU format for visualization (ParaView/VisIt)
template<unsigned int dim>
void FemTD<dim>::outputVtu(const std::string& filename) const {
    std::ofstream vtu(filename);
    if (!vtu) { std::cerr << "Cannot open " << filename << "\n"; return; }

    vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtu << "<UnstructuredGrid>\n";
    vtu << "<Piece NumberOfPoints=\"" << mesh_.getNumNodes()
        << "\" NumberOfCells=\"" << mesh_.getNumElements() << "\">\n";

    vtu << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Point<dim>& node : mesh_.getUniqueNodes()) {
        for (unsigned i=0; i<3; ++i) vtu << (i<dim ? node[i] : 0.0) << " ";
    }
    vtu << "\n</DataArray>\n</Points>\n";

    vtu << "<PointData Scalars=\"solution\">\n";
    vtu << "<DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n";
    for (unsigned i=0; i<mesh_.getNumNodes(); ++i) vtu << u_[i] << "\n";
    vtu << "</DataArray>\n</PointData>\n";

    vtu << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const Cell<dim>& cell : mesh_.getCells()) {
        for (unsigned i=0; i<cell.getN(); ++i) vtu << cell.getNodeIndex(i) << " ";
        vtu << "\n";
    }
    vtu << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    unsigned offset = 0;
    for (const Cell<dim>& cell : mesh_.getCells()) { offset += cell.getN(); vtu << offset << "\n"; }
    vtu << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    const unsigned vtkType = (dim==1)?3 : (dim==2)?5 : 10;
    for (unsigned i=0; i<mesh_.getNumElements(); ++i) vtu << vtkType << "\n";
    vtu << "</DataArray>\n</Cells>\n";

    vtu << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

// Explicit instantiations
template class FemTD<1>;
template class FemTD<2>;
template class FemTD<3>;
