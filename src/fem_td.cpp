#include "fem_td.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

// ===== ctor / set =====
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
    const int N = mesh_.getNumNodes();
    M_.resize(N,N);
    K_.resize(N,N);
    A_.resize(N,N);
    rhs_.resize(N); rhs_.setZero();
    u_.resize(N);   u_.setZero();
    u_old_.resize(N); u_old_.setZero();
}

template<unsigned int dim>
void FemTD<dim>::set_forcing(ForcingTD f_td) { forcing_td_ = std::move(f_td); }

template<unsigned int dim>
void FemTD<dim>::set_initial_condition(Function<dim,1> u0) { u0_ = std::move(u0); }

// ===== assemble M & K (time-invariant) =====
template<unsigned int dim>
void FemTD<dim>::assemble_time_invariant() {
    const int ne = mesh_.getNumElements();
    const unsigned matSize = dim + 1;
    const unsigned expPerElem = matSize * matSize;

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    std::vector<std::vector<Triplet>> tripM_thr(nthreads), tripK_thr(nthreads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<Triplet>& tM = tripM_thr[tid];
        std::vector<Triplet>& tK = tripK_thr[tid];
        tM.reserve(expPerElem * (ne/nthreads + 1));
        tK.reserve(expPerElem * (ne/nthreads + 1));

        #pragma omp for schedule(static)
        for (int e=0; e<ne; ++e) assemble_M_and_K_element(e, tM, tK);
    }
    std::vector<Triplet> tripM, tripK;
    for (std::vector<Triplet>& v: tripM_thr) tripM.insert(tripM.end(), v.begin(), v.end());
    for (std::vector<Triplet>& v: tripK_thr) tripK.insert(tripK.end(), v.begin(), v.end());
#else
    std::vector<Triplet> tripM, tripK;
    tripM.reserve(expPerElem * ne);
    tripK.reserve(expPerElem * ne);
    for (int e=0; e<ne; ++e) assemble_M_and_K_element(e, tripM, tripK);
#endif

    M_.setFromTriplets(tripM.begin(), tripM.end());
    K_.setFromTriplets(tripK.begin(), tripK.end());
}

template<unsigned int dim>
void FemTD<dim>::assemble_M_and_K_element(int elem,
                                          std::vector<Triplet>& tripM,
                                          std::vector<Triplet>& tripK)
{
    const Cell<dim>& cell = mesh_.getCell(elem);
    const unsigned matSize = dim + 1;

    std::vector<Point<dim>> grad_phi;
    std::vector<Point<dim>> xq;
    std::vector<std::vector<double>> phi;
    std::vector<double> w;
    quad_.getQuadratureData(cell, grad_phi, xq, phi, w);

    MatrixXd Mloc = MatrixXd::Zero(matSize, matSize);
    MatrixXd Kloc = MatrixXd::Zero(matSize, matSize);

    for (int q=0; q<(int)xq.size(); ++q) {
        const Point<dim>& p = xq[q];
        const double ww = w[q];

        const double mu = diffusion_.value(p);
        const Point<dim> b = transport_.value(p);
        const double r = reaction_.value(p);

        for (unsigned i=0; i<matSize; ++i) {
            for (unsigned j=0; j<matSize; ++j) {
                // M
                Mloc(i,j) += ww * phi[q][i] * phi[q][j];

                // K = μ ∇φ_i·∇φ_j + (b·∇φ_i) φ_j + r φ_i φ_j
                const double diff_c = mu * (grad_phi[i] * grad_phi[j]);
                const double adv_c  = (b * grad_phi[i]) * phi[q][j];
                const double reac_c = r * phi[q][i] * phi[q][j];
                Kloc(i,j) += ww * (diff_c + adv_c + reac_c);
            }
        }
    }

    for (unsigned i=0; i<matSize; ++i) {
        const int I = cell.getNodeIndex(i);
        for (unsigned j=0; j<matSize; ++j) {
            const int J = cell.getNodeIndex(j);
            if (std::abs(Mloc(i,j)) > 1e-14) tripM.emplace_back(I,J,Mloc(i,j));
            if (std::abs(Kloc(i,j)) > 1e-14) tripK.emplace_back(I,J,Kloc(i,j));
        }
    }
}

// ===== RHS load =====
template<unsigned int dim>
void FemTD<dim>::build_load(VectorXd& F, double t) const {
    F.setZero();
    
    const unsigned matSize = dim + 1;

    #pragma omp parallel for reduction(+:F)
    for (int i = 0; i < mesh_.getNumElements(); ++i) {
        const Cell<dim>& cell = mesh_.getCell(i);

        std::vector<Point<dim>> grad_phi;
        std::vector<Point<dim>> xq;
        std::vector<std::vector<double>> phi;
        std::vector<double> w;
        // Rimuovi const usando un cast
        quad_.getQuadratureData(cell, grad_phi, xq, phi, w);

        for (int q=0; q<(int)xq.size(); ++q) {
            const double fval = forcing_td_(xq[q], t);
            const double ww = w[q];
            for (unsigned i=0; i<matSize; ++i) {
                const int I = cell.getNodeIndex(i);
                F[I] += ww * fval * phi[q][i];
            }
        }
    }
}\

// ===== IC =====
template<unsigned int dim>
void FemTD<dim>::apply_initial_condition() {
    #pragma omp parallel for
    for (int i=0; i<mesh_.getNumNodes(); ++i) {
        const Point<dim>& X = mesh_.getNode(i);
        u_[i] = u0_.value(X);
    }
    u_old_ = u_;
}

// ===== singolo passo =====
template<unsigned int dim>
double FemTD<dim>::step(double t_new, double dt, double theta) {
    const int N = mesh_.getNumNodes();

    VectorXd f_new(N), f_old(N);
    build_load(f_new, t_new);
    build_load(f_old, t_new - dt);

    rhs_ = ( (1.0/dt) * (M_ * u_old_) )
         - ( (1.0 - theta) * (K_ * u_old_) )
         + theta * f_new + (1.0 - theta) * f_old;

    // BC al passo
    bc_.apply(mesh_, A_, rhs_, t_new);

    // solve
    if (A_.nonZeros() < 10000) {
        std::cout << "[TD] Solving small system with SparseLU\n";
        Eigen::SparseLU<SparseMat> solver;
        solver.analyzePattern(A_);
        solver.factorize(A_);
        u_ = solver.solve(rhs_);
    } else {
        Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double>> solver;
        std::cout << "[TD] Solving large system with BiCGSTAB\n";
        solver.setMaxIterations(10000);
        solver.setTolerance(1e-8);
        solver.compute(A_);
        u_ = solver.solve(rhs_);
    }

    const double diff = (u_ - u_old_).norm();
    u_old_ = u_;
    return diff;
}

// ===== run =====
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
    while (t < T - 1e-14) {
        t += dt;
        ++step_id;
        const double change = step(t, dt, theta);

        if (!vtu_prefix.empty()) outputVtu(vtu_prefix + "_t" + std::to_string(step_id) + ".vtu");
        if (!csv_prefix.empty()) outputCsv(csv_prefix + "_t" + std::to_string(step_id) + ".csv");

        std::cout << "[TD] step " << step_id << " t=" << t
                  << "  ||u^{n+1}-u^{n}||=" << change << "\n";
    }
}

// ===== output =====
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

// ===== ISTANZIAZIONI ESPLICITE =====
template class FemTD<1>;
template class FemTD<2>;
template class FemTD<3>;
