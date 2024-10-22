// This class numerically solves for the position of a flexible
// magnetic tentacle under homogenous fields

#include "SolverTwoD.h"

// Constructor without external forces
SolverTwoD::SolverTwoD(int n, double L, double Bxin, double Bzin, const Eigen::MatrixXd& magMoments, int orientation)
    : n_(n), L_(L), Bxin_(Bxin), Bzin_(Bzin), orientation_(orientation), magMoments_(magMoments),
    use_external_forces_(false) {
    initialise();
    solve();
}

// Constructor with external forces
SolverTwoD::SolverTwoD(int n, double L, double Bxin, double Bzin, const Eigen::MatrixXd& magMoments, int orientation,
    const Eigen::VectorXd& Fx, const Eigen::VectorXd& Fy)
    : n_(n), L_(L), Bxin_(Bxin), Bzin_(Bzin), orientation_(orientation), magMoments_(magMoments),
    use_external_forces_(true), Fx_(Fx), Fz_(Fy) {
    if (Fx_.size() != n_ || Fz_.size() != n_) {
        throw std::invalid_argument("Size of external force vectors must match the number of nodes.");
    }
    initialise();
    solve();
}

// Destructor
SolverTwoD::~SolverTwoD() {
}

// Initialise variables and precompute constants
void SolverTwoD::initialise() {
    // Resize vectors and matrices
    xpos_.resize(n_);
    zpos_.resize(n_);
    theta_.resize(n_);
    theta_p_.resize(n_);
    theta_p_old_.resize(n_);
    r_.resize(n_, 2);
    m_def_.resize(n_, 2);
    m_hat_.resize(n_, 2);

    // Set delta_s_ as the link length between nodes
    delta_s_ = L_ / ((static_cast<double>(n_)) - 1);

    // Initialize theta_p_old_ to a small value
    theta_p_old_.setConstant(1e-6);

    // Geometric Parameters
    double radius = 2e-3;                      // Radius (m)
    double A = M_PI * radius * radius;         // Cross-sectional area (m^2)
    double V = A;                              // Volume per unit length (m^2)
    double I = 0.25 * M_PI * std::pow(radius, 4); // Second moment of area (m^4)

    // Material Properties
    double E = 1e5;                            // Young's modulus (N/m^2)
    double rho = 1100;                         // Density (kg/m^3)
    double g = -9.81;                          // Acceleration due to gravity (m/s^2)

    // Magnetic Properties
    double mu0 = 4 * M_PI * 1e-7;              // Magnetic permeability (H/m)
    double remanence = 0.05;                   // Saturation remanent magnetic field (T)
    double M = remanence * V / mu0;            // Magnitude of magnetic moment (A·m^2)

    // Scale magMoments_ by M (element-wise multiplication)
    magMoments_ = magMoments_ * M;

    // Precompute constants
    const3_ = E * I;                         // Bending stiffness
    const1_ = (rho * g * A) / const3_;
    const2_ = 1.0 / const3_;

    // Set initial orientation
    if (orientation_ == 0) {
        theta_.setZero();
    }
    else if (orientation_ == 1) {
        theta_.setConstant(-M_PI / 2);
    }
    else {
        throw std::invalid_argument("Invalid orientation. Use 0 for horizontal, 1 for vertically downwards.");
    }

    // Initialize first position
    r_.row(0) << 0.0, 0.0;
}

// Iteratively solves for the tentacle positions
void SolverTwoD::solve() {
    int max_loops = 100;        // Maximum number of iterations
    double mu = 0.5;            // Numerical damping parameter
    double err_tol = 1e-6;      // Error tolerance for convergence

    for (int j = 0; j < max_loops; ++j) {
        // Compute cos(theta) and sin(theta)
        Eigen::VectorXd cos_theta = theta_.array().cos();
        Eigen::VectorXd sin_theta = theta_.array().sin();

        // Compute deformed magnetic moments (m_def_) and their derivatives (m_hat_)
        m_def_.col(0) = cos_theta.array() * magMoments_.col(0).array() - sin_theta.array() * magMoments_.col(1).array();
        m_def_.col(1) = sin_theta.array() * magMoments_.col(0).array() + cos_theta.array() * magMoments_.col(1).array();

        m_hat_.col(0) = -sin_theta.array() * magMoments_.col(0).array() - cos_theta.array() * magMoments_.col(1).array();
        m_hat_.col(1) = cos_theta.array() * magMoments_.col(0).array() - sin_theta.array() * magMoments_.col(1).array();

        // Compute positions of nodes
        Eigen::VectorXd dx = delta_s_ * cos_theta;
        Eigen::VectorXd dz = delta_s_ * sin_theta;

        // Compute cumulative sums for positions
        Eigen::VectorXd cum_dx = Eigen::VectorXd::Zero(n_);
        Eigen::VectorXd cum_dz = Eigen::VectorXd::Zero(n_);

        cum_dx.segment(1, n_ - 1) = SolverTwoD::cumulativeSum(dx.head(n_ - 1));
        cum_dz.segment(1, n_ - 1) = SolverTwoD::cumulativeSum(dz.head(n_ - 1));

        // Set positions
        r_.col(0) = cum_dx;
        r_.col(1) = cum_dz;

        // Compute dot product of m_hat_ and magnetic field components
        Eigen::VectorXd dot_product = m_hat_.col(0) * Bxin_ + m_hat_.col(1) * Bzin_;

        // Initialize expr with gravitational term
        Eigen::VectorXd expr = const1_ * cos_theta;

        // Add magnetic term
        expr.array() += theta_p_old_.array().square() * (const2_ * dot_product.array());

        // Add external forces if applicable
        if (use_external_forces_) {

            // External forces per unit length divided by EI
            Eigen::VectorXd Fx_term = Fx_ / (const3_);
            Eigen::VectorXd Fy_term = Fz_ / (const3_);

            // Project external forces onto the local coordinate system
            Eigen::VectorXd external_force_term = cos_theta.array() * Fy_term.array() - sin_theta.array() * Fx_term.array();

            // Add to expr
            expr += external_force_term;
        }

        // Apply cube root function
        theta_p_ = (1.0 - mu) * theta_p_old_ + mu * expr.unaryExpr([](double x) {
            return std::copysign(std::pow(std::abs(x), 1.0 / 3.0), x);
            });
        // Eigen::Matrix.unaryExpr applies an operation to each element in the array.
        // The operation is a lambda
        // The lambda ensures that the cube root doesn't produce complex numbers

        // Integrate theta_p_ to get theta_
        theta_.segment(1, n_ - 1) = theta_.segment(0, n_ - 1) + ((theta_p_.segment(1, n_ - 1) + theta_p_.segment(0, n_ - 1)) / 2.0) * delta_s_;

        // Check for convergence
        double max_diff = (theta_p_ - theta_p_old_).array().abs().maxCoeff();
        if (max_diff < err_tol) {
            break;
        }

        // Update theta_p_old_ for next iteration
        theta_p_old_ = theta_p_;
    }

    // Store final positions
    xpos_ = r_.col(0);
    zpos_ = r_.col(1);
}

// cumulative sum helper function
Eigen::VectorXd SolverTwoD::cumulativeSum(const Eigen::VectorXd& v) {
    Eigen::VectorXd result(v.size());
    if (v.size() == 0) return result;
    result[0] = v[0];
    for (Eigen::Index i = 1; i < v.size(); ++i) {
        result[i] = result[i - 1] + v[i];
    }
    return result;
}

// Accessor for X positions
const Eigen::VectorXd& SolverTwoD::getXPositions() const {
    return xpos_;
}

// Accessor for Y positions
const Eigen::VectorXd& SolverTwoD::getZPositions() const {
    return zpos_;
}
