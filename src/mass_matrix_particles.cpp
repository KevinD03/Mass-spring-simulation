#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    M.resize(q.size(), q.size());
    // make M a mass "identity" matrix
    M.setIdentity();
    for (int i = 0; i < M.rows(); i++) {
        M.insert(i, i) = mass;
    }
    
}
