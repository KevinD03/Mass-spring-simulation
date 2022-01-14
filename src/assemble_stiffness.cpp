#include <assemble_stiffness.h>

//Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  K - the 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function. 

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 

    K.resize(q.rows(), q.rows());
    K.setZero();
    typedef Eigen::Triplet<double> T;
    std::vector<T> coordinate_triplets;
    coordinate_triplets.reserve(q.rows()*q.rows());

    Eigen::Matrix66d H_local_spring;

    Eigen::Vector3d q0;
    Eigen::Vector3d q1;

    for (int spring_i = 0; spring_i < E.rows(); spring_i++) {
        
        // find q0 and q1
        auto q0_start_coordinate_x0 = 3 * E(spring_i, 0);
        auto q0_start_coordinate_y0 = 3 * E(spring_i, 0) + 1;
        auto q0_start_coordinate_z0 = 3 * E(spring_i, 0) + 2;
        auto q1_start_coordinate_x1 = 3 * E(spring_i, 1);
        auto q1_start_coordinate_y1 = 3 * E(spring_i, 1) + 1;
        auto q1_start_coordinate_z1 = 3 * E(spring_i, 1) + 2;

        q0 << q(q0_start_coordinate_x0), q(q0_start_coordinate_y0), q(q0_start_coordinate_z0);
        q1 << q(q1_start_coordinate_x1), q(q1_start_coordinate_y1), q(q1_start_coordinate_z1);


        d2V_spring_particle_particle_dq2(H_local_spring, q0, q1, l0(spring_i), k);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coordinate_triplets.push_back(T(3 * E(spring_i, 0) + i, 3 * E(spring_i, 0) + j, H_local_spring(i, j)));
                coordinate_triplets.push_back(T(3 * E(spring_i, 0) + i, 3 * E(spring_i, 1) + j, H_local_spring(i, j + 3.0)));
                coordinate_triplets.push_back(T(3 * E(spring_i, 1) + i, 3 * E(spring_i, 0) + j, H_local_spring(i + 3.0, j)));
                coordinate_triplets.push_back(T(3 * E(spring_i, 1) + i, 3 * E(spring_i, 1) + j, H_local_spring(i + 3.0, j + 3.0)));
            }
        }
    }
    K.setFromTriplets(coordinate_triplets.begin(), coordinate_triplets.end());
};