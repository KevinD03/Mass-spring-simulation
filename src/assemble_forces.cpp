#include <assemble_forces.h>
#include <iostream>

//Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains two indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  m - the mass of each particle in the mass-spring system
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 

    f.resize(q.rows());
    f.setZero();

    Eigen::Vector3d q0;
    Eigen::Vector3d q1;

    Eigen::Vector3d g(0, 0, 0);
    Eigen::Vector3d gravity;
    dV_gravity_particle_dq(gravity, mass, g);

    Eigen::Vector6d force_per_spring;

    for (int spring_i = 0; spring_i < E.rows(); spring_i++) {


        // g(0, 0, 9.8)?

        // find q0 and q1
        // *3 because every coordinate comes in x,y,z , so we jump every 3 value
        auto q0_start_coordinate_x0 = 3 * E(spring_i, 0);
        auto q0_start_coordinate_y0 = 3 * E(spring_i, 0) + 1;
        auto q0_start_coordinate_z0 = 3 * E(spring_i, 0) + 2;
        auto q1_start_coordinate_x1 = 3 * E(spring_i, 1);
        auto q1_start_coordinate_y1 = 3 * E(spring_i, 1) + 1;
        auto q1_start_coordinate_z1 = 3 * E(spring_i, 1) + 2;

        q0 << q(q0_start_coordinate_x0), q(q0_start_coordinate_y0), q(q0_start_coordinate_z0);
        q1 << q(q1_start_coordinate_x1), q(q1_start_coordinate_y1), q(q1_start_coordinate_z1);

        // Calculate f force matrix

        dV_spring_particle_particle_dq(force_per_spring, q0, q1, l0(spring_i), k);

        //Distribute the force back to globel force vector
        f(q0_start_coordinate_x0) -= (force_per_spring(0) + gravity(0));
        f(q0_start_coordinate_y0) -= (force_per_spring(1) + gravity(1));
        f(q0_start_coordinate_z0) -= (force_per_spring(2) + gravity(2));
        f(q1_start_coordinate_x1) -= (force_per_spring(3) + gravity(0));
        f(q1_start_coordinate_y1) -= (force_per_spring(4) + gravity(1));
        f(q1_start_coordinate_z1) -= (force_per_spring(5) + gravity(2));
    }
};