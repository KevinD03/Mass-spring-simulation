#include <V_gravity_particle.h>

void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    // -g for downward gravitional force
    // F = mgh
    V = mass * q.transpose() * -g;
}