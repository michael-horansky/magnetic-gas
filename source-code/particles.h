#ifndef C_PARTICLE
#define C_PARTICLE

#include "container.h"

class Container;

class Particle {

  public:
    Particle();
    Particle(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_inertia);
    ~Particle();

    void randomize_particle(double max_vel);

    virtual std::vector<double> force_imposed_by(Particle *other) {return std::vector<double>(6,0.0);}
    virtual std::vector<double> force_imposed_by_container()      {return std::vector<double>(6,0.0);}

    virtual void rotate_self() {}

    double p_x;
    double p_y;
    double p_z;
    double v_x;
    double v_y;
    double v_z;
    double a_x;
    double a_y;
    double a_z;
    double next_a_x;
    double next_a_y;
    double next_a_z;
    double m;
    double r;
    double inertia;

    double omega_x;
    double omega_y;
    double omega_z;
    double eps_x;
    double eps_y;
    double eps_z;

    std::vector<double> cloud_p_x;
    std::vector<double> cloud_p_y;

    Container *parent_container;

};

#endif // C_PARTICLE
