#ifndef C_ELECTRIC_PARTICLE
#define C_ELECTRIC_PARTICLE

#include "particles.cpp"

class ElectricParticle : public Particle {

  public:
    ElectricParticle();
    ElectricParticle(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_charge, double start_radius);
    ~ElectricParticle();

    std::vector<double> electromagnetic_field_induced_by(ElectricParticle* other);

    virtual std::vector<double> force_imposed_by(Particle *other);
    virtual std::vector<double> force_imposed_by_container();

    double q;

};

#endif // C_ELECTRIC_PARTICLE
