#ifndef C_MAGNETIC_DIPOLE
#define C_MAGNETIC_DIPOLE

#include "electric_particle.cpp"

class MagneticDipole : public ElectricParticle {

  public:
    MagneticDipole();
    MagneticDipole(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_charge, double start_magnetic_momentum_x, double start_magnetic_momentum_y, double start_magnetic_momentum_z, double start_radius);
    ~MagneticDipole();

    void randomize_particle(double max_vel, double max_magnetic_moment);

    virtual std::vector<double> force_imposed_by(Particle *other);
    virtual std::vector<double> force_imposed_by_container();

    void rotate_self();

    double mag_m_x;
    double mag_m_y;
    double mag_m_z;

};

#endif // C_MAGNETIC_DIPOLE
