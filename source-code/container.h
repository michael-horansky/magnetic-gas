#ifndef C_CONTAINER
#define C_CONTAINER

#include <vector>
#include <chrono>
#include <random>
#include <math.h>
#include "Journal.cpp"

double permitivity = 0.000000000008854187817;
double permeability = 0.0000012566370614;
double c = sqrt(1.0/(permitivity*permeability));
double boltzmann_constant = 0.00000000000000000000001380650424;
double avogadro_constant = 602214150000000000000000.0;
double elementary_charge = 0.00000000000000000016021766208;

std::uniform_real_distribution<double> rand_phi(0.0, 2.0*M_PIl);
std::uniform_real_distribution<double> rand_theta(0.0, M_PIl);
std::uniform_real_distribution<double> rand_def(0.0, 1.0);
std::default_random_engine random_engine;

std::vector<double> vector_direction(double v_x, double v_y, double v_z) {
    std::vector<double> result(2,0.0);
    if (v_y == 0 && v_x == 0) {
      result[1] = 0;
    } else {
      result[1] = M_PIl / 2 - atan(v_z / sqrt(v_y*v_y + v_x*v_x));
    }
    if (v_y == 0) {
      result[0] = 0;
    } else {
      result[0] = atan(v_x / v_y);
    }
    return result;
}
std::vector<double> vector_components(double r, double phi, double theta) {
  std::vector<double> result(3,0.0);
  result[0] = r * sin(theta) * cos(phi);
  result[1] = r * sin(theta) * sin(phi);
  result[2] = r * cos(theta);
  return result;
}
std::vector<double> cross_product(double a_x, double a_y, double a_z, double b_x, double b_y, double b_z) {
  std::vector<double> result(3,0.0);
  result[0] = a_y * b_z - a_z * b_y;
  result[1] = a_z * b_x - a_x * b_z;
  result[2] = a_x * b_y - a_y * b_x;
  return result;
}
std::vector<double> rotate_vector(double x, double y, double z, double alpha_x, double alpha_y, double alpha_z) {
  std::vector<double> result(3,0.0);
  double sinx = sin(alpha_x);
  double cosx = cos(alpha_x);
  double siny = sin(alpha_y);
  double cosy = cos(alpha_y);
  double sinz = sin(alpha_z);
  double cosz = cos(alpha_z);
  result[0] = x * cosz * cosy + y * (cosz * siny * sinx - sinz * cosx) + z * (sinz * sinx + cosz * siny * cosx);
  result[1] = x * sinz * cosy + y * (cosz * cosx + sinz * siny * sinx) + z * (sinz * siny * cosx - cosz * sinx);
  result[2] = -x * siny       + y * cosy * sinx                        + z * cosy * cosx;
  return result;
}

class Particle;

struct ExternalMagneticDipole {

  double p_x;
  double p_y;
  double p_z;
  double mag_m_x;
  double mag_m_y;
  double mag_m_z;

};

class Container {

  public:
    Container();
    Container(double my_dt);
    Container(double my_dt, std::string journal_name);
    ~Container();

    void cuboid_shape(double my_d_x, double my_d_y, double my_d_z);
    void cylindrical_shape(double my_radius, double my_length);

    void add_external_dipole(double my_p_x, double my_p_y, double my_p_z, double my_mag_m_x, double my_mag_m_y, double my_mag_m_z);
    void external_magnetic_field(double m_B_x, double m_B_y, double m_B_z);

    void add_determined_particle(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_r);
    void add_particle(double m_v, double m_m, double m_r);
    void add_determined_electric_particle(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_q, double m_r);
    void add_electric_particle(double m_v, double m_m, double m_q, double m_r);
    void add_determined_magnetic_dipole(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_q, double m_mag_x, double m_mag_y, double m_mag_z, double m_r);
    void add_magnetic_dipole(double m_v, double m_m, double m_q, double m_mag, double m_r);

    void tick();
    void Simulate(double max_t, double check_time, double data_points);
    void Simulate(double max_t, double check_time, double data_points, double omega, long int magnet_distance, long int permanent_magnet_force_coef);

    std::uniform_real_distribution<double> rand_x;
    std::uniform_real_distribution<double> rand_y;
    std::uniform_real_distribution<double> rand_z;

    double dt;
    double d_x;
    double d_y;
    double d_z;
    double radius;
    double length;

    int shape; //0 = cuboid, 1 = cylindrical
    double surface;

    double t;
    double kinetic_energy;
    double temperature;
    double surface_impulse;
    double cylindrical_impulse;
    double xy1_impulse;
    double ext_B_x;
    double ext_B_y;
    double ext_B_z;
    std::vector<ExternalMagneticDipole*> ext_dipoles;

    std::vector<Particle*> particles;

    Journal *my_journal;
    int is_journal;
    std::ofstream output_file;
    void print_particles();

    std::ofstream velocity_distribution_output_file;
    void SaveVelocityDistribution();

    std::ofstream energy_output_file;
    std::ofstream pressure_output_file;
    std::ofstream temperature_output_file;
    std::ofstream cloud_chamber_output_file;
    std::ofstream enthropy_output_file;
    int save_energy;
    int save_thermodynamics;
    int save_cloud_chamber;
    void SaveEnergy();
    void SaveThermodynamics();
    void CloudChamber();
    double delta_alpha;
};

#endif // C_CONTAINER
