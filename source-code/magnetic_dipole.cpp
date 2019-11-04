#include "magnetic_dipole.h"


MagneticDipole :: MagneticDipole() {

}

MagneticDipole :: MagneticDipole(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_charge, double start_magnetic_momentum_x, double start_magnetic_momentum_y, double start_magnetic_momentum_z, double start_radius)
                : ElectricParticle(my_parent_container, start_pos_x, start_pos_y, start_pos_z, start_vel_x, start_vel_y, start_vel_z, start_mass, start_charge, start_radius) {

  this->mag_m_x = start_magnetic_momentum_x;
  this->mag_m_y = start_magnetic_momentum_y;
  this->mag_m_z = start_magnetic_momentum_z;

}

MagneticDipole :: ~MagneticDipole() {

}

void MagneticDipole :: randomize_particle(double max_vel, double max_magnetic_momentum) {

  double phi, theta;
  this->Particle::randomize_particle(max_vel);
  phi = rand_phi(random_engine);
  theta = rand_theta(random_engine);
  this->mag_m_x = max_magnetic_momentum * sin(theta) * cos(phi);
  this->mag_m_y = max_magnetic_momentum * sin(theta) * sin(phi);
  this->mag_m_z = max_magnetic_momentum * cos(theta);

}

std::vector<double> MagneticDipole :: force_imposed_by(Particle *other) {

  std::vector<double> result(6,0.0), fields(6,0.0), magnetic_torque(3, 0.0), magnetic_field(3,0.0), lorenz_force(3,0.0);
  MagneticDipole *r = (MagneticDipole*) other;
  double safe_distance = 0.01;
  double B_x, B_y, B_z, dx, dy, dz, d, d_5, m_1_r, m_1_x, m_1_y, m_1_z, m_2_r, m_1_m_2, m_2_x, m_2_y, m_2_z;
  double magnetic_coeffitient;

  dx = this->p_x - other->p_x;
  dy = this->p_y - other->p_y;
  dz = this->p_z - other->p_z;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  d_5 = d * d * d * d * d;
  if (d < safe_distance) d = safe_distance;

  fields = this->electromagnetic_field_induced_by((ElectricParticle *) other);

  m_1_x = r->mag_m_x;
  m_1_y = r->mag_m_y;
  m_1_z = r->mag_m_z;
  m_2_x = this->mag_m_x;
  m_2_y = this->mag_m_y;
  m_2_z = this->mag_m_z;
  m_1_r = m_1_x * dx + m_1_y * dy + m_1_z * dz;
  m_2_r = m_2_x * dx + m_2_y * dy + m_2_z * dz;
  m_1_m_2 = m_1_x * m_2_x + m_1_y * m_2_y + m_1_z * m_2_z;

  magnetic_field[0] = fields[3] + permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dx / (d_5) - m_1_x / (d * d * d));
  magnetic_field[1] = fields[4] + permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dy / (d_5) - m_1_y / (d * d * d));
  magnetic_field[2] = fields[5] + permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dz / (d_5) - m_1_z / (d * d * d));

  lorenz_force = cross_product(this->v_x, this->v_y, this->v_z, magnetic_field[0], magnetic_field[1], magnetic_field[2]);
  result[0] += this->q * (fields[0] + lorenz_force[0]);
  result[1] += this->q * (fields[1] + lorenz_force[1]);
  result[2] += this->q * (fields[2] + lorenz_force[2]);

  magnetic_coeffitient = 3.0 * permeability / (4.0 * M_PIl * (d_5));
  result[0] += magnetic_coeffitient * (m_1_r * m_2_x + m_2_r * m_1_x + m_1_m_2 * dx - (5.0 * m_1_r * m_2_r) / (d * d) * dx);
  result[1] += magnetic_coeffitient * (m_1_r * m_2_y + m_2_r * m_1_y + m_1_m_2 * dy - (5.0 * m_1_r * m_2_r) / (d * d) * dy);
  result[2] += magnetic_coeffitient * (m_1_r * m_2_z + m_2_r * m_1_z + m_1_m_2 * dz - (5.0 * m_1_r * m_2_r) / (d * d) * dz);

  magnetic_torque = cross_product(m_2_x, m_2_y, m_2_z, magnetic_field[0], magnetic_field[1], magnetic_field[2]);
  result[3] = magnetic_torque[0];
  result[4] = magnetic_torque[1];
  result[5] = magnetic_torque[2];

  return result;

}


std::vector<double> MagneticDipole :: force_imposed_by_container() {

  std::vector<double> result(6,0.0), lorenz_force(3,0.0), magnetic_torque(3, 0.0), v_x_B(3,0.0), v_perp_B(3,0.0), radius(3,0.0), new_v(3,0.0), new_p(3,0.0);
  double B_x, B_y, B_z, B, B_sq, dx, dy, dz, d, d_5, m_1_r, m_1_x, m_1_y, m_1_z, m_2_r, m_1_m_2, m_2_x, m_2_y, m_2_z, r, F_L, v;
  double magnetic_coeffitient, magnetic_force_x, magnetic_force_y, magnetic_force_z;
  int i;

  m_2_x = this->mag_m_x;
  m_2_y = this->mag_m_y;
  m_2_z = this->mag_m_z;

  B_x = this->parent_container->ext_B_x;
  B_y = this->parent_container->ext_B_y;
  B_z = this->parent_container->ext_B_z;
  magnetic_force_x = 0;
  magnetic_force_y = 0;
  magnetic_force_z = 0;

  for(i = 0; i < this->parent_container->ext_dipoles.size(); i++) {

    dx = this->p_x - this->parent_container->ext_dipoles[i]->p_x;
    dy = this->p_y - this->parent_container->ext_dipoles[i]->p_y;
    dz = this->p_z - this->parent_container->ext_dipoles[i]->p_z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    d_5 = d * d * d * d * d;

    m_1_x = this->parent_container->ext_dipoles[i]->mag_m_x;
    m_1_y = this->parent_container->ext_dipoles[i]->mag_m_y;
    m_1_z = this->parent_container->ext_dipoles[i]->mag_m_z;
    m_1_r = m_1_x * dx + m_1_y * dy + m_1_z * dz;
    m_2_r = m_2_x * dx + m_2_y * dy + m_2_z * dz;
    m_1_m_2 = m_1_x * m_2_x + m_1_y * m_2_y + m_1_z * m_2_z;

    B_x += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dx / (d_5) - m_1_x / (d * d * d));
    B_y += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dy / (d_5) - m_1_y / (d * d * d));
    B_z += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dz / (d_5) - m_1_z / (d * d * d));

    magnetic_coeffitient = 3.0 * permeability / (4.0 * M_PIl * (d_5));
    magnetic_force_x += magnetic_coeffitient * (m_1_r * m_2_x + m_2_r * m_1_x + m_1_m_2 * dx - (5.0 * m_1_r * m_2_r) / (d * d) * dx);
    magnetic_force_y += magnetic_coeffitient * (m_1_r * m_2_y + m_2_r * m_1_y + m_1_m_2 * dy - (5.0 * m_1_r * m_2_r) / (d * d) * dy);
    magnetic_force_z += magnetic_coeffitient * (m_1_r * m_2_z + m_2_r * m_1_z + m_1_m_2 * dz - (5.0 * m_1_r * m_2_r) / (d * d) * dz);

  }
  magnetic_torque = cross_product(m_2_x, m_2_y, m_2_z, B_x, B_y, B_z);

  B_sq = B_x * B_x + B_y * B_y + B_z * B_z;
  B = sqrt(B_sq);

  v_perp_B[0] = this->v_x - B_x * this->v_x * B_x / B_sq;
  v_perp_B[1] = this->v_y - B_y * this->v_y * B_y / B_sq;
  v_perp_B[2] = this->v_z - B_z * this->v_z * B_z / B_sq;
  v = sqrt(v_perp_B[0] * v_perp_B[0] + v_perp_B[1] * v_perp_B[1] + v_perp_B[2] * v_perp_B[2]);
  v_x_B = cross_product(this->v_x, this->v_y, this->v_z, B_x, B_y, B_z);
  F_L = sqrt(v_x_B[0] * v_x_B[0] + v_x_B[1] * v_x_B[1] + v_x_B[2] * v_x_B[2]);
  if(F_L != 0.0) {
    r = this->m * (v * v) / (this->q * F_L);
    radius[0] = r * v_x_B[0] / F_L;
    radius[1] = r * v_x_B[1] / F_L;
    radius[2] = r * v_x_B[2] / F_L;
    new_v = rotate_vector(this->v_x, this->v_y, this->v_z, B_x / B * v / r * this->parent_container->dt, B_y / B * v / r * this->parent_container->dt, B_z / B * v / r * this->parent_container->dt);
    new_p = rotate_vector(radius[0], radius[1], radius[2], B_x / B * v / r * this->parent_container->dt, B_y / B * v / r * this->parent_container->dt, B_z / B * v / r * this->parent_container->dt);
    this->p_x = this->p_x - radius[0] + new_p[0];
    this->p_y = this->p_y - radius[1] + new_p[1];
    this->p_z = this->p_z - radius[2] + new_p[2];
    this->v_x = new_v[0];
    this->v_y = new_v[1];
    this->v_z = new_v[2];
  }
  result[0] = magnetic_force_x;
  result[1] = magnetic_force_y;
  result[2] = magnetic_force_z;
  result[3] = magnetic_torque[0];
  result[4] = magnetic_torque[1];
  result[5] = magnetic_torque[2];
  return result;

}

void MagneticDipole :: rotate_self() {

  double dt = this->parent_container->dt;
  std::vector<double> new_mag_m = rotate_vector(this->mag_m_x,this->mag_m_y,this->mag_m_z,this->omega_x * dt + 1.0 / 2.0 * this->eps_x * dt * dt, this->omega_y * dt + 1.0 / 2.0 * this->eps_y * dt * dt, this->omega_z * dt + 1.0 / 2.0 * this->eps_z * dt * dt);
  this->mag_m_x = new_mag_m[0];
  this->mag_m_y = new_mag_m[1];
  this->mag_m_z = new_mag_m[2];
  this->omega_x += this->eps_x * dt;
  this->omega_y += this->eps_y * dt;
  this->omega_z += this->eps_z * dt;

}
