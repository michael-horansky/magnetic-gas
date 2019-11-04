#include "electric_particle.h"

ElectricParticle :: ElectricParticle() {

}

ElectricParticle :: ElectricParticle(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_charge, double start_radius)
                  : Particle(my_parent_container, start_pos_x, start_pos_y, start_pos_z, start_vel_x, start_vel_y, start_vel_z, start_mass, start_radius)
{

  this->q = start_charge;

}

ElectricParticle :: ~ElectricParticle() {

}

std::vector<double> ElectricParticle :: electromagnetic_field_induced_by(ElectricParticle *other) {

  if(this == other) return std::vector<double>(6,0.0);

  double safe_distance = 0.01;
  double dx, dy, dz, d, magnetic_coeffitient;
  std::vector<double> result(6,0.0);
  std::vector<double> mag_field;
  dx = this->p_x - other->p_x;
  dy = this->p_y - other->p_y;
  dz = this->p_z - other->p_z;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d < safe_distance) d = safe_distance;
  result[0] = other->q / (4 * M_PIl * permitivity) * dx/(d * d * d);
  result[1] = other->q / (4 * M_PIl * permitivity) * dy/(d * d * d);
  result[2] = other->q / (4 * M_PIl * permitivity) * dz/(d * d * d);

  mag_field = cross_product(other->v_x, other->v_y, other->v_z, dx, dy, dz);
  magnetic_coeffitient = other->q * permeability / (4.0 * M_PIl) / (d * d * d);
  result[3] = mag_field[0] * magnetic_coeffitient;
  result[4] = mag_field[1] * magnetic_coeffitient;
  result[5] = mag_field[2] * magnetic_coeffitient;
  return result;

}

std::vector<double> ElectricParticle :: force_imposed_by(Particle *other) {

  if(this == other) return std::vector<double>(6,0.0);
  ElectricParticle* r = (ElectricParticle*) other;
  std::vector<double> result(6,0.0), field, mag_force(3,0.0);
  field = this->electromagnetic_field_induced_by(r);
  result[0] = field[0] * this->q;
  result[1] = field[1] * this->q;
  result[2] = field[2] * this->q;
  mag_force = cross_product(this->v_x, this->v_y, this->v_z, field[3], field[4], field[5]);
  result[3] = mag_force[0] * this->q;
  result[4] = mag_force[1] * this->q;
  result[5] = mag_force[2] * this->q;
  return result;

}

std::vector<double> ElectricParticle :: force_imposed_by_container() {

  std::vector<double> result(6,0.0), lorenz_force(3,0.0), v_x_B(3,0.0), v_perp_B(3,0.0), radius(3,0.0), new_v(3,0.0), new_p(3,0.0);
  double B_x, B_y, B_z, B, B_sq, dx, dy, dz, d, m_1_r, m_1_x, m_1_y, m_1_z, r, F_L, v;
  int i;

  B_x = this->parent_container->ext_B_x;
  B_y = this->parent_container->ext_B_y;
  B_z = this->parent_container->ext_B_z;

  for(i = 0; i < this->parent_container->ext_dipoles.size(); i++) {

    dx = this->p_x - this->parent_container->ext_dipoles[i]->p_x;
    dy = this->p_y - this->parent_container->ext_dipoles[i]->p_y;
    dz = this->p_z - this->parent_container->ext_dipoles[i]->p_z;
    d = sqrt(dx * dx + dy * dy + dz * dz);

    m_1_x = this->parent_container->ext_dipoles[i]->mag_m_x;
    m_1_y = this->parent_container->ext_dipoles[i]->mag_m_y;
    m_1_z = this->parent_container->ext_dipoles[i]->mag_m_z;
    m_1_r = m_1_x * dx + m_1_y * dy + m_1_z * dz;

    B_x += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dx / (d * d * d * d * d) - m_1_x / (d * d * d));
    B_y += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dy / (d * d * d * d * d) - m_1_y / (d * d * d));
    B_z += permeability / (4.0 * M_PIl) * (3.0 * m_1_r * dz / (d * d * d * d * d) - m_1_z / (d * d * d));

  }
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

  return result;

}
