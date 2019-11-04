#include "container.h"
#include "magnetic_dipole.cpp"

Container :: Container() {

  this->dt = 0.01;
  this->output_file.open("magnetic_gas_output.txt");
  this->is_journal = 0;
  this->save_energy = 0;
  this->save_thermodynamics = 0;
  this->save_cloud_chamber = 0;

  this->cylindrical_impulse = 0.0;
  this->surface_impulse = 0.0;
  this->xy1_impulse = 0.0;
  this->ext_B_x = 0.0;
  this->ext_B_y = 0.0;
  this->ext_B_z = 0.0;

}

Container :: Container(double my_dt) {

  this->dt = my_dt;
  this->output_file.open("magnetic_gas_output.txt");
  this->is_journal = 0;
  this->save_energy = 0;
  this->save_thermodynamics = 0;
  this->save_cloud_chamber = 0;

  this->cylindrical_impulse = 0.0;
  this->surface_impulse = 0.0;
  this->xy1_impulse = 0.0;
  this->ext_B_x = 0.0;
  this->ext_B_y = 0.0;
  this->ext_B_z = 0.0;

}

Container :: Container(double my_dt, std::string journal_name) {

  this->dt = my_dt;
  this->output_file.open("magnetic_gas_output.txt");
  this->is_journal = 1;
  this->my_journal = new Journal(journal_name);
  if(this->is_journal) this->my_journal->JournalW("Container created.");
  this->save_energy = 0;
  this->save_thermodynamics = 0;
  this->save_cloud_chamber = 0;

  this->cylindrical_impulse = 0.0;
  this->surface_impulse = 0.0;
  this->xy1_impulse = 0.0;
  this->ext_B_x = 0.0;
  this->ext_B_y = 0.0;
  this->ext_B_z = 0.0;

}

Container :: ~Container() {

  int i;
  if(this->is_journal) this->my_journal->JournalIn("Delete container...");
  if(this->save_energy) {
    if(this->is_journal) this->my_journal->JournalW("Closing energy output file...");
    this->energy_output_file.close();
  }
  if(this->save_thermodynamics) {
    if(this->is_journal) this->my_journal->JournalW("Closing thermodynamics output files...");
    this->pressure_output_file.close();
    this->temperature_output_file.close();
  }
  if(this->save_cloud_chamber) {
    if(this->is_journal) this->my_journal->JournalW("Closing cloud chamber output file...");
    this->cloud_chamber_output_file.close();
  }

  if(this->is_journal) this->my_journal->JournalW("Closing default output file...");
  this->output_file.close();
  if(this->is_journal) this->my_journal->JournalW("Deleting particles...");
  for(i = 0; i < this->particles.size(); i++) delete this->particles[i];
  if(this->is_journal) this->my_journal->JournalOut();
  delete this->my_journal;

}

void Container :: cuboid_shape(double my_d_x, double my_d_y, double my_d_z) {

  if(this->is_journal) this->my_journal->JournalIn("Setting shape to cuboid...");
  this->d_x = my_d_x;
  this->d_y = my_d_y;
  this->d_z = my_d_z;
  this->shape = 0;
  this->surface = 2 * this->d_x * this->d_y + 2 * this->d_x * this->d_z + 2 * this->d_y * this->d_z;

  if(this->is_journal) this->my_journal->JournalW("Initializing random generators...");

  this->rand_x = std::uniform_real_distribution<double>(0.0, this->d_x);
  this->rand_y = std::uniform_real_distribution<double>(0.0, this->d_y);
  this->rand_z = std::uniform_real_distribution<double>(0.0, this->d_z);

  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: cylindrical_shape(double my_radius, double my_length) {

  if(this->is_journal) this->my_journal->JournalIn("Setting shape to cylindrical...");

  this->radius = my_radius;
  this->length = my_length;
  this->shape = 1;
  this->surface = 2.0 * M_PIl * this->radius * this->length + 2.0 * M_PIl * this->radius * this->radius;

  if(this->is_journal) this->my_journal->JournalW("Initializing random generators...");

  this->rand_x = std::uniform_real_distribution<double>(0.0, this->radius);
  this->rand_z = std::uniform_real_distribution<double>(0.0, this->length);

  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: add_external_dipole(double my_p_x, double my_p_y, double my_p_z, double my_mag_m_x, double my_mag_m_y, double my_mag_m_z) {

  ExternalMagneticDipole* new_ext_mag_dipole = new ExternalMagneticDipole;
  new_ext_mag_dipole->p_x = my_p_x;
  new_ext_mag_dipole->p_y = my_p_y;
  new_ext_mag_dipole->p_z = my_p_z;
  new_ext_mag_dipole->mag_m_x = my_mag_m_x;
  new_ext_mag_dipole->mag_m_y = my_mag_m_y;
  new_ext_mag_dipole->mag_m_z = my_mag_m_z;
  this->ext_dipoles.push_back(new_ext_mag_dipole);

}

void Container :: external_magnetic_field(double m_B_x, double m_B_y, double m_B_z) {

  this->ext_B_x = m_B_x;
  this->ext_B_y = m_B_y;
  this->ext_B_z = m_B_z;

}

void Container :: add_determined_particle(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding determined particle...");
  this->particles.push_back(new Particle(this, m_p_x, m_p_y, m_p_z, m_v_x, m_v_y, m_v_z, m_m, m_r));
  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: add_particle(double m_v, double m_m, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding random particle...");
  this->particles.push_back(new Particle(this, 0,0,0,0,0,0, m_m, m_r));
  if(this->is_journal) this->my_journal->JournalW("Randomizing my particle");
  this->particles[this->particles.size()-1]->randomize_particle(m_v);
  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: add_determined_electric_particle(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_q, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding determined electric particle...");
  this->particles.push_back(new ElectricParticle(this, m_p_x, m_p_y, m_p_z, m_v_x, m_v_y, m_v_z, m_m, m_q, m_r));
  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: add_electric_particle(double m_v, double m_m, double m_q, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding random electric particle...");
  this->particles.push_back(new ElectricParticle(this, 0,0,0,0,0,0, m_m, m_q, m_r));
  if(this->is_journal) this->my_journal->JournalW("Randomizing my particle");
  this->particles[this->particles.size()-1]->randomize_particle(m_v);
  if(this->is_journal) this->my_journal->JournalOut();

}

void Container :: add_determined_magnetic_dipole(double m_p_x, double m_p_y, double m_p_z, double m_v_x, double m_v_y, double m_v_z, double m_m, double m_q, double m_mag_x, double m_mag_y, double m_mag_z, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding determined magnetic dipole...");
  this->particles.push_back(new MagneticDipole(this, m_p_x, m_p_y, m_p_z, m_v_x, m_v_y, m_v_z, m_m, m_q, m_mag_x, m_mag_y, m_mag_z, m_r));
  if(this->is_journal) this->my_journal->JournalOut();

}


void Container :: add_magnetic_dipole(double m_v, double m_m, double m_q, double m_mag, double m_r) {

  if(this->is_journal) this->my_journal->JournalIn("Adding random magnetic dipole...");
  this->particles.push_back(new MagneticDipole(this, 0,0,0,0,0,0, m_m, m_q, 0,0,0, m_r));
  if(this->is_journal) this->my_journal->JournalW("Randomizing my particle");
  ((MagneticDipole*) this->particles[this->particles.size()-1])->randomize_particle(m_v, m_mag);
  if(this->is_journal) this->my_journal->JournalOut();

}


void Container :: tick() {

  int i, j;
  double a, dis, my_x, my_y, my_v_x, my_v_y;
  std::vector<double> response;
  if(this->is_journal) this->my_journal->JournalIn("Tick");
  this->kinetic_energy = 0;
  if(this->is_journal) this->my_journal->JournalW("Calculating kinetic energy...");
  for(i = 0; i < this->particles.size(); i++) {
    this->kinetic_energy += 0.5 * (this->particles[i]->v_x * this->particles[i]->v_x + this->particles[i]->v_y * this->particles[i]->v_y + this->particles[i]->v_z * this->particles[i]->v_z) * this->particles[i]->m;
    this->particles[i]->next_a_x = 0;
    this->particles[i]->next_a_y = 0;
    this->particles[i]->next_a_z = 0;
    this->particles[i]->eps_x = 0;
    this->particles[i]->eps_y = 0;
    this->particles[i]->eps_z = 0;
  }
  if(this->is_journal) this->my_journal->JournalW("Calculating interactions...");
  for(i = 0; i < this->particles.size(); i++) {
    if(this->is_journal) this->my_journal->JournalIn("This particle...");
    for(j = 0; j < this->particles.size(); j++) {
      if(i == j) continue;
      if(this->is_journal) this->my_journal->JournalW("Getting response...");
      if(this->is_journal) this->my_journal->JournalW(std::to_string(j));
      response = this->particles[i]->force_imposed_by(this->particles[j]);
      if(this->is_journal) this->my_journal->JournalW("Updating accelaration...");
      this->particles[i]->next_a_x += response[0] / this->particles[i]->m;
      this->particles[i]->next_a_y += response[1] / this->particles[i]->m;
      this->particles[i]->next_a_z += response[2] / this->particles[i]->m;
      this->particles[i]->eps_x    += response[3] / this->particles[i]->inertia;
      this->particles[i]->eps_y    += response[4] / this->particles[i]->inertia;
      this->particles[i]->eps_z    += response[5] / this->particles[i]->inertia;
    }
    if(this->is_journal) this->my_journal->JournalW("Container: Getting response...");
    response = this->particles[i]->force_imposed_by_container();
    if(this->is_journal) this->my_journal->JournalW("Container: Updating acceleration...");
    this->particles[i]->next_a_x += response[0] / this->particles[i]->m;
    this->particles[i]->next_a_y += response[1] / this->particles[i]->m;
    this->particles[i]->next_a_z += response[2] / this->particles[i]->m;
    this->particles[i]->eps_x    += response[3] / this->particles[i]->inertia;
    this->particles[i]->eps_y    += response[4] / this->particles[i]->inertia;
    this->particles[i]->eps_z    += response[5] / this->particles[i]->inertia;
    if(this->is_journal) this->my_journal->JournalOut();
  }
  if(this->is_journal) this->my_journal->JournalW("Moving particles...");
  for(i = 0; i < this->particles.size(); i++) {
    this->particles[i]->a_x = this->particles[i]->next_a_x;
    this->particles[i]->a_y = this->particles[i]->next_a_y;
    this->particles[i]->a_z = this->particles[i]->next_a_z;
    this->particles[i]->v_x += this->particles[i]->a_x * this->dt;
    this->particles[i]->v_y += this->particles[i]->a_y * this->dt;
    this->particles[i]->v_z += this->particles[i]->a_z * this->dt;
    this->particles[i]->p_x += this->particles[i]->v_x * this->dt + 0.5 * this->particles[i]->a_x * this->dt * this->dt;
    this->particles[i]->p_y += this->particles[i]->v_y * this->dt + 0.5 * this->particles[i]->a_y * this->dt * this->dt;
    this->particles[i]->p_z += this->particles[i]->v_z * this->dt + 0.5 * this->particles[i]->a_z * this->dt * this->dt;
    this->particles[i]->rotate_self();
    switch(this->shape) {
      case 0:
        if(this->particles[i]->p_x < 0) {
          this->particles[i]->p_x *= -1.0;
          this->particles[i]->v_x *= -1.0;
          this->surface_impulse += this->particles[i]->v_x * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_x > this->d_x) {
          this->particles[i]->p_x = this->d_x - (this->particles[i]->p_x - this->d_x);
          this->particles[i]->v_x *= -1.0;
          this->surface_impulse -= this->particles[i]->v_x * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_y < 0) {
          this->particles[i]->p_y *= -1.0;
          this->particles[i]->v_y *= -1.0;
          this->surface_impulse += this->particles[i]->v_y * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_y > this->d_y) {
          this->particles[i]->p_y = this->d_y - (this->particles[i]->p_y - this->d_y);
          this->particles[i]->v_y *= -1.0;
          this->surface_impulse -= this->particles[i]->v_y * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_z < 0) {
          this->particles[i]->p_z *= -1.0;
          this->particles[i]->v_z *= -1.0;
          this->xy1_impulse += this->particles[i]->v_z * 2.0 * this->particles[i]->m; //here be dragons
          this->surface_impulse += this->particles[i]->v_z * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_z > this->d_z) {
          this->particles[i]->p_z = this->d_z - (this->particles[i]->p_z - this->d_z);
          this->particles[i]->v_z *= -1.0;
          this->surface_impulse -= this->particles[i]->v_z * 2.0 * this->particles[i]->m;
        }
        break;
      case 1:
        a = (this->particles[i]->p_x * this->particles[i]->p_x) + (this->particles[i]->p_y * this->particles[i]->p_y);
        if(this->particles[i]->p_z < 0) {
          this->particles[i]->p_z *= -1.0;
          this->particles[i]->v_z *= -1.0;
          this->xy1_impulse += this->particles[i]->v_z * 2.0 * this->particles[i]->m; //here be dragons
          this->surface_impulse += this->particles[i]->v_z * 2.0 * this->particles[i]->m;
        }
        if(this->particles[i]->p_z > this->length) {
          this->particles[i]->p_z = this->length - (this->particles[i]->p_z - this->length);
          this->particles[i]->v_z *= -1.0;
          this->surface_impulse -= this->particles[i]->v_z * 2.0 * this->particles[i]->m;
        }
        if(a >= this->radius * this->radius) {
          dis = sqrt(a);
          my_x = this->particles[i]->p_x * this->radius / dis;
          my_y = this->particles[i]->p_y * this->radius / dis;
          my_v_x = this->particles[i]->v_x;
          my_v_y = this->particles[i]->v_y;
          this->particles[i]->p_x = this->particles[i]->p_x * (2.0 - dis / this->radius);
          this->particles[i]->p_y = this->particles[i]->p_y * (2.0 - dis / this->radius);
          this->particles[i]->v_x = (my_y * (my_v_x * my_y - my_v_y * my_x) - my_x * (my_v_x * my_x + my_v_y * my_y)) / (this->radius * this->radius);
          this->particles[i]->v_y = (my_x * (my_v_y * my_x - my_v_x * my_y) - my_y * (my_v_x * my_x + my_v_y * my_y)) / (this->radius * this->radius);
          this->surface_impulse += (1.0 / this->radius) * (my_v_x * my_x + my_v_y * my_y) * 2.0 * this->particles[i]->m;
          this->cylindrical_impulse += (1.0 / this->radius) * (my_v_x * my_x + my_v_y * my_y) * 2.0 * this->particles[i]->m;
        }
        break;
    }
    if(this->save_cloud_chamber) {
      this->cloud_chamber_output_file << " " << this->particles[i]->p_x << " " << this->particles[i]->p_y;
    }
  }
  if(this->is_journal) this->my_journal->JournalW("Calculating temperature...");
  this->temperature = this->kinetic_energy / this->particles.size() * 2.0 / 3.0 / boltzmann_constant;
  if(this->is_journal) this->my_journal->JournalW("Updating cloud chamber...");
  if(this->save_cloud_chamber) {
    this->cloud_chamber_output_file << std::endl;
  }
  this->t += this->dt;
  if(this->is_journal) this->my_journal->JournalOut();
}

void Container :: Simulate(double max_t, double check_time, double data_points) {

  double checkpoint = check_time;
  double delta_t = max_t / data_points;
  double data_checkpoint = delta_t;
  std::cout << "Simulation started" << std::endl;


  while(this->t < max_t) {

    this->tick();
    if (this->t > checkpoint) {
      std::cout << "t = " << checkpoint << "s" << std::endl;
      checkpoint += check_time;
    }

    if (this->t > data_checkpoint) {
      if(this->save_energy) this->energy_output_file << this->t << " " << this->kinetic_energy << std::endl;
      if(this->save_thermodynamics) {
        this->temperature_output_file << this->t << " " << this->temperature << std::endl;
        this->pressure_output_file    << this->t << " " << this->surface_impulse / this->t / this->surface << std::endl;
      }
      data_checkpoint += delta_t;
    }

  }
  std::cout << "Simulation ended" << std::endl;

}

void Container :: Simulate(double max_t, double check_time, double data_points, double omega, long int magnet_distance, long int permanent_magnet_force_coef) {

  double checkpoint = check_time;
  double delta_t = max_t / data_points;
  double data_checkpoint = delta_t;
  double positions_z;
  std::vector<double> juj;
  double f_mag;
  int i, no_in_row = 0;
  std::vector<double> smer, smer_dipolu;
  std::vector<double> pozicia_dipolu(3,0.0);
  std::ofstream default_file_output;
  this->delta_alpha = 0.0;
  default_file_output.open("smer_mag_momentu.dat");
  pozicia_dipolu[0] = ((float) magnet_distance) / 10.0;
  for(positions_z = 0.0; positions_z <= this->length; positions_z += 0.1) {
    this->add_external_dipole(0.0, 0.0, positions_z, 0.0, 0.0, 0.0);
    this->add_external_dipole(0.0, 0.0, positions_z, 0.0, 0.0, 0.0);
    no_in_row++;
  }

  while(this->t < max_t) {

    for(i = 0; i < no_in_row; i++) {
      this->ext_dipoles[i]->p_x = pozicia_dipolu[0];
      this->ext_dipoles[i]->p_y = pozicia_dipolu[1];
      this->ext_dipoles[i]->mag_m_x = pozicia_dipolu[0] * ((float) permanent_magnet_force_coef);
      this->ext_dipoles[i]->mag_m_y = pozicia_dipolu[1] * ((float) permanent_magnet_force_coef);
    }
    for(i = no_in_row; i < this->ext_dipoles.size(); i++) {
      this->ext_dipoles[i]->p_x = -1 * pozicia_dipolu[0];
      this->ext_dipoles[i]->p_y = -1 * pozicia_dipolu[1];
      this->ext_dipoles[i]->mag_m_x = pozicia_dipolu[0] * ((float) permanent_magnet_force_coef);
      this->ext_dipoles[i]->mag_m_y = pozicia_dipolu[1] * ((float) permanent_magnet_force_coef);
    }

    this->tick();
    this->delta_alpha += this->particles[0]->omega_z * this->dt + 0.5 * this->particles[0]->eps_z * this->dt * this->dt;
    if (this->t > checkpoint) {
      std::cout << "t = " << checkpoint << "s" << std::endl;
      checkpoint += check_time;
    }
    smer = vector_direction(((MagneticDipole*) this->particles[0])->mag_m_x, ((MagneticDipole*) this->particles[0])->mag_m_y, ((MagneticDipole*) this->particles[0])->mag_m_z);
    smer_dipolu = vector_direction(pozicia_dipolu[0], pozicia_dipolu[1], pozicia_dipolu[2]);
    pozicia_dipolu = rotate_vector(pozicia_dipolu[0], pozicia_dipolu[1], pozicia_dipolu[2], 0, 0, omega * this->dt);

    if (this->t > data_checkpoint) {
      if(this->save_energy) this->energy_output_file << this->t << " " << this->kinetic_energy << std::endl;
      if(this->save_thermodynamics) {
        this->temperature_output_file << this->t << " " << this->temperature << std::endl;
        this->pressure_output_file    << this->t << " " << this->surface_impulse / this->t / this->surface << std::endl;
      }
      data_checkpoint += delta_t;
    }

  }
  default_file_output.close();
  std::cout << "Simulation ended" << std::endl;

}

void Container :: print_particles() {

  int i;
  for (i = 0; i < this->particles.size(); i++) this->output_file << "Particle " << i << ":" << std::endl << "  p_x = " << this->particles[i]->p_x << std::endl << "  p_y = " << this->particles[i]->p_y << std::endl << "  p_z = " << this->particles[i]->p_z << std::endl << "  v_x = " << this->particles[i]->v_x << std::endl << "  v_y = " << this->particles[i]->v_y << std::endl << "  v_z = " << this->particles[i]->v_z << std::endl << "  v = " << sqrt(this->particles[i]->v_x * this->particles[i]->v_x + this->particles[i]->v_y * this->particles[i]->v_y + this->particles[i]->v_z * this->particles[i]->v_z) << std::endl;

}

void Container :: SaveVelocityDistribution() {
  double my_v, max_v = 0.0;
  int i, d_size;
  for(i = 0; i < this->particles.size(); i++) {
    my_v = sqrt(this->particles[i]->v_x * this->particles[i]->v_x + this->particles[i]->v_y * this->particles[i]->v_y + this->particles[i]->v_z * this->particles[i]->v_z);
    if(my_v > max_v) max_v = my_v;
  }
  d_size = ceil(max_v);
  std::vector<int> v_dis(d_size, 0);
  for(i = 0; i < this->particles.size(); i++) {
    my_v = sqrt(this->particles[i]->v_x * this->particles[i]->v_x + this->particles[i]->v_y * this->particles[i]->v_y + this->particles[i]->v_z * this->particles[i]->v_z);
    v_dis[floor(my_v)]++;
  }
  this->velocity_distribution_output_file.open("magnetic_gas_velocity_distribution.dat");
  for(i = 0; i < d_size; i++) this->velocity_distribution_output_file << i << " " << v_dis[i] << std::endl;
  this->velocity_distribution_output_file.close();
}

void Container :: SaveEnergy() {
  this->save_energy = 1;
  this->energy_output_file.open("magnetic_gas_energy.dat");
}

void Container :: SaveThermodynamics() {
  this->save_thermodynamics = 1;
  this->pressure_output_file.open("magnetic_gas_pressure.dat");
  this->temperature_output_file.open("magnetic_gas_temperature.dat");
}

void Container :: CloudChamber() {
  this->save_cloud_chamber = 1;
  this->cloud_chamber_output_file.open("magnetic_gas_cloud_chamber.dat");
}
