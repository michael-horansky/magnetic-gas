#include "particles.h"

Particle :: Particle() {


}

Particle :: Particle(Container *my_parent_container, double start_pos_x, double start_pos_y, double start_pos_z, double start_vel_x, double start_vel_y, double start_vel_z, double start_mass, double start_inertia) {

  this->parent_container = my_parent_container;
  if(this->parent_container->is_journal) this->parent_container->my_journal->JournalIn("Initializing particle...");
  this->p_x = start_pos_x;
  this->p_y = start_pos_y;
  this->p_z = start_pos_z;
  this->v_x = start_vel_x;
  this->v_y = start_vel_y;
  this->v_z = start_vel_z;
  this->a_x = 0;
  this->a_y = 0;
  this->a_z = 0;
  this->next_a_x = 0;
  this->next_a_y = 0;
  this->next_a_z = 0;

  this->m = start_mass;
  this->r = 0.0000001;
  this->inertia = start_inertia;

  this->omega_x = 0;
  this->omega_y = 0;
  this->omega_z = 0;
  this->eps_x = 0;
  this->eps_y = 0;
  this->eps_z = 0;
  if(this->parent_container->is_journal) this->parent_container->my_journal->JournalOut();

}

Particle :: ~Particle() {

}

void Particle :: randomize_particle(double max_vel) {

  if(this->parent_container->is_journal) this->parent_container->my_journal->JournalIn("Randomizing particle...");
  double phi;
  if(this->parent_container->shape == 0) {
    if(this->parent_container->is_journal) this->parent_container->my_journal->JournalW("Randomizing cuboid position...");
    this->p_x = this->parent_container->rand_x(random_engine);
    this->p_y = this->parent_container->rand_y(random_engine);
    this->p_z = this->parent_container->rand_z(random_engine);
  } else if(this->parent_container->shape == 1) {
    if(this->parent_container->is_journal) this->parent_container->my_journal->JournalW("Randomizing cylindrical position...");
    double radius_distance = this->parent_container->rand_x(random_engine);
    phi = rand_phi(random_engine);
    this->p_x = radius_distance * sin(phi);
    this->p_y = radius_distance * cos(phi);
    this->p_z = this->parent_container->rand_z(random_engine);
  }
  if(this->parent_container->is_journal) this->parent_container->my_journal->JournalW("Randomizing velocity...");
  double velocity = rand_def(random_engine) * max_vel;
  phi = rand_phi(random_engine);
  double theta = rand_theta(random_engine);
  this->v_x = velocity * sin(theta) * cos(phi);
  this->v_y = velocity * sin(theta) * sin(phi);
  this->v_z = velocity * cos(theta);
  this->a_x = 0;
  this->a_y = 0;
  this->a_z = 0;
  this->next_a_x = 0;
  this->next_a_y = 0;
  this->next_a_z = 0;
  if(this->parent_container->is_journal) this->parent_container->my_journal->JournalOut();

}

