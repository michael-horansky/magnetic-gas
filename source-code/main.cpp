#include "container.cpp"
#include <omp.h>

int main(int argc, char *argv[])
{

    long int i, n_part, n_sim, step_coef, magnet_distance, dt_coef, max_vel, permanent_magnet_force_coef, my_inertia, vlastny_moment, radius_coef, omega_coef;

    if (argc > 1 )
        n_sim = atol(argv[1]); // Reading number of particles
    else
        n_sim = 300;

    if (argc > 2 )
        n_part = atol(argv[2]); // Reading number of particles
    else
        n_part = 300;

    if (argc > 3 )
        step_coef = atol(argv[3]); // Reading step coef
    else
        step_coef = 2;

    if (argc > 4 )
        magnet_distance = atol(argv[4]); // Reading magnet distance
    else
        magnet_distance = 10;

    if (argc > 5 )
        dt_coef = atol(argv[5]); // Reading dt coefficient (0.0001 * c)
    else
        dt_coef = 5;

    if (argc > 6 )
        max_vel = atol(argv[6]); // Reading maximal velocity
    else
        max_vel = 0;

    if (argc > 7 )
        permanent_magnet_force_coef = atol(argv[7]); // Reading permanent magnet force coef
    else
        permanent_magnet_force_coef = 2000000;

    if (argc > 8 )
        my_inertia = atol(argv[8]); // Reading particle inertia
    else
        my_inertia = 10;

    if (argc > 9 )
        vlastny_moment = atol(argv[9]); // Reading particle inertia
    else
        vlastny_moment = 10;

    if (argc > 10 )
        radius_coef = atol(argv[10]); // Reading particle inertia
    else
        radius_coef = 10;

    if (argc > 11 )
        omega_coef = atol(argv[11]); // Reading particle inertia
    else
        omega_coef = 12;

    std::cout << "Number of simulations:  " << n_sim << std::endl;
    std::cout << "Number of particles:  " << n_part << std::endl;
    std::cout << "Step:  " << step_coef << std::endl;
    std::cout << "Magnet distance:  " << 0.1 * ((float) magnet_distance) << std::endl;
    std::cout << "dt:  " << 0.00001 * ((float) dt_coef) << std::endl;
    std::cout << "Maximal velocity:  " << (float) max_vel << std::endl;
    std::cout << "Lambda:  " << 1.0 * ((float) magnet_distance) * (float) permanent_magnet_force_coef << std::endl;
    std::cout << "Particle inertia:  " << 1e-31 * ((float) my_inertia) << std::endl;
    std::cout << "Particle magnetic moment:  " << ((float) vlastny_moment) * 1.5e-28 << std::endl;
    std::cout << "Cylinder radius:  " << 0.005 * ((float) radius_coef) << std::endl;
    std::cout << "Omega (radial velocity):  " << ((float) omega_coef) * 5.0 << std::endl;

    std::ofstream default_file_output;
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 = d.count();
    random_engine.seed(seed2);
    default_file_output.open("omega_of_spikes_phase_shift_real.dat");
    double actual_inertia;

    double average_omega;
    double average_phase_shift;
    double average_x, average_y;

    #pragma omp parallel
    {

      Container *la;
      int N;

      #pragma omp for schedule(static,1)
      for(i = 1; i <= n_sim; i++) {
        std::cout << "Simulation no. " << i << " started." << std::endl;
        la = new Container(0.00001 * ((float) dt_coef));
        la->cylindrical_shape(0.005 * ((float) radius_coef), 5.0);

        actual_inertia = 1e-32 * ((float) my_inertia) * (((float) i) * ((float) step_coef));

        for(N = 0; N < n_part; N++) la->add_magnetic_dipole((float) max_vel, 0.004 / avogadro_constant, elementary_charge * 5.0e-06, 0.0, actual_inertia);
        for(N = 0; N < n_part; N++) ((MagneticDipole*)(la->particles[N]))->mag_m_x = ((float) vlastny_moment) * 1.5e-28;

        la->Simulate(3.0, 0.5, 500, ((float) omega_coef) * 5.0, magnet_distance, permanent_magnet_force_coef);

        #pragma omp critical
        //std::cout << la->particles[0]->v_x << std::endl;

        average_omega = 0.0;
        average_phase_shift = 0.0;
        average_x = 0.0;
        average_y = 0.0;
        for(N = 0; N < n_part; N++) {
          average_omega += la->particles[N]->omega_z;
          average_x += ((MagneticDipole*) la->particles[N])->p_x;
          average_y += ((MagneticDipole*) la->particles[N])->p_y;
        }
        average_omega /= n_part;
        average_phase_shift = atan(la->ext_dipoles[0]->mag_m_y / la->ext_dipoles[0]->mag_m_x) - atan(average_y / average_x);
        std::cout << average_omega <<  " " << average_phase_shift << std::endl;

        default_file_output << /*1e-32 * ((float) my_inertia) * (((float) i) * ((float) step_coef))*/actual_inertia << " " << /*la->delta_alpha / la->t*/average_omega << " " << average_phase_shift << " " << la->cylindrical_impulse / la->t / (2.0 * M_PIl * la->radius * la->length) << " " << la->temperature << std::endl;
      }
      delete la;
    }
    default_file_output.close();
    std::cout << "Simulation ended " << avogadro_constant << std::endl;
}
