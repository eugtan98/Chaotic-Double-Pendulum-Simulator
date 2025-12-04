import pendulum_model
import chaos_tools
import plotter

# Set default physical parameters
deflt_m1 = 1.0 #kg
deflt_m2 = 1.0 #kg
deflt_l1 = 1.0 #m
deflt_l2 = 1.0 #m
deflt_g = 9.81 #kgm/s^2
#Set up defaults in a tuple
defaults = (deflt_m1, deflt_m2, deflt_l1, deflt_l2, deflt_g)

def get_parameters():
    print('Use default parameters or custom?')
    choice_para = input('Enter d for default and c for custom:')

    if choice_para == 'd':
        m1, m2, l1, l2 ,g = defaults
    else:
        m1 = input("Enter m1:")
        m2 = input("Enter m2:")
        l1 = input("Enter l1:")
        l2 = input("Enter l2")
        g = input("Enter g")
    para = {"m1": m1, "m2":m2, "l1": l1, "l2":l2, "g":g}
    return para

def get_initial_state():
    print("Enter initial angles (in radians) and angular velocities.")
    theta1_0 = input("Enter first initial angle:")
    omega1_0 = input("Enter first angular velocity: ")
    theta2_0 = input("Enter second initial angle:")
    omega2_0 = input("Enter second angular velocity: ")
    initial_state = [theta1_0, omega1_0, theta2_0, omega2_0]
    return initial_state

def get_time_settings():
    t_max = input("Enter total simulation time:")
    dt = input("Enter time step:")
    return t_max, dt

def run_single_simulation():
    print("=== Single Trajectory Simulation ===")
    #Get parameters, initial conditions, time settings
    params = get_parameters()
    initial_state = get_initial_state()
    t_max, dt = get_time_settings()

    #Run simulation
    (t_array, theta1_array, omega1_array, theta2_array, omega2_array) = pendulum_model.simulate(initial_state, t_max, dt, params)

    #Print brief numeric summary
    print("Simulation finished.")
    print("Number of steps =", len(t_array))
    print(f'Initial theta1 = {theta1_array[0]}, Final theta1 = {theta1_array[-1]}')

    #Optional: Print max angle
    max_theta1 = max(theta1_array)
    min_theta1 = min(theta1_array)
    print("Max theta1 =", max_theta1, ", Min theta1 =", min_theta1)

    #Plots / animation
    ans = input("Animate the pendulum? (y/n): ").strip().lower()
    if ans == 'y':
        plotter.animate_pendulum(theta1_array, theta2_array, params)

    ans = input("Plot angles vs time? (y/n): ").strip().lower()
    if ans == 'y':
        plotter.plot_angles_vs_time(t_array, theta1_array, theta2_array)

def run_chaos_demo():
    #Get parameters, initial conditions, time settings
    print("=== Chaos Demo: Two Nearby Trajectories ===")
    params = get_parameters()
    initial_state = get_initial_state()
    t_max, dt = get_time_settings()

    #Ask for small perturbation
    delta_theta1 = input("Enter small change in theta1 for the second trajectory (e.g. 0.001 rad):")

    #Run two simulations                                                                                 
    (t_array, theta1_a, omega1_a, theta2_a, omega2_a, 
     theta1_b, omega1_b, theta2_b, omega2_b) = chaos_tools.simulate_two_trajectories(initial_state, 
                                                                                     delta_theta1, t_max, dt, params)
    #Compute separation in state space
    separation_array = chaos_tools.compute_separation(theta1_a, theta2_a, theta1_b, theta2_b)

    #Print info
    print(f'Initial angle difference = {delta_theta1}')
    print(f'Final seperation in angle space = {separation_array[-1]}')

    #Plot/ animations
    ans = input("Plot separation vs time? (y/n)").strip().lower()
    if ans == 'y':
        plotter.plot_separation_vs_time(t_array, separation_array)
    ans = input("Animate both pendulums together? (y/n)")
    if ans == 'y':
        plotter.animate_two_pendulums(theta1_a, theta2_a, theta1_b, theta2_b, params)

def run_phase_space_plot():
    #Get parameters, initial conditions, time settings
    print("=== Phase-Space Plot (theta vs omega) ===")
    params = get_parameters()
    initial_state = get_initial_state
    t_max, dt = get_time_settings()

    #Simulate once
    (t_array, theta1_array, omega1_array, theta2_array, omega2_array) = pendulum_model.simulate(initial_state, t_max, dt, params)

    #Ask which pendulum to plot
    choice_plot = input("Which pendulum to plot in phase-space? (1 or 2)")
    if choice_plot == '1':
        plotter.plot_phase_space(theta1_array, omega1_array)
    else:
        plotter.plot_phase_space(theta2_array, omega2_array)

def main():
    print("Welcome")
    print("A Python-based simulator that solves and " \
    "visualizes the chaotic motion of a double pendulum, highlighting extreme sensitivity to initial conditions.")

    

    while True:
        print("1 for Single trajectory simulation & animation" \
        "2 for Chaos demo: two nearby trajectories" \
        "3 for Phase-space plot (theta vs omega)" \
        "4 for Exit")
        choice = input('Enter choice (1,2,3,4):')
        if choice == '1':
            run_single_simulation()
        elif choice == '2':
            run_chaos_demo()
        elif choice == '3':
            run_phase_space_plot()
        elif choice == '4':
            print("Goodbye!")
            break
    

if __name__ == '__main__':
    main()