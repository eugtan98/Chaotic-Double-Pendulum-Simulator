# 2025/12/3
Ask ChatGPT to provides a list of scientific/engineering projects that are possible to be done by python coding.
**Chaotic Double Pendulum Simulator & Chaos Visualizer** was chosen

Ask ChatGPT to introduce the concept of double pendulum, then teach me to do the project by providing me a **pseudocode** and finish the coding program on my own.

Ask ChatGPT when the programming process is stuck, where the code is difficult to write it out.
**main.py** coding program was done

# 2025/12/4
Finished **pendulum_model.py**

# 2025/12/8
Finished **chaos_tools.py**, Writen some of the code for **plotter.py**

# 2025/12/10
Finished up **plotter.py** (angles-time plot, phase-space plot, sepration-time plot)

Run the program to DEBUG:
**1. Fixed Runtime Error in main.py**
- Issue: The function 
"get_initial_state" function was assigned to a variable instead of being called, causing a ZeroDivisionError later because "initial_state" was a function object, not a list of values.

- Fix: Added parentheses to call the function: initial_state = get_initial_state().

**2. Corrected Physics in pendulum_model.py**
- Issue: The equations of motion for the second pendulum (domega2_dt) contained copy-paste errors: den2 was calculated using l1 instead of l2.

num2 used l1 instead of l2 in the velocity-dependent term and had incorrect signs.

- Fix: Updated den2 and num2 to match standard Lagrangian mechanics for a double pendulum.

**3. Verified plotter.py**
- Action: Temporarily modified  plotter.py to save plots to files for automated verification.

- Result: Confirmed that plotting functions work correctly.

- Cleanup: Reverted all changes to plotter.py to restore the original interactive plt.show() behavior.

- Verification Results
    
    - Execution: main.py runs without errors.
    
    - Physics: A test simulation run confirmed that the corrected equations produce valid numerical output without crashing.

# 2025/12/11
**1. The Reported Error: TypeError: 'Line2D' object is not iterable**
- Location: plotter.py, inside animate_pendulum function.
- Cause: Matplotlib's FuncAnimation with blit=True requires the updateand init functions to return an iterable of artists (like a tuple or list). The code was returning a single Line2D object.

- Fix: Changed the return statements to return a tuple.

(Before) return rod_line

(After) return rod_line, [note the comma]

**2. Mathematical Error in pendulum_model.py**
- Location: pendulum_model.py, inside rk4_step function.

- Cause: Noticed an error in the Runge-Kutta 4 (RK4) implementation. The calculation for state_k4 was using 0.5 * dt instead of dt. This would have caused simulation inaccuracies.

- Fix:

(Before) state_k4 = [state[i] + 0.5 * dt *k3[i] for i in range(4)]

(After) state_k4 = [state[i] + dt *k3[i] for i in range(4)]

**3. Missing Arguments in plotter.py**
- Location: plotter.py, inside animate_two_pendulums function.

- Cause: The function animate_two_pendulums calls compute_positions, but compute_positions expects 4 arguments (theta1, theta2, l1, l2). The code was only passing 2 arguments (theta1, theta2), which would have caused a runtime error when running the Chaos Demo.

- Fix: Passed the missing l1 and l2 arguments.

(Before) (x1_a, y1_a, x2_a, y2_a) = compute_positions(theta1_a, theta2_a)
(After) (x1_a, y1_a, x2_a, y2_a) = compute_positions(theta1_a, theta2_a, l1, l2)

**4. Minor Fixes**
chaos_tools.py: Fixed a typo in the argument name intiial_state -> initial_state.

**5. Verification Process**

Steps Verified:

Single Simulation logic: Verified that **run_single_simulation** calls **pendulum_model.simulate** and **plotter.animate_pendulum** correctly.

Result: PASSED. (The TypeError is gone).

Chaos Demo logic: Verified that **run_chaos_demo** calls **chaos_tools** and **plotter.animate_two_pendulums**.

Result: PASSED. (Caught the missing l1, l2 bug).