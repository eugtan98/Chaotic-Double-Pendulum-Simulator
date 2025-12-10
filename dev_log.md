2025/12/3
Ask ChatGPT to provides a list of scientific/engineering projects that are possible to be done by python coding.
Chaotic Double Pendulum Simulator & Chaos Visualizer was chosen
Ask ChatGPT to introduce the concept of double pendulum, then teach me to do the project by providing me a pseudocode and finish the coding program on my own.
Ask ChatGPT when the programming process is stuck, where the code is difficult to write it out.
main.py coding program was done

2025/12/4
Finished pendulum_model.py

2025/12/8
Finished chaos_tools.py, Writen some of the code for plotter.py

2025/12/10
Finished up plotter.py (angles-time plot, phase-space plot, sepration-time plot)
Run the program to DEBUG:
1. Fixed Runtime Error in 
main.py
- Issue: The function 
get_initial_state
 was assigned to a variable instead of being called, causing a ZeroDivisionError later because 
initial_state
 was a function object, not a list of values.
- Fix: Added parentheses to call the function: initial_state = get_initial_state().

2. Corrected Physics in 
pendulum_model.py
- Issue: The equations of motion for the second pendulum (domega2_dt) contained copy-paste errors:
den2 was calculated using l1 instead of l2.
num2 used l1 instead of l2 in the velocity-dependent term and had incorrect signs.
- Fix: Updated den2 and num2 to match standard Lagrangian mechanics for a double pendulum.

3. Verified 
plotter.py
- Action: Temporarily modified 
plotter.py
 to save plots to files for automated verification.
- Result: Confirmed that plotting functions work correctly.
- Cleanup: Reverted all changes to 
plotter.py
 to restore the original interactive plt.show() behavior.

Verification Results
- Execution: 
main.py
 runs without errors.
- Physics: A test simulation run confirmed that the corrected equations produce valid numerical output without crashing.
