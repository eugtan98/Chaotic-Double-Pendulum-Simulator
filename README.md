# Chaotic Double Pendulum Simulator & Dashboard

## 1. Principles and Functions
**Scientific Background:**
The double pendulum is a classic example of a physical system that exhibits **deterministic chaos**. While governed by known laws of physics (Newtonian/Lagrangian mechanics), its motion is highly sensitive to initial conditions.

**Physics & Mathematics:**
* **Lagrangian Mechanics:** The equations of motion are derived using the Lagrangian $L = T - V$, where $T$ is kinetic energy and $V$ is potential energy.
* **System State:** The system is described by four variables:
    $$y(t) = [\theta_1, \omega_1, \theta_2, \omega_2]$$
    Where $\theta$ is the angle and $\omega$ is the angular velocity.
* **Numerical Integration (RK4):**
    The program uses the **4th-Order Runge-Kutta (RK4)** method to numerically approximate the solution. This method estimates the slope at four points within a time step $dt$ to provide a highly accurate update:
    $$y_{n+1} = y_n + \frac{dt}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
    This logic is encapsulated in `pendulum_model.py`.
* **Lyapunov Exponent ($\lambda$):**
    To quantify chaos, the dashboard estimates the Lyapunov exponent in real-time using the divergence rate between two trajectories:
    $$\lambda \approx \frac{1}{t} \ln\left(\frac{|\delta Z(t)|}{|\delta Z_0|}\right)$$
    A positive $\lambda$ indicates chaotic behavior.
* **Energy Conservation ($E$):**
    In an ideal frictionless system, the total mechanical energy $E$ must remain constant. The dashboard calculates this in real-time to verify the accuracy of the numerical simulation:
    $$E(t) = K(t) + U(t)$$

**Key Functions:**
1.  **Real-Time Simulation:** Solves the equations of motion on-the-fly as users adjust parameters.
2.  **Chaos Visualization:** Simulates two pendulums with a microscopic difference ($0.001$ rad) to visualize divergence.
3.  **Phase Space Analysis:** Plots $\omega$ vs $\theta$ to analyze the system's dynamic stability.
4.  **Energy Analysis:** Plots Kinetic, Potential, and Total Energy over time to visualize energy transfer and numerical drift.

## 2. Usage
**Requirements:**
* Python 3.8+
* Libraries: `matplotlib`, `numpy`

**Execution:**
1.  Install dependencies: `pip install matplotlib numpy`
2.  Run the main dashboard: `python gui_dashboard.py`

**Interaction:**
* **Control Panel:** Use the "Smart Sliders" on the left to adjust Mass, Length, Gravity, and Initial Conditions.
* **Playback:**
    * **Pause/Resume:** Freeze the simulation to analyze specific moments.
    * **Replay:** Restart the animation from $t=0$ with current settings.
    * **Reset:** Restore all parameters to default values.
* **Zooming:** You can use the Matplotlib zoom/lasso tool on any graph. Pressing "Reset" or "Replay" will automatically reset the view to fit the data.

## 3. Program Structure
The project follows a **Model-View-Controller (MVC)** pattern to separate physics from the interface.

* **`pendulum_model.py` (The Model):**
    * Contains the physics constants and differential equations (`derivatives` function).
    * Implements the `rk4_step` and `simulate` functions for numerical integration.
    * *Relation to Principles:* This file strictly handles the math described in Section 1.

* **`main.py` (The View & Controller):**
    * **SimulationManager Class:** Manages the state, runs the physics engine, and calculates separation, energy, and Lyapunov metrics.
    * **SmartSlider Class:** A custom UI widget combining Matplotlib sliders with increment buttons.
    * **Layout Logic:** Handles the 6-panel dashboard arrangement (Single Sim, Chaos Demo, Angle Plot, Energy Plot, Phase Space, Separation Plot).

## 4. Development Process
**Process & Challenges:**
1.  **Initial CLI Phase:** The project started as a command-line tool (`main.py`) where users had to manually input numbers.
2.  **Debugging Physics:**
    * *Issue:* The second pendulum's motion was physically incorrect.
    * *Solution:* Found a copy-paste error in `pendulum_model.py` where `l1` was used instead of `l2` in the `num2` equation. Corrected the signs to match standard Lagrangian mechanics.
3.  **Debugging Math:**
    * *Issue:* Simulation inaccuracies.
    * *Solution:* Identified an error in the RK4 implementation where `0.5 * dt` was used incorrectly in the final step. Corrected to standard RK4 form.
4.  **UI Pivot:** Realized that text input was poor for exploring chaos. Decided to build a GUI Dashboard inspired by ideal gas simulators shown by professor.
5.  **Layout Challenges:** The labels in Matplotlib overlapped with sliders. I resolved this by manually calculating coordinate offsets and creating a custom Layout Manager.

## 5. References
* **Physics Equations:** Derived from standard Classical Mechanics textbooks (Lagrangian dynamics of a double pendulum).
* **AI Assistance (ChatGPT、Gemini、Antigravity):**
    * Used for generating the initial pseudocode and skeleton code for the CLI version.(GPT)
    * Consulted for debugging the `ZeroDivisionError` in the initial `main.py`.(Antigravity)
    * Assisted in writing the `matplotlib.widgets` code for the GUI transition and designing the "Smart Slider" class.(Gemini)

## 6. Modifications & Enhancements
**From CLI to GUI Dashboard:**
The original reference was a simple script that ran one simulation. I have significantly enhanced this into a scientific exploration tool.

**Self Contributions:**
1.  **Smart Sliders (Custom Class):**
    * I wrote a custom `SmartSlider` class that adds `<` and `>` buttons to every slider. This allows for precise, discrete increments (e.g., steps of 0.001 for angles), which is critical for testing sensitivity to initial conditions.
2.  **Advanced Scientific Metrics:**
    * **Separation vs. Time:** Quantifies the "Butterfly Effect" by plotting the distance between chaotic trajectories.
    * **Lyapunov Exponent:** Implemented real-time estimation of the Lyapunov exponent to scientifically confirm chaotic behavior.
    * **Energy Verification:** Added an **Energy vs. Time** graph plotting $K$, $U$, and Total Energy ($E$). This allows users to verify the conservation of energy and check the accuracy of the numerical method.
3.  **Educational Tools:**
    * Added a **Pause/Resume** feature to allow user to freeze the simulation for explanation.
    * Added an **Angle vs. Time** graph to visualize periodicity.
    * **Smart View Reset:** Implemented logic to automatically reset zoomed graph views only when necessary (Replay/Reset), improving the user experience.
4.  **Merged Phase Space:**
    * Combined the phase space of both pendulums into one graph for easier comparison of their dynamic states.