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
    Because the coupled differential equations are non-linear and cannot be solved analytically, the program uses the **4th-Order Runge-Kutta (RK4)** method to numerically approximate the solution. This method estimates the slope at four points within a time step $dt$ to provide a highly accurate update:
    $$y_{n+1} = y_n + \frac{dt}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
    This logic is encapsulated in `pendulum_model.py`.

**Key Functions:**
1.  **Real-Time Simulation:** Solves the equations of motion on-the-fly as users adjust parameters.
2.  **Chaos Visualization:** Simulates two pendulums with a microscopic difference ($0.001$ rad) to visualize divergence.
3.  **Phase Space Analysis:** Plots $\omega$ vs $\theta$ to analyze the system's dynamic stability (limit cycles vs. chaos).

## 2. Usage
**Requirements:**
* Python 3.8+
* Libraries: `matplotlib`, `numpy`

**Execution:**
1.  Install dependencies: `pip install matplotlib numpy`
2.  Run the main dashboard: `python gui_dashboard.py`

**Interaction:**
* **Control Panel:** Use the "Smart Sliders" on the left.
    * **Click Arrows:** For precise increments (e.g., $\pm 0.001$ rad for angles).
    * **Drag Slider:** For coarse adjustments.
* **Simulation:** The plots update automatically. Use **"Replay"** to watch the animation again or **"Reset"** to restore default values.

## 3. Program Structure
The project follows a **Model-View-Controller (MVC)** pattern to separate physics from the interface.

* **`pendulum_model.py` (The Model):**
    * Contains the physics constants and differential equations (`derivatives` function).
    * Implements the `rk4_step` and `simulate` functions for numerical integration.
    * *Relation to Principles:* This file strictly handles the math described in Section 1.

* **`main.py` (The View & Controller):**
    * **SimulationManager Class:** Manages the state and calls the physics engine.
    * **SmartSlider Class:** A custom UI widget combining Matplotlib sliders with buttons.
    * **Layout Logic:** Handles the "Dashboard" arrangement (Plots, Graphs, and Controls).

## 4. References
* **Physics Equations:** Derived from standard Classical Mechanics textbooks (Lagrangian dynamics of a double pendulum).
* **AI Assistance (ChatGPT):**
    * Used for generating the initial skeleton code for the CLI version.
    * Consulted for debugging the `ZeroDivisionError` in the initial `main.py`.
    * Assisted in writing the `matplotlib.widgets` code for the GUI transition and designing the "Smart Slider" class.

## 5. Modifications & Enhancements
**From CLI to GUI Dashboard:**
The initial program was a simple script that ran one simulation and showed one plot. I have significantly enhanced this to be a scientific exploration tool.

**Specific Contributions:**
1.  **Smart Sliders (Custom Class):**
    * *Problem:* Standard sliders are hard to set to exact values (e.g., exactly 0.000 velocity).
    * *My Contribution:* I wrote a custom `SmartSlider` class that adds `<` and `>` buttons to every slider. This allows for precise, discrete increments (e.g., steps of 0.001 for angles), which is critical for testing sensitivity to initial conditions.
2.  **Dashboard Architecture:**
    * *Enhancement:* Instead of popup windows, I designed a single-window dashboard containing 4 simultaneous panels: Control, Trajectory, Chaos Comparison, and Data Analysis.
3.  **Real-Time Metrics:**
    * *Enhancement:* Added a "Separation vs. Time" graph that calculates the Euclidean distance between the two chaotic trajectories in real-time, providing a quantitative measure of the "Butterfly Effect" that wasn't present in the original code.
4.  **Merged Phase Space:**
    * *Enhancement:* Combined the phase space of both pendulums into one graph for easier comparison of their dynamic states.