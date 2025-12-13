# Chaotic Double Pendulum Simulator & Dashboard

## 1. Principles and Functions
**Scientific Background:**
The double pendulum is a classic example of a physical system that exhibits **deterministic chaos**. While governed by known laws of physics (Newtonian/Lagrangian mechanics), its motion is highly sensitive to initial conditions.

**Physics & Mathematics:**
* **Lagrangian Mechanics:** The equations of motion are derived using the Lagrangian $L = T - V$.
* **Numerical Integration (RK4):** The program uses the **4th-Order Runge-Kutta** method to numerically approximate the solution, minimizing error drift over time.
* **Rotational Inertia ($I$):** The simulator calculates the Moment of Inertia in real-time. It tracks the inertia of individual masses and the dynamic **System Inertia** ($I_{sys}$) relative to the pivot, which fluctuates as the pendulum changes shape.
* **Energy Conservation:** Ideally, Total Energy $E = K + U$ remains constant. The dashboard tracks this to verify simulation accuracy.

**Key Functions:**
1.  **3D Interactive Simulation:** Visualizes the pendulum in a 3D space that users can rotate 360° to view motion from any angle.
2.  **Chaos Visualization:** Simulates two double pendulums with a microscopic difference ($0.001$ rad) to visualize the "Butterfly Effect."
3.  **Real-Time Physics Metrics:** Displays live calculations for Energy, Inertia, and Lyapunov Exponents.

## 2. Usage
**Requirements:**
* Python 3.8+
* Libraries: `matplotlib`, `numpy`

**Execution:**
1.  Install dependencies: `pip install matplotlib numpy`
2.  Run the dashboard: `python main.py`

**Interaction:**
* **Control Panel (Left):** Use "Smart Sliders" to adjust Mass, Length, Gravity, and Initial Conditions with 0.01 precision.
* **3D Plots (Right):** Click and drag the "Single Trajectory" or "Chaos Demo" plots to rotate the view.
* **Playback:** Use **Pause**, **Replay**, and **Reset** to control the simulation flow for teaching or analysis.

## 3. Program Structure
The project follows a **Model-View-Controller (MVC)** pattern.

* **`pendulum_model.py` (The Model):**
    * Contains the physics constants and differential equations.
    * Implements the RK4 integration logic.
* **`main.py` (The View & Controller):**
    * **SimulationManager Class:** Manages state, physics calculations (Energy, Inertia), and threading.
    * **Layout Logic:** Handles the split-screen design (Left Controls vs. Right 2x3 Plot Grid).
    * **SmartSlider Class:** Custom widget combining sliders with precise increment buttons.

## 4. Development Process
**Process & Challenges:**
1.  **Transition from CLI to GUI:** Originally a command-line tool which is a limited usage of program that needed user to input all data manually, the project was then refactored into a real-time interactive dashboard to improve the interactions with user.
2.  **Physics Verification:** Initially, the simulation showed energy drift. This was addressed by refining the time-step logic (`t = i * dt`) and adding an Energy vs. Time graph for visual verification.
3.  **Layout Optimization:** Arranging 6 plots and detailed controls without overlapping was difficult. I solved this by implementing a strict grid layout and moving all controls to a dedicated left-hand column.

## 5. References
* **Physics Equations:** Standard Lagrangian dynamics of a double pendulum.
* **AI Assistance:** Used for debugging `matplotlib` layout issues and generating the boilerplate for the `Axes3D` integration. (Google Gemini, Google Antigravity)

## 6. Modifications & Enhancements
**Self Contributions:**
1.  **3D Visualization:**
    * Upgraded the standard 2D plots to **interactive 3D axes**, allowing users to rotate the view and understand the spatial context of the motion.
2.  **Advanced Physics Metrics:**
    * **Inertia:** Added real-time calculation of $I_1$, $I_2$, and the complex $I_{sys}$ (System Inertia).
    * **Lyapunov & Separation:** Quantified chaos by plotting the divergence distance and estimating the Lyapunov exponent.
3.  **Smart UI Design:**
    * Designed a custom **"Smart Slider"** that includes `<` and `>` buttons for precise 0.01 adjustments, solving the difficulty of setting exact initial conditions on standard sliders.
    * Created a **2x3 Grid Dashboard** that displays Trajectory, Angles, Energy, Phase Space, and Separation simultaneously.