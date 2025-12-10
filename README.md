# Chaotic Double Pendulum Simulator

## 1. Project Overview

**Short Description**
A Python-based simulator that solves and visualizes the **chaotic motion of a double pendulum**, demonstrating extreme sensitivity to initial conditions.

**Motivation**
The double pendulum is a simple mechanical system that can show chaos: two pendulum with almost identical starting angles quickly diverge into very different trajectories.
This project uses **numerical ODE methods** to solve the equations of motion for a double pendulum and **visualization** to explore chaos in classical mechanics.

## 2. Features
### Core Features
- Simulate the motion of a double pendulum with configurable:
    -masses \ (m1, m2\)
    -lengths \ (l1, l2\)
    -gravity \ (g\)
    -initial conditions \ (theta1_0, omega1_0, theta2_0, omega2_0\)
    -simulation time \ (t_max\)
    -time step \ (dt\)
- Numerically solve the equations of motion  using a **4th order Runge-Kutta method**
- Animate the double pendulum motion using 'matplotlib'
Plot angles vs time for both pendulum arms.

### Chaos & Analysis Features
- **Chaos Demo:** 
    - Run two simulations with a tiny difference in initial angle (e.g., theta_1 + 0.001 rad)
    - Compute and plot the separation between the two trajectories as a function of time
- **Phase-space plot:**
    - Plot the phase-space graph for one pendulum arm
    - (Optional) Estimate a simple **divergence rate** to illustrate Lyapunov-like behavior

## 3. Project Structure
project_root/
├── main.py              # CLI menu and user interaction
├── pendulum_model.py    # Equations of motion + RK4 integrator + simulation
├── chaos_tools.py       # Two-trajectory simulation + separation analysis
├── plotter.py           # Plotting and animation routines
├── README.md            # Describe Function, Usage, Design, Development process, Sources, Updates
└── dev_log.md           # Development log and notes (including AI usage)

## 4. Physics & Mathematical Background
**4.1 System Description**
The double pendulum constis of :
- First pendulum: mass m1 at the end of a massless rod of length l1
- Second pendulum: mass m2 attached to the end of the first mass by a massless rod of length l2
- Angle θ1​(t) and θ2​(t) are mesured from the vertical
- The motion is driven by gravity g, and constrained by the rods
The state vector is :

y(t) = [θ1​, ω1​, θ2​, ω2​]

where wi = dθi/dt is angular velocity.

**4.2 Equations of Motion (Conceptual)**
The equations of motion come from Lagrangian mechanics (kinetic and potential energy of the two masses).
The form a pair of nonlinear, coupled second-order ODEs in θ1 and θ2, often written in the form:
d2θ1/dt2 = f(θ1, θ2, ω1, ω2)
d2θ2/dt2 = g(θ1, θ2, ω1, ω2)

In the code, these are rewritten as a first-order system:
d/dt [θ1, ω1, θ2, ω2] = [ω1, f(θ1, θ2, ω1, ω2), ω2, g(θ1, θ2, ω1, ω2)]

**4.3 Numerical Methods: Runge-Kutta 4th order (RK4)**
The program uses a fixed time-step RK4 integrator to solve the equations of motion.
Given current state yn at time tn:
1. Compute k1 = f(tn, yn)
2. Compute k2 = f(tn + dt/2, yn + dt/2 * k1)
3. Compute k3 = f(tn + dt/2, yn + dt/2 * k2)
4. Compute k4 = f(tn + dt, yn + dt * k3)
Update:
yn+1 = yn + dt/6 * (k1 + 2k2 + 2k3 + k4)
This method provides good accurancy for sufficient small time steps dt.


## 5. Ways to Run
**5.1 Requirements**
- Python 3.8 or higher
- Recommended packeages: matplotlib, numpy

**5.2 Setup and Execution**
- clone or download the project
git clone https://github.com/eugtan98/Chaotic-Double-Pendulum-Simulator.git

- install dependencies
pip install matplotlib numpy

- run the program
python main.py





