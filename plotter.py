import math
from typing import List, Dict

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def compute_positions(theta1_array: List[float], 
                      theta2_array: List[float],
                      l1: float,
                      l2: float,
):
    """
    Helper function to compute (x1, y1), (x2, y2) positions from angles
    """
    x1_array, y1_array = [], []
    x2_array, y2_array = [], []

    for theta1, theta2 in zip(theta1_array, theta2_array):
        x1 = l1 * math.sin(theta1)
        y1 = -l1 * math.cos(theta1)
        x2 = x1 + l2 * math.sin(theta2)
        y2 = y1 - l2 * math.cos(theta2)

        x1_array.append(x1)
        y1_array.append(y1)
        x2_array.append(x2)
        y2_array.append(y2)

    return (x1_array, y1_array, x2_array, y2_array)

def animate_pendulum(theta1_array, theta2_array, params):
    """
    Animate the motion of a single double pendulum
    """
    l1 = params["l1"]
    l2 = params["l2"]

    x1_array = []
    y1_array = []
    x2_array = []
    y2_array = []

    for i in range(len(theta1_array)):
        theta1 = theta1_array[i]
        theta2 = theta2_array[i]

        x1 = l1 * math.sin(theta1)
        y1 = -l1 * math.cos(theta1)
        x2 = x1 + l2 * math.sin(theta2)
        y2 = y1 - l2 * math.cos(theta2)

        x1_array.append(x1)
        y1_array.append(y1)
        x2_array.append(x2)
        y2_array.append(y2)

    fig, ax = plt.subplots()

    ax.set_aspect('equal', adjustable='box')

    limit = l1 +l2
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Double Pendulum Animation")

    rod_line, = ax.plot([], [], "-o")

    def init():
        rod_line.set_data([], [])
        return rod_line,
    
    def update(frame_index):
        x1 = x1_array[frame_index]
        y1 = y1_array[frame_index]
        x2 = x2_array[frame_index]
        y2 = y2_array[frame_index]

        xs = [0.0, x1, x2]
        ys = [0.0, y1, y2]
        rod_line.set_data(xs, ys)
        return rod_line,
    
    ani = FuncAnimation(fig, update, frames=len(theta1_array), init_func=init, interval=20, blit = True)
    plt.show()

def animate_two_pendulums(theta1_a: List[float], 
                          theta2_a: List[float], 
                          theta1_b: List[float], 
                          theta2_b: List[float], 
                          params: Dict[str, float], 
                          interval: int=20,
                          ):
    """
    Animate the motion of two double pendulums side by side
    """
    l1 = params["l1"]
    l2 = params["l2"]

    # Ensure same length (use the shorter)
    n = min(len(theta1_a), len(theta2_a), len(theta1_b), len(theta2_b))
    theta1_a = theta1_a[:n]
    theta2_a = theta2_a[:n]
    theta1_b = theta1_b[:n]
    theta2_b = theta2_b[:n]

    #Precompute position
    (x1_a, y1_a, x2_a, y2_a) = compute_positions(theta1_a, theta2_a, l1, l2)
    (x1_b, y1_b, x2_b, y2_b) = compute_positions(theta1_b, theta2_b, l1, l2)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    for ax in (ax1, ax2):
        ax.set_aspect('equal', adjustable='box')
        limit = l1 + l2
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    ax1.set_title("Pendulum A")
    ax2.set_title("Pendulum B")

    line_a, = ax1.plot([], [], "-o")
    line_b, = ax2.plot([], [], "-o")

    def init():
        line_a.set_data([], [])
        line_b.set_data([], [])
        return line_a, line_b
    
    def update(frame_index):
        # Pendulum A
        x1a = x1_a[frame_index]
        y1a = y1_a[frame_index]
        x2a = x2_a[frame_index]
        y2a = y2_a[frame_index]

        # Pendulum B
        x1b = x1_b[frame_index]
        y1b = y1_b[frame_index]
        x2b = x2_b[frame_index]
        y2b = y2_b[frame_index]

        line_a.set_data([0.0, x1a, x2a], [0.0, y1a, y2a])
        line_b.set_data([0.0, x1b, x2b], [0.0, y1b, y2b])
        return line_a, line_b
    ani = FuncAnimation(
        fig,
        update,
        frames=n,
        init_func = init,
        interval=interval,
        blit=True,  
    )

    plt.show()

def plot_angles_vs_time(t_array: List[float], 
                        theta1_array: List[float], 
                        theta2_array: List[float],
                        ):
    """
    Plot angles vs time
    """
    fig, ax = plt.subplots()
    ax.plot(t_array, theta1_array, label="Theta 1")
    ax.plot(t_array, theta2_array, label="Theta 2")
    ax.set_xlabel("Time (s)")      
    ax.set_ylabel("Angle (rad)")
    ax.set_title("Double Pendulum Angles vs Time")
    ax.legend()
    ax.grid(True)
    ax.grid(True)
    plt.show()

def plot_phase_space(theta_array: List[float], 
                     omega_array: List[float],
                     ):
    """
    Plot phase-space (theta vs omega)
    """
    fig, ax = plt.subplots()
    ax.plot(theta_array, omega_array)
    ax.set_xlabel("Angle (rad)")
    ax.set_ylabel("Angular Velocity (rad/s)")
    ax.set_title("Phase-Space Plot")
    ax.grid(True)
    ax.grid(True)
    plt.show()

def plot_separation_vs_time(t_array: List[float],
                            separation_array: List[float],
                            ):
    """
    Plot separation vs time
    """
    fig, ax = plt.subplots()
    ax.plot(t_array, separation_array)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Separation in Angle Space")
    ax.set_title("Separation vs Time")
    ax.grid(True)
    ax.grid(True)
    plt.show()
