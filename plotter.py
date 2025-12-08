import math
from typing import List, Dict

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


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
        return rod_line
    
    ani = FuncAnimation(fig, update, frames=len(theta1_array), init_func=init, interval=20, blit = True)
    plt.show