from typing import Dict, List, Tuple
import math
import pendulum_model   

def simulate_two_trajectories(intiial_state: List[float], 
                              delta_theta1: float, 
                              t_max: float, 
                              dt: float, 
                              params: Dict[str, float],
                              )-> Tuple[List[float],  # t_array
                                        List[float], List[float], List[float], List[float],  # traj A: theta1, omega1, theta2, omega2
                                        List[float], List[float], List[float], List[float],  # traj B: theta1, omega1, theta2, omega2
                                        ]:
    state_a = list(intiial_state)
    state_b = list(intiial_state)

    state_b[0] += delta_theta1

    (t_array, theta1_a, omega1_a, 
     theta2_a, omega2_a,) = pendulum_model.simulate(state_a, t_max, dt, params)
    
    (t_array_b, theta1_b, omega1_b, 
     theta2_b, omega2_b,) = pendulum_model.simulate(state_b, t_max, dt, params)
    
    n = min(len(t_array), len(theta1_a), len(theta1_b), len(theta2_a), len(theta2_b),)

    t_array = t_array[:n]
    theta1_a = theta1_a[:n]
    omega1_a = omega1_a[:n]
    theta2_a = theta2_a[:n]
    omega2_a = omega2_a[:n]

    theta1_b = theta1_b[:n]
    omega1_b = omega1_b[:n]
    theta2_b = theta2_b[:n]
    omega2_b = omega2_b[:n]

    return (t_array,
            theta1_a,
            omega1_a,
            theta2_a,
            omega2_a,
            theta1_b,
            omega1_b,
            theta2_b,
            omega2_b,)

def compute_separation(theta1_a : List[float], 
                       theta2_a: List[float], 
                       theta1_b: List[float], 
                       theta2_b: List[float],
                       )-> List[float]:
    n = min(len(theta1_a), len(theta1_b), len(theta2_a), len(theta2_b))
    separation: List[float] = []

    for i in range(n):
        d_theta1 = theta1_a[i] - theta1_b[i]
        d_theta2 = theta2_a[i] - theta2_b[i]

        d = math.sqrt(d_theta1 * d_theta1 + d_theta2 * d_theta2)
        separation.append(d)

    return separation