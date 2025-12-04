import math

def derivatives(state, params):
    """
    Compute time derivatives for double pendulum
    """
    theta1, omega1, theta2, omega2 = state

    m1 = params["m1"]
    m2 = params["m2"]
    l1 = params["l1"]
    l2 = params["l2"]
    g = params["g"]

    #Angle difference
    delta = theta2 - theta1
    sin_delta = math.sin(delta)
    cos_delta = math.cos(delta)

    #Denominators from standard double pendulum equation
    den1 = (m1+m2)*l1 - m2*l1*(cos_delta**2)
    den2 = (l1/l1)*den1

    #Angular accelerations (aplha1, alpha2)
    #alpha1 = domega1_dt
    num1 = (m2 * l1 * (omega1**2) * sin_delta * cos_delta
            + m2 * g * math.sin(theta2) * cos_delta
            + m2 * l2 * (omega2**2) * sin_delta
            - (m1 + m2) * g * math.sin(theta1))
    domega1_dt = num1 / den1

    #aplha2 = domega2_dt
    num2 = (m2 * l1 * (omega2**2) * sin_delta * cos_delta
            + (m1 + m2) * (g *math.sin(theta1) * cos_delta 
            - l1 * (omega1**2) * sin_delta - g * math.sin(theta2)))
    domega2_dt = num2 / den2

    dtheta1_dt = omega1
    dtheta2_dt = omega2

    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt]

def rk4_step(state, dt, params):
    k1 = derivatives(state, params)

    state_k2 = [state[i] + 0.5 * dt *k1[i] for i in range(4)]
    k2 = derivatives(state_k2, params)

    state_k3 = [state[i] + 0.5 * dt *k2[i] for i in range(4)]
    k3 = derivatives(state_k3, params)

    state_k4 = [state[i] + 0.5 * dt *k3[i] for i in range(4)]
    k4 = derivatives(state_k4, params)
    
    new_state = []
    for i in range(4):
        increment = (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        new_state.append(state[i] + increment)
    return new_state

def simulate(initial_state, t_max, dt, params):
    n_steps = int(t_max / dt) + 1

    t_array = []
    theta1_array = []
    omega1_array = []
    theta2_array = []
    omega2_array = []

    #Initial conditions
    t = 0.0
    state = list(initial_state) # make copy

    for _ in range (n_steps):
        t_array.append(t)
        theta1_array.append(state[0])
        omega1_array.append(state[1])
        theta2_array.append(state[2])
        omega2_array.append(state[3])

        state = rk4_step(state, dt, params)
        t+=dt
    return t_array, theta1_array, omega1_array, theta2_array, omega2_array