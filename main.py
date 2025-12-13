import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.patches as patches # No longer needed
import pendulum_model

# --- 1. Simulation Manager ---
class SimulationManager:
    def __init__(self):
        # Default Parameters
        self.params = {
            "m1": 1.0, "m2": 1.0, 
            "l1": 1.0, "l2": 1.0, 
            "g": 9.81
        }
        # Initial State [theta1, omega1, theta2, omega2]
        self.initial_state = [np.pi/2, 0.0, np.pi/2, 0.0] 
        self.t_max = 10.0 
        self.dt = 0.05
        
        # Playback control
        self.current_frame = 0 
        self.is_running = True 
        self.is_paused = False
        
        self.sol_a = {} 
        self.sol_b = {} 
        self.separation = [] 
        self.perturbation = 0.001 
        
        self.run_sim()

    def calculate_energy(self, params, theta1, omega1, theta2, omega2):
        """Calculates Kinetic (K), Potential (U), and Total (E) Energy."""
        m1, m2 = params['m1'], params['m2']
        l1, l2 = params['l1'], params['l2']
        g = params['g']
        
        t1, w1 = np.array(theta1), np.array(omega1)
        t2, w2 = np.array(theta2), np.array(omega2)
        
        # Potential Energy (Reference y=0 at pivot)
        y1 = -l1 * np.cos(t1)
        y2 = y1 - l2 * np.cos(t2)
        U = m1 * g * y1 + m2 * g * y2
        
        # Kinetic Energy
        v1_sq = (l1 * w1)**2
        v2_sq = (l1 * w1)**2 + (l2 * w2)**2 + 2 * l1 * l2 * w1 * w2 * np.cos(t1 - t2)
        K = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq
        
        return K, U, K + U

    def calculate_trajectory(self, state, params):
        (t_arr, theta1, omega1, theta2, omega2) = pendulum_model.simulate(
            state, self.t_max, self.dt, params
        )
        l1 = params['l1']
        l2 = params['l2']
        t1, t2 = np.array(theta1), np.array(theta2)
        
        x1 = l1 * np.sin(t1)
        y1 = -l1 * np.cos(t1)
        x2 = x1 + l2 * np.sin(t2)
        y2 = y1 - l2 * np.cos(t2)
        
        K, U, E = self.calculate_energy(params, theta1, omega1, theta2, omega2)
        
        return {
            "t": t_arr, "x1": x1, "y1": y1, "x2": x2, "y2": y2,
            "theta1": theta1, "omega1": omega1,
            "theta2": theta2, "omega2": omega2,
            "K": K, "U": U, "E": E
        }

    def run_sim(self):
        # Sim A (Standard)
        self.sol_a = self.calculate_trajectory(self.initial_state, self.params)
        
        # Sim B (Chaos)
        state_b = list(self.initial_state)
        state_b[0] += self.perturbation
        self.sol_b = self.calculate_trajectory(state_b, self.params)
        
        # Separation
        x2a = self.sol_a['x2']
        y2a = self.sol_a['y2']
        x2b = self.sol_b['x2']
        y2b = self.sol_b['y2']
        self.separation = np.sqrt((x2a - x2b)**2 + (y2a - y2b)**2)
        
        # Reset playhead
        self.current_frame = 0
        self.is_running = True
        self.is_paused = False

sim = SimulationManager()

# --- 2. Dashboard Layout ---
fig = plt.figure(figsize=(18, 10))

# === LEFT COLUMN: Controls & Results (Shifted Left) ===
fig.text(0.02, 0.96, "CONTROL PANEL", fontsize=16, fontweight='bold')
fig.text(0.02, 0.93, "Parameters:", fontsize=11, fontweight='bold')

# === RIGHT SIDE: PLOTS (2x3 Grid) ===
# Row 1: 3D Visuals
ax_main = plt.axes([0.36, 0.68, 0.28, 0.28], projection='3d') 
ax_main.set_title("Single Trajectory")
ax_main.set_xlabel("x (m)") 
ax_main.set_ylabel("z (m)") 
ax_main.set_zlabel("y (m)") 

ax_chaos = plt.axes([0.69, 0.68, 0.28, 0.28], projection='3d')
ax_chaos.set_title(f"Chaos Demo (Diff: {sim.perturbation} rad)")
ax_chaos.set_xlabel("x (m)") 
ax_chaos.set_ylabel("z (m)") 
ax_chaos.set_zlabel("y (m)") 

# Row 2: Time Series
ax_angle = plt.axes([0.36, 0.38, 0.28, 0.20])
ax_angle.set_title("Angles - Time")
ax_angle.set_xlabel("Time (s)")
ax_angle.set_ylabel("Angle (rad)")
ax_angle.grid(True)

ax_energy = plt.axes([0.69, 0.38, 0.28, 0.20])
ax_energy.set_title("Energy - Time")
ax_energy.set_xlabel("Time (s)")
ax_energy.set_ylabel("Energy (J)")
ax_energy.grid(True)

# Row 3: Analysis
ax_phase = plt.axes([0.36, 0.08, 0.28, 0.20])
ax_phase.set_title("Phase Space")
ax_phase.set_xlabel("Angle (rad)")
ax_phase.set_ylabel("Velocity (rad/s)")
ax_phase.grid(True)

ax_sep = plt.axes([0.69, 0.08, 0.28, 0.20])
ax_sep.set_title("Separation - Time")
ax_sep.set_xlabel("Time (s)")
ax_sep.set_ylabel("Distance (m)")
ax_sep.grid(True)

# === VISUAL ELEMENTS ===
line_main, = ax_main.plot([0], [0], [0], 'o-', lw=3, color='black')
line_chaos_a, = ax_chaos.plot([0], [0], [0], 'o-', lw=2, color='blue', alpha=0.6, label='A')
line_chaos_b, = ax_chaos.plot([0], [0], [0], 'o-', lw=2, color='red', alpha=0.6, label='B')
ax_chaos.legend(loc='upper right', fontsize='x-small')
line_angle1, = ax_angle.plot([], [], '-', lw=2, color='green', label='Theta 1')
line_angle2, = ax_angle.plot([], [], '-', lw=2, color='purple', label='Theta 2')
ax_angle.legend(loc='upper right', fontsize='x-small')
line_E, = ax_energy.plot([], [], '-', lw=2, color='black', label='Total (E)')
line_K, = ax_energy.plot([], [], '--', lw=1, color='blue', label='Kinetic (K)')
line_U, = ax_energy.plot([], [], '--', lw=1, color='red', label='Potential (U)')
ax_energy.legend(loc='upper right', fontsize='x-small')
point_phase1, = ax_phase.plot([], [], 'o', color='green', label='P1')
point_phase2, = ax_phase.plot([], [], 'o', color='purple', label='P2')
ax_phase.legend(loc='upper right', fontsize='x-small')
line_sep, = ax_sep.plot([], [], '-', lw=2, color='orange')


# --- 3. Control Panel Implementation (Standard Slider) ---

class SmartSlider:
    def __init__(self, fig, label, y_pos, vmin, vmax, valinit, step_val):
        self.step_val = step_val
        self.valinit = valinit
        
        # Decrease Button (<)
        ax_dec = plt.axes([0.07, y_pos, 0.015, 0.025]) 
        self.b_dec = Button(ax_dec, '<', color='lightgray', hovercolor='0.9')
        self.b_dec.on_clicked(self.decrement)
        
        # Slider (Middle) - Using default style
        ax_slide = plt.axes([0.09, y_pos, 0.10, 0.025])
        self.slider = Slider(ax_slide, "", vmin, vmax, valinit=valinit) 
        self.slider.valtext.set_visible(False)
        
        # No custom handle code here anymore.

        # Label (Left of button)
        fig.text(0.06, y_pos + 0.005, label, ha='right', va='center', fontsize=10)
        
        # Increase Button (>)
        ax_inc = plt.axes([0.20, y_pos, 0.015, 0.025]) 
        self.b_inc = Button(ax_inc, '>', color='lightgray', hovercolor='0.9')
        self.b_inc.on_clicked(self.increment)
        
        # Value Text (Right of button)
        self.val_text = fig.text(0.22, y_pos + 0.005, f"{valinit:.2f}", 
                                 ha='left', va='center', fontsize=10, color='blue', fontweight='bold')

        # Connect updates
        self.slider.on_changed(self.update_text)
        
    def update_text(self, val):
        self.val_text.set_text(f"{val:.2f}")

    def increment(self, event):
        new_val = self.slider.val + self.step_val
        if new_val <= self.slider.valmax:
            self.slider.set_val(new_val)
            
    def decrement(self, event):
        new_val = self.slider.val - self.step_val
        if new_val >= self.slider.valmin:
            self.slider.set_val(new_val)

    def on_changed(self, func):
        self.slider.on_changed(func)

    def reset(self):
        self.slider.reset()
        
    @property
    def val(self):
        return self.slider.val

# Create Sliders
y_start = 0.90
y_step = 0.035
s_time = SmartSlider(fig, 'Time (s)',      y_start, 5.0, 40.0, sim.t_max, 0.01)
s_g    = SmartSlider(fig, 'Gravity',       y_start - 1*y_step, 1.0, 20.0, sim.params['g'], 0.01)
s_m1   = SmartSlider(fig, 'Mass 1',        y_start - 2*y_step, 0.1, 5.0,  sim.params['m1'], 0.01)
s_m2   = SmartSlider(fig, 'Mass 2',        y_start - 3*y_step, 0.1, 5.0,  sim.params['m2'], 0.01)
s_l1   = SmartSlider(fig, 'Len 1',         y_start - 4*y_step, 0.5, 3.0,  sim.params['l1'], 0.01)
s_l2   = SmartSlider(fig, 'Len 2',         y_start - 5*y_step, 0.5, 3.0,  sim.params['l2'], 0.01)

# Initial Conditions Group
y_ic_start = y_start - 7*y_step
# Moved text up slightly to fix overlap
fig.text(0.02, y_ic_start + 0.035, "Initial Conditions:", fontsize=11, fontweight='bold')

s_th1  = SmartSlider(fig, 'Theta 1',       y_ic_start, -np.pi, np.pi, sim.initial_state[0], 0.01)
s_th2  = SmartSlider(fig, 'Theta 2',       y_ic_start - 1*y_step, -np.pi, np.pi, sim.initial_state[2], 0.01)
s_om1  = SmartSlider(fig, 'Omega 1',       y_ic_start - 2*y_step, -10.0, 10.0,   sim.initial_state[1], 0.01)
s_om2  = SmartSlider(fig, 'Omega 2',       y_ic_start - 3*y_step, -10.0, 10.0,   sim.initial_state[3], 0.01)

# Buttons (Shifted Left)
btn_y = y_ic_start - 5*y_step
ax_replay = plt.axes([0.02, btn_y, 0.09, 0.04])
b_replay = Button(ax_replay, 'Replay', hovercolor='0.9')
ax_reset = plt.axes([0.12, btn_y, 0.09, 0.04])
b_reset = Button(ax_reset, 'Reset', hovercolor='0.9', color='mistyrose')
ax_pause = plt.axes([0.22, btn_y, 0.09, 0.04])
b_pause = Button(ax_pause, 'Pause', hovercolor='0.9', color='lightcyan')


# === RESULTS PANEL (Shifted Left) ===
results_y_top = btn_y - 0.05
fig.text(0.02, results_y_top, "--- SIMULATION RESULTS ---", fontsize=12, fontweight='bold', ha='left')

# Main Stats Text
stats_text = fig.text(0.02, results_y_top - 0.03, '', 
                      va='top', ha='left', fontsize=10, 
                      fontfamily='monospace', 
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))


# --- 4. Update Logic & View Reset ---
smart_sliders = [s_time, s_g, s_m1, s_m2, s_l1, s_l2, s_th1, s_th2, s_om1, s_om2]

def reset_graph_views():
    """Forces all axes to reset limits based on current data."""
    # 1. Spatial Limits (3D)
    limit = sim.params['l1'] + sim.params['l2'] + 0.2
    for ax in [ax_main, ax_chaos]:
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit) 
        ax.set_zlim(-limit, limit) 
        ax.set_box_aspect((1, 1, 1)) 
    
    all_theta = np.concatenate([sim.sol_a['theta1'], sim.sol_a['theta2']])
    all_omega = np.concatenate([sim.sol_a['omega1'], sim.sol_a['omega2']])
    
    # 2. Angles vs Time
    ax_angle.set_xlim(0, sim.t_max)
    if len(all_theta) > 0:
        ax_angle.set_ylim(np.min(all_theta) - 1.0, np.max(all_theta) + 1.0)
    
    # 3. Energy vs Time
    ax_energy.set_xlim(0, sim.t_max)
    if len(sim.sol_a['E']) > 0:
        vals = np.concatenate([sim.sol_a['K'], sim.sol_a['U'], sim.sol_a['E']])
        ax_energy.set_ylim(np.min(vals) - 5, np.max(vals) + 5)

    # 4. Phase Space
    if len(all_theta) > 0:
        ax_phase.set_xlim(np.min(all_theta) - 0.5, np.max(all_theta) + 0.5)
        ax_phase.set_ylim(np.min(all_omega) - 0.5, np.max(all_omega) + 0.5)
    
    # 5. Separation
    ax_sep.set_xlim(0, sim.t_max)
    if len(sim.separation) > 0:
        ax_sep.set_ylim(0, np.max(sim.separation) * 1.1)

def update_sim(val=None):
    sim.t_max = s_time.val
    sim.params['m1'] = s_m1.val
    sim.params['m2'] = s_m2.val
    sim.params['l1'] = s_l1.val
    sim.params['l2'] = s_l2.val
    sim.params['g']  = s_g.val
    sim.initial_state = [s_th1.val, s_om1.val, s_th2.val, s_om2.val]
    
    sim.run_sim()
    reset_graph_views()

    sim.is_paused = False
    b_pause.label.set_text("Pause")

def replay(event):
    sim.current_frame = 0
    sim.is_running = True
    sim.is_paused = False
    b_pause.label.set_text("Pause")
    reset_graph_views() 

def reset_all(event):
    for s in smart_sliders:
        s.reset()
    replay(event)

def toggle_pause(event):
    sim.is_paused = not sim.is_paused
    if sim.is_paused:
        b_pause.label.set_text("Resume")
    else:
        b_pause.label.set_text("Pause")
        reset_graph_views()

for s in smart_sliders:
    s.on_changed(update_sim)

b_replay.on_clicked(replay)
b_reset.on_clicked(reset_all)
b_pause.on_clicked(toggle_pause)

# Initial run
update_sim()

# --- 5. Animation Loop ---
def animate(frame):
    if not sim.is_running:
        return line_main, line_chaos_a, line_chaos_b, line_angle1, line_angle2, line_E, line_K, line_U, point_phase1, point_phase2, line_sep, stats_text
    
    if sim.is_paused:
        return line_main, line_chaos_a, line_chaos_b, line_angle1, line_angle2, line_E, line_K, line_U, point_phase1, point_phase2, line_sep, stats_text

    i = sim.current_frame
    
    if i >= len(sim.sol_a['x1']) - 1:
        sim.is_running = False 
        i = len(sim.sol_a['x1']) - 1 
    
    # 1. Main Sim (3D UPDATE)
    x1, y1 = sim.sol_a['x1'][i], sim.sol_a['y1'][i]
    x2, y2 = sim.sol_a['x2'][i], sim.sol_a['y2'][i]
    
    line_main.set_data([0, x1, x2], [0, 0, 0]) 
    line_main.set_3d_properties([0, y1, y2])
    
    # 2. Chaos Sim (3D UPDATE)
    line_chaos_a.set_data([0, x1, x2], [0, 0, 0])
    line_chaos_a.set_3d_properties([0, y1, y2])
    
    xb1, yb1 = sim.sol_b['x1'][i], sim.sol_b['y1'][i]
    xb2, yb2 = sim.sol_b['x2'][i], sim.sol_b['y2'][i]
    line_chaos_b.set_data([0, xb1, xb2], [0, 0, 0])
    line_chaos_b.set_3d_properties([0, yb1, yb2])
    
    # 3. Angle vs Time
    line_angle1.set_data(sim.sol_a['t'][:i+1], sim.sol_a['theta1'][:i+1])
    line_angle2.set_data(sim.sol_a['t'][:i+1], sim.sol_a['theta2'][:i+1])
    
    # 4. Energy vs Time
    line_E.set_data(sim.sol_a['t'][:i+1], sim.sol_a['E'][:i+1])
    line_K.set_data(sim.sol_a['t'][:i+1], sim.sol_a['K'][:i+1])
    line_U.set_data(sim.sol_a['t'][:i+1], sim.sol_a['U'][:i+1])
    
    # 5. Phase Space
    point_phase1.set_data([sim.sol_a['theta1'][i]], [sim.sol_a['omega1'][i]])
    point_phase2.set_data([sim.sol_a['theta2'][i]], [sim.sol_a['omega2'][i]])
    
    # 6. Separation
    line_sep.set_data(sim.sol_a['t'][:i+1], sim.separation[:i+1])
    
    # Results & INERTIA CALCULATION
    diff_pos = sim.separation[i]
    t_curr = sim.sol_a['t'][i]
    current_E = sim.sol_a['E'][i]
    
    # Parameters for Inertia
    m1, m2 = sim.params['m1'], sim.params['m2']
    l1, l2 = sim.params['l1'], sim.params['l2']
    th1 = sim.sol_a['theta1'][i]
    th2 = sim.sol_a['theta2'][i]
    
    # 1. Inertia 1 (About Pivot)
    I1 = m1 * (l1**2)
    # 2. Inertia 2 (About Joint)
    I2 = m2 * (l2**2)
    # 3. Total System Inertia (About Pivot) - Dynamic
    dist_sq = l1**2 + l2**2 + 2*l1*l2*np.cos(th1 - th2)
    I_sys = m1 * (l1**2) + m2 * dist_sq
    
    # Lyapunov Estimate
    lyap_est = 0.0
    if t_curr > 0.1 and sim.separation[0] > 0:
        lyap_est = (1.0 / t_curr) * np.log(diff_pos / sim.separation[0])

    # Formatting Results Text
    results_str = (
        f"Time:      {t_curr:.2f} / {sim.t_max:.1f} s\n"
        f"---------------------------\n"
        f"State:     Pendulum 1   Pendulum 2\n"
        f"Theta:    {sim.sol_a['theta1'][i]:7.2f} rad {sim.sol_a['theta2'][i]:7.2f} rad\n"
        f"Omega:    {sim.sol_a['omega1'][i]:7.2f} r/s {sim.sol_a['omega2'][i]:7.2f} r/s\n"
        f"Inertia:  {I1:7.2f} kgm²{I2:7.2f} kgm²\n"
        f"---------------------------\n"
        f"Sys Inertia (I_sys): {I_sys:.2f} kg·m²\n"
        f"---------------------------\n"
        f"Chaos Separation:    {diff_pos:.4f} m\n"
        f"Lyapunov Est:        {lyap_est:.3f}\n"
        f"---------------------------\n"
        f"Mech. Energy (E):    {current_E:.2f} J\n"
        f"(E = K + U)"
    )
    
    if not sim.is_running:
        results_str += "\n[FINISHED]"
    elif sim.is_paused:
        results_str += "\n[PAUSED]"
        
    stats_text.set_text(results_str)
    
    if sim.is_running and not sim.is_paused:
        sim.current_frame += 1
    
    return line_main, line_chaos_a, line_chaos_b, line_angle1, line_angle2, line_E, line_K, line_U, point_phase1, point_phase2, line_sep, stats_text

ani = FuncAnimation(fig, animate, interval=20, blit=False, cache_frame_data=False)
plt.show()