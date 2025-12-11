import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
import pendulum_model  # Uses your existing physics engine

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
        
        self.sol_a = {} 
        self.sol_b = {} 
        self.separation = [] 
        self.perturbation = 0.001 
        
        self.run_sim()

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
        
        return {
            "t": t_arr, "x1": x1, "y1": y1, "x2": x2, "y2": y2,
            "theta1": theta1, "omega1": omega1,
            "theta2": theta2, "omega2": omega2
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

sim = SimulationManager()

# --- 2. Dashboard Layout ---
fig = plt.figure(figsize=(16, 9))

# Adjusted Plot Positions to prevent overlap
# Left plots moved from 0.33 -> 0.38
# Right plots moved from 0.69 -> 0.70
# Width reduced slightly 0.28 -> 0.27

# A. TOP ROW
ax_main = plt.axes([0.38, 0.55, 0.27, 0.40]) 
ax_main.set_title("Single Trajectory Simulation")
ax_main.set_xlabel("x (m)") 
ax_main.set_ylabel("y (m)") 
ax_main.set_aspect('equal')
ax_main.grid(True)

ax_chaos = plt.axes([0.70, 0.55, 0.27, 0.40])
ax_chaos.set_title(f"Chaos Demo (Diff: {sim.perturbation} rad)")
ax_chaos.set_xlabel("x (m)") 
ax_chaos.set_ylabel("y (m)") 
ax_chaos.set_aspect('equal')
ax_chaos.grid(True)

# B. BOTTOM ROW
# Merged Phase Space
ax_phase = plt.axes([0.38, 0.08, 0.27, 0.35])
ax_phase.set_title("Phase Space (Merged)")
ax_phase.set_xlabel("Angle (rad)")
ax_phase.set_ylabel("Angular Velocity (rad/s)")
ax_phase.grid(True)

# Separation Graph
ax_sep = plt.axes([0.70, 0.08, 0.27, 0.35])
ax_sep.set_title("Separation vs Time")
ax_sep.set_xlabel("Time (s)")
ax_sep.set_ylabel("Distance (m)")
ax_sep.grid(True)

# C. LEFT COLUMN (Results)
stats_text = fig.text(0.05, 0.30, '', va='top', ha='left', fontsize=11, 
                      fontfamily='monospace', bbox=dict(facecolor='white', alpha=0.5))
fig.text(0.05, 0.32, "--- SIMULATION RESULTS ---", fontsize=12, fontweight='bold')

# -- Visual Elements --
line_main, = ax_main.plot([], [], 'o-', lw=3, color='black')
line_chaos_a, = ax_chaos.plot([], [], 'o-', lw=2, color='blue', alpha=0.6, label='A')
line_chaos_b, = ax_chaos.plot([], [], 'o-', lw=2, color='red', alpha=0.6, label='B')
ax_chaos.legend(loc='upper right', fontsize='small')

point_phase1, = ax_phase.plot([], [], 'o', color='green', label='Pendulum 1')
point_phase2, = ax_phase.plot([], [], 'o', color='purple', label='Pendulum 2')
ax_phase.legend(loc='upper right', fontsize='small')

line_sep, = ax_sep.plot([], [], '-', lw=2, color='orange')

# --- 3. Control Panel with Layout Shifted Left ---
fig.text(0.05, 0.96, "CONTROL PANEL", fontsize=14, fontweight='bold')
fig.text(0.05, 0.93, "Parameters:", fontsize=10, fontweight='bold')

class SmartSlider:
    def __init__(self, fig, label, y_pos, vmin, vmax, valinit, step_val):
        self.step_val = step_val
        self.valinit = valinit
        
        # Coordinates shifted LEFT to create gap between sliders and plots
        # Previous X for dec button was 0.13 -> Now 0.10
        
        # 1. Decrease Button (<)
        ax_dec = plt.axes([0.10, y_pos, 0.015, 0.025]) 
        self.b_dec = Button(ax_dec, '<', color='lightgray', hovercolor='0.9')
        self.b_dec.on_clicked(self.decrement)
        
        # 2. Main Slider (Middle)
        ax_slide = plt.axes([0.12, y_pos, 0.10, 0.025], facecolor='lightgoldenrodyellow')
        
        # Hide default text
        self.slider = Slider(ax_slide, "", vmin, vmax, valinit=valinit)
        self.slider.valtext.set_visible(False)
        
        # Label (Left of button)
        fig.text(0.095, y_pos + 0.005, label, ha='right', va='center', fontsize=9)
        
        # 3. Increase Button (>)
        ax_inc = plt.axes([0.225, y_pos, 0.015, 0.025]) 
        self.b_inc = Button(ax_inc, '>', color='lightgray', hovercolor='0.9')
        self.b_inc.on_clicked(self.increment)
        
        # 4. Custom Value Text (Right of button)
        self.val_text = fig.text(0.25, y_pos + 0.005, f"{valinit:.3f}", 
                                 ha='left', va='center', fontsize=9, color='blue')

        self.slider.on_changed(self.update_text)
        
    def update_text(self, val):
        self.val_text.set_text(f"{val:.3f}")

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

# Create Smart Sliders
# Time: step 0.01
s_time = SmartSlider(fig, 'Time (s)',      0.90, 5.0, 40.0, sim.t_max, 0.01)

# Physics: step 0.1
s_g    = SmartSlider(fig, 'Gravity',       0.87, 1.0, 20.0, sim.params['g'], 0.1)
s_m1   = SmartSlider(fig, 'Mass 1',        0.84, 0.1, 5.0,  sim.params['m1'], 0.1)
s_m2   = SmartSlider(fig, 'Mass 2',        0.81, 0.1, 5.0,  sim.params['m2'], 0.1)
s_l1   = SmartSlider(fig, 'Len 1',         0.78, 0.5, 3.0,  sim.params['l1'], 0.1)
s_l2   = SmartSlider(fig, 'Len 2',         0.75, 0.5, 3.0,  sim.params['l2'], 0.1)

# Initial Conditions Group
fig.text(0.05, 0.72, "Initial Conditions:", fontsize=10, fontweight='bold')

# Angles: step 0.001
s_th1  = SmartSlider(fig, 'Theta 1',       0.69, -np.pi, np.pi, sim.initial_state[0], 0.001)
s_th2  = SmartSlider(fig, 'Theta 2',       0.66, -np.pi, np.pi, sim.initial_state[2], 0.001)

# Omega: step 0.01
s_om1  = SmartSlider(fig, 'Omega 1',       0.63, -10.0, 10.0,   sim.initial_state[1], 0.01)
s_om2  = SmartSlider(fig, 'Omega 2',       0.60, -10.0, 10.0,   sim.initial_state[3], 0.01)

# Buttons
ax_replay = plt.axes([0.05, 0.50, 0.09, 0.05])
b_replay = Button(ax_replay, 'Replay', hovercolor='0.975')

ax_reset = plt.axes([0.16, 0.50, 0.09, 0.05])
b_reset = Button(ax_reset, 'Reset', hovercolor='0.975', color='mistyrose')

# --- 4. Update Logic ---
smart_sliders = [s_time, s_g, s_m1, s_m2, s_l1, s_l2, s_th1, s_th2, s_om1, s_om2]

def update_sim(val=None):
    # 1. Update Params
    sim.t_max = s_time.val
    sim.params['m1'] = s_m1.val
    sim.params['m2'] = s_m2.val
    sim.params['l1'] = s_l1.val
    sim.params['l2'] = s_l2.val
    sim.params['g']  = s_g.val
    
    # 2. Update Initial State
    sim.initial_state = [s_th1.val, s_om1.val, s_th2.val, s_om2.val]
    
    # 3. Re-calculate
    sim.run_sim()
    
    # 4. Update Visual Limits
    limit = sim.params['l1'] + sim.params['l2'] + 0.2
    for ax in [ax_main, ax_chaos]:
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        
    all_theta = np.concatenate([sim.sol_a['theta1'], sim.sol_a['theta2']])
    all_omega = np.concatenate([sim.sol_a['omega1'], sim.sol_a['omega2']])
    
    # Safety check for flat-line signals
    if len(all_theta) > 0:
        ax_phase.set_xlim(np.min(all_theta) - 0.5, np.max(all_theta) + 0.5)
        ax_phase.set_ylim(np.min(all_omega) - 0.5, np.max(all_omega) + 0.5)
    
    ax_sep.set_xlim(0, sim.t_max)
    if len(sim.separation) > 0:
        ax_sep.set_ylim(0, np.max(sim.separation) * 1.1)

def replay(event):
    sim.current_frame = 0
    sim.is_running = True

def reset_all(event):
    for s in smart_sliders:
        s.reset()
    replay(event)

# Connect
for s in smart_sliders:
    s.on_changed(update_sim)

b_replay.on_clicked(replay)
b_reset.on_clicked(reset_all)

# Initial run
update_sim()

# --- 5. Animation Loop ---
def animate(frame):
    if not sim.is_running:
        return line_main, line_chaos_a, line_chaos_b, point_phase1, point_phase2, line_sep, stats_text

    i = sim.current_frame
    
    if i >= len(sim.sol_a['x1']) - 1:
        sim.is_running = False 
        i = len(sim.sol_a['x1']) - 1 
    
    # Main Sim
    x1, y1 = sim.sol_a['x1'][i], sim.sol_a['y1'][i]
    x2, y2 = sim.sol_a['x2'][i], sim.sol_a['y2'][i]
    line_main.set_data([0, x1, x2], [0, y1, y2])
    
    # Chaos Sim
    line_chaos_a.set_data([0, x1, x2], [0, y1, y2])
    xb1, yb1 = sim.sol_b['x1'][i], sim.sol_b['y1'][i]
    xb2, yb2 = sim.sol_b['x2'][i], sim.sol_b['y2'][i]
    line_chaos_b.set_data([0, xb1, xb2], [0, yb1, yb2])
    
    # Phase Space
    point_phase1.set_data([sim.sol_a['theta1'][i]], [sim.sol_a['omega1'][i]])
    point_phase2.set_data([sim.sol_a['theta2'][i]], [sim.sol_a['omega2'][i]])
    
    # Separation
    line_sep.set_data(sim.sol_a['t'][:i+1], sim.separation[:i+1])
    
    # Results
    diff_pos = sim.separation[i]
    results_str = (
        f"Time:       {sim.sol_a['t'][i]:.2f} / {sim.t_max:.1f} s\n\n"
        f"Theta 1:    {sim.sol_a['theta1'][i]:.2f} rad\n"
        f"Omega 1:    {sim.sol_a['omega1'][i]:.2f} rad/s\n\n"
        f"Theta 2:    {sim.sol_a['theta2'][i]:.2f} rad\n"
        f"Omega 2:    {sim.sol_a['omega2'][i]:.2f} rad/s\n\n"
        f"Chaos Sep:  {diff_pos:.4f} m\n"
    )
    if not sim.is_running:
        results_str += "\n[FINISHED]"
        
    stats_text.set_text(results_str)
    
    if sim.is_running:
        sim.current_frame += 1
    
    return line_main, line_chaos_a, line_chaos_b, point_phase1, point_phase2, line_sep, stats_text

ani = FuncAnimation(fig, animate, interval=20, blit=False, cache_frame_data=False)
plt.show()