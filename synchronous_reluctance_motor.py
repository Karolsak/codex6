"""
Synchronous Reluctance Motor Analysis and Simulation
Advanced Tkinter GUI Application for Electrical Engineering

Features:
- Motor parameter calculations
- Dynamic simulation with ODE solvers (RK45, Euler)
- Real-time visualization
- Interactive parameter adjustment
- Auto-scaling responsive design
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.integrate import solve_ivp
import math


class AlternatorAnalysis:
    """Utility class to study synchronous generator behavior."""

    def __init__(self, frequency=50.0):
        self.frequency = frequency

    def compute_parameters(self, e_open_circuit, i_short_circuit, r_a,
                            line_voltage, load_current, power_factor):
        """
        Calculate synchronous impedance, reactance and no-load voltage.

        Parameters
        ----------
        e_open_circuit : float
            Open-circuit generated line voltage (V).
        i_short_circuit : float
            Short-circuit armature current (A).
        r_a : float
            Armature resistance per phase (Ohm).
        line_voltage : float
            Rated line voltage at load (V).
        load_current : float
            Line current drawn by the load (A).
        power_factor : float
            Lagging power factor of the load (0 to 1).

        Returns
        -------
        dict
            Dictionary containing synchronous impedance, reactance and
            the expected no-load terminal voltage when the load is suddenly
            removed (same excitation).
        """
        z_s = e_open_circuit / i_short_circuit
        x_s = max(0.0, math.sqrt(max(z_s ** 2 - r_a ** 2, 0.0)))

        v_phase = line_voltage / math.sqrt(3)
        phi = math.acos(power_factor)
        i_phase = load_current
        i_complex = i_phase * (math.cos(-phi) + 1j * math.sin(-phi))

        e_phase = v_phase + (r_a + 1j * x_s) * i_complex
        v_no_load_line = abs(e_phase) * math.sqrt(3)
        voltage_regulation = (v_no_load_line - line_voltage) / line_voltage * 100

        return {
            "z_s": z_s,
            "x_s": x_s,
            "e_phase": e_phase,
            "no_load_voltage_line": v_no_load_line,
            "voltage_regulation_pct": voltage_regulation,
            "power_factor_angle_deg": math.degrees(phi)
        }

    def simulate_load_rejection(self, e_phase, r_a, x_s, i_initial,
                                phi_deg, duration, step, method="RK45"):
        """Simulate current decay and terminal voltage when the load is opened."""
        omega = 2 * math.pi * self.frequency
        l_s = x_s / omega if omega else 0.0

        def deriv(t, i):
            return [-i[0] * (r_a / l_s)] if l_s > 0 else [0]

        times = []
        currents = []
        voltages = []

        if method == "RK45":
            sol = solve_ivp(
                deriv, [0, duration], [i_initial], max_step=step, method="RK45"
            )
            i_values = sol.y[0]
            times = sol.t.tolist()
        else:
            i_values = []
            t = 0.0
            i_val = i_initial
            while t <= duration:
                i_values.append(i_val)
                times.append(t)
                if l_s > 0:
                    di_dt = -i_val * (r_a / l_s)
                    i_val += di_dt * step
                t += step

        angle = math.radians(-phi_deg)
        for i_val, t in zip(i_values, times):
            i_complex = i_val * (math.cos(angle) + 1j * math.sin(angle))
            v_phase = e_phase - (r_a + 1j * x_s) * i_complex
            voltages.append(abs(v_phase) * math.sqrt(3))
            currents.append(abs(i_complex))

        return times, currents, voltages


class SynchronousReluctanceMotor:
    """Mathematical model of a synchronous reluctance motor"""

    def __init__(self, V_line=380, poles=4, freq=50, P_rated=5500,
                 R1=0.49, Xsd=15.5, Xsq=5.7, Prot=110, Pstr_factor=0.06):
        """
        Initialize motor parameters

        Parameters:
        -----------
        V_line : float
            Line-to-line voltage (V)
        poles : int
            Number of poles
        freq : float
            Frequency (Hz)
        P_rated : float
            Rated power (W)
        R1 : float
            Stator resistance per phase (Ohm)
        Xsd : float
            d-axis synchronous reactance (Ohm)
        Xsq : float
            q-axis synchronous reactance (Ohm)
        Prot : float
            Rotational losses (W)
        Pstr_factor : float
            Stray losses factor (proportion of output power)
        """
        self.V_line = V_line
        self.poles = poles
        self.freq = freq
        self.P_rated = P_rated
        self.R1 = R1
        self.Xsd = Xsd
        self.Xsq = Xsq
        self.Prot = Prot
        self.Pstr_factor = Pstr_factor

        # Calculate derived parameters
        self.V_phase = V_line / np.sqrt(3)  # Phase voltage (Y-connected)
        self.omega_s = 2 * np.pi * freq  # Synchronous angular velocity (rad/s)
        self.n_s = 120 * freq / poles  # Synchronous speed (rpm)

        # Motor state variables
        self.delta = 0  # Power angle (degrees)
        self.I_a = 0  # Armature current (A)
        self.I_d = 0  # d-axis current (A)
        self.I_q = 0  # q-axis current (A)
        self.P_out = 0  # Output power (W)
        self.T_shaft = 0  # Shaft torque (N.m)
        self.efficiency = 0  # Efficiency (%)
        self.power_factor = 0  # Power factor
        self.phi = 0  # Phase angle (degrees)

    def calculate_steady_state(self, delta_deg):
        """
        Calculate steady-state operating point for a given power angle

        Parameters:
        -----------
        delta_deg : float
            Power angle in degrees

        Returns:
        --------
        dict : Dictionary containing all calculated parameters
        """
        self.delta = delta_deg
        delta_rad = np.deg2rad(delta_deg)

        # Voltage phasor components (assuming V is on q-axis reference)
        V_d = -self.V_phase * np.sin(delta_rad)
        V_q = self.V_phase * np.cos(delta_rad)

        # Calculate currents using dq-axis equations
        # Solving the equations:
        # V_d = -R1*I_d + omega_s*Lsq*I_q
        # V_q = -R1*I_q - omega_s*Lsd*I_d

        Lsd = self.Xsd / self.omega_s
        Lsq = self.Xsq / self.omega_s

        # Matrix solution for currents
        det = self.R1**2 + self.omega_s**2 * Lsd * Lsq

        self.I_d = (-self.R1 * V_d - self.omega_s * Lsq * V_q) / det
        self.I_q = (-self.R1 * V_q + self.omega_s * Lsd * V_d) / det

        # Total armature current
        self.I_a = np.sqrt(self.I_d**2 + self.I_q**2)

        # Electromagnetic torque
        T_em = 1.5 * self.poles / 2 * (Lsd - Lsq) * self.I_d * self.I_q

        # Copper losses
        P_cu = 3 * self.R1 * self.I_a**2

        # Air gap power
        P_ag = 3 * (V_q * self.I_q + V_d * self.I_d) - P_cu

        # Output power (after rotational losses)
        self.P_out = P_ag - self.Prot

        # Stray losses
        P_str = self.Pstr_factor * self.P_out if self.P_out > 0 else 0

        # Correct output power
        self.P_out = self.P_out - P_str

        # Shaft torque
        omega_m = self.omega_s * 2 / self.poles  # Mechanical angular velocity
        self.T_shaft = self.P_out / omega_m if omega_m > 0 else 0

        # Input power
        P_in = 3 * self.V_phase * self.I_a * np.cos(np.arctan2(self.I_d, self.I_q) - delta_rad)

        # Efficiency
        self.efficiency = (self.P_out / P_in * 100) if P_in > 0 else 0

        # Power factor
        self.phi = np.arctan2(self.I_d, self.I_q) - delta_rad
        self.power_factor = np.cos(self.phi)

        return {
            'delta': self.delta,
            'I_a': self.I_a,
            'I_d': self.I_d,
            'I_q': self.I_q,
            'P_out': self.P_out,
            'T_shaft': self.T_shaft,
            'efficiency': self.efficiency,
            'power_factor': self.power_factor,
            'phi_deg': np.rad2deg(self.phi),
            'T_em': T_em,
            'P_cu': P_cu,
            'P_in': P_in
        }

    def dynamic_model(self, t, y, T_load, V_phase):
        """
        Dynamic equations for the synchronous reluctance motor

        State variables:
        y[0] = i_d (d-axis current)
        y[1] = i_q (q-axis current)
        y[2] = omega_r (rotor angular velocity)
        y[3] = theta_r (rotor angle)

        Parameters:
        -----------
        t : float
            Time
        y : array
            State vector
        T_load : float
            Load torque (N.m)
        V_phase : float
            Phase voltage (V)

        Returns:
        --------
        dydt : array
            State derivatives
        """
        i_d, i_q, omega_r, theta_r = y

        Lsd = self.Xsd / self.omega_s
        Lsq = self.Xsq / self.omega_s

        # Voltage equations in dq frame
        v_d = V_phase * np.sin(self.omega_s * t - theta_r)
        v_q = V_phase * np.cos(self.omega_s * t - theta_r)

        # Current derivatives
        di_d_dt = (v_d - self.R1 * i_d + omega_r * Lsq * i_q) / Lsd
        di_q_dt = (v_q - self.R1 * i_q - omega_r * Lsd * i_d) / Lsq

        # Electromagnetic torque
        T_em = 1.5 * self.poles / 2 * (Lsd - Lsq) * i_d * i_q

        # Mechanical equation (simplified, assuming small inertia)
        J = 0.01  # Moment of inertia (kg.m^2) - typical value
        B = 0.001  # Friction coefficient

        domega_r_dt = (T_em - T_load - B * omega_r) / J
        dtheta_r_dt = omega_r

        return [di_d_dt, di_q_dt, domega_r_dt, dtheta_r_dt]


class MotorSimulationGUI:
    """Advanced Tkinter GUI for motor simulation"""

    def __init__(self, root):
        """Initialize the GUI application"""
        self.root = root
        self.root.title("Synchronous Reluctance Motor - Advanced Analysis & Simulation")
        self.root.geometry("1400x900")

        # Initialize motor
        self.motor = SynchronousReluctanceMotor()
        self.alternator = AlternatorAnalysis(frequency=self.motor.freq)

        # Simulation parameters
        self.is_running = False
        self.simulation_time = 0
        self.time_step = 0.001
        self.simulation_data = {
            'time': [],
            'current': [],
            'torque': [],
            'speed': [],
            'power': []
        }

        # ODE solver selection
        self.solver_method = tk.StringVar(value="RK45")

        # Initialize GUI components
        self.setup_menu()
        self.setup_gui()

        # Bind resize event for auto-scaling
        self.root.bind('<Configure>', self.on_window_resize)

        # Initial calculation
        self.calculate_steady_state()

    def setup_menu(self):
        """Create the main menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Save Results", command=self.save_results)
        file_menu.add_command(label="Export Data", command=self.export_data)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

        # Simulation menu
        sim_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Simulation", menu=sim_menu)
        sim_menu.add_command(label="Start", command=self.start_simulation)
        sim_menu.add_command(label="Stop", command=self.stop_simulation)
        sim_menu.add_command(label="Reset", command=self.reset_simulation)

        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_command(label="Steady-State Analysis", command=self.show_steady_state)
        view_menu.add_command(label="Dynamic Simulation", command=self.show_dynamic_sim)
        view_menu.add_command(label="Characteristics Curves", command=self.show_characteristics)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self.show_about)
        help_menu.add_command(label="User Guide", command=self.show_guide)

    def setup_gui(self):
        """Setup the main GUI layout"""
        # Create main container with grid weight for auto-resize
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        main_frame.grid_rowconfigure(1, weight=1)
        main_frame.grid_columnconfigure(1, weight=1)

        # Left panel - Input parameters and controls
        self.setup_left_panel(main_frame)

        # Right panel - Visualization
        self.setup_right_panel(main_frame)

        # Bottom panel - Results and status
        self.setup_bottom_panel(main_frame)

    def setup_left_panel(self, parent):
        """Setup the left panel with inputs and controls"""
        left_frame = ttk.LabelFrame(parent, text="Motor Parameters & Control", padding="10")
        left_frame.grid(row=0, column=0, rowspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)

        # Motor Parameters Section
        params_frame = ttk.LabelFrame(left_frame, text="Motor Parameters", padding="10")
        params_frame.pack(fill=tk.BOTH, expand=True, pady=5)

        # Create parameter entries with labels
        self.param_entries = {}
        parameters = [
            ("V_line", "Line Voltage (V)", self.motor.V_line, 200, 600),
            ("poles", "Number of Poles", self.motor.poles, 2, 12),
            ("freq", "Frequency (Hz)", self.motor.freq, 30, 100),
            ("P_rated", "Rated Power (W)", self.motor.P_rated, 1000, 15000),
            ("R1", "Stator Resistance (Ω)", self.motor.R1, 0.1, 2.0),
            ("Xsd", "Xsd - d-axis (Ω)", self.motor.Xsd, 5, 30),
            ("Xsq", "Xsq - q-axis (Ω)", self.motor.Xsq, 2, 20),
            ("Prot", "Rotational Losses (W)", self.motor.Prot, 0, 500)
        ]

        for i, (key, label, default, min_val, max_val) in enumerate(parameters):
            # Label
            ttk.Label(params_frame, text=label).grid(row=i, column=0, sticky=tk.W, pady=2)

            # Entry
            var = tk.DoubleVar(value=default)
            entry = ttk.Entry(params_frame, textvariable=var, width=15)
            entry.grid(row=i, column=1, padx=5, pady=2)
            self.param_entries[key] = var

            # Update button
            btn = ttk.Button(params_frame, text="Update",
                           command=lambda k=key: self.update_parameter(k))
            btn.grid(row=i, column=2, padx=5, pady=2)

        # Operating Point Section
        operating_frame = ttk.LabelFrame(left_frame, text="Operating Point", padding="10")
        operating_frame.pack(fill=tk.BOTH, expand=True, pady=5)

        # Power angle slider
        ttk.Label(operating_frame, text="Power Angle δ (degrees)").pack(anchor=tk.W)
        self.delta_var = tk.DoubleVar(value=21)
        self.delta_slider = ttk.Scale(operating_frame, from_=0, to=90,
                                     variable=self.delta_var, orient=tk.HORIZONTAL,
                                     command=self.on_delta_change)
        self.delta_slider.pack(fill=tk.X, pady=5)

        self.delta_label = ttk.Label(operating_frame, text="δ = 21.0°")
        self.delta_label.pack()

        # Load torque slider (for dynamic simulation)
        ttk.Label(operating_frame, text="Load Torque (N.m)").pack(anchor=tk.W, pady=(10, 0))
        self.torque_var = tk.DoubleVar(value=10)
        self.torque_slider = ttk.Scale(operating_frame, from_=0, to=50,
                                      variable=self.torque_var, orient=tk.HORIZONTAL)
        self.torque_slider.pack(fill=tk.X, pady=5)

        self.torque_label = ttk.Label(operating_frame, text="T_load = 10.0 N.m")
        self.torque_label.pack()

        # Solver Selection
        solver_frame = ttk.LabelFrame(left_frame, text="ODE Solver", padding="10")
        solver_frame.pack(fill=tk.X, pady=5)

        ttk.Radiobutton(solver_frame, text="RK45 (Runge-Kutta 4-5)",
                       variable=self.solver_method, value="RK45").pack(anchor=tk.W)
        ttk.Radiobutton(solver_frame, text="Euler Method",
                       variable=self.solver_method, value="Euler").pack(anchor=tk.W)

        # Control Buttons
        control_frame = ttk.LabelFrame(left_frame, text="Simulation Control", padding="10")
        control_frame.pack(fill=tk.X, pady=5)

        btn_frame = ttk.Frame(control_frame)
        btn_frame.pack(fill=tk.X)

        self.start_btn = ttk.Button(btn_frame, text="▶ Start", command=self.start_simulation)
        self.start_btn.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        self.stop_btn = ttk.Button(btn_frame, text="⏸ Stop", command=self.stop_simulation)
        self.stop_btn.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        self.reset_btn = ttk.Button(btn_frame, text="⟲ Reset", command=self.reset_simulation)
        self.reset_btn.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Calculate button for steady-state
        ttk.Button(control_frame, text="Calculate Steady-State",
                  command=self.calculate_steady_state).pack(fill=tk.X, pady=5)

    def setup_right_panel(self, parent):
        """Setup the right panel with visualization"""
        right_frame = ttk.LabelFrame(parent, text="Visualization", padding="10")
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        right_frame.grid_rowconfigure(0, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)

        # Create notebook for different views
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1: Real-time plots
        self.setup_realtime_plots()

        # Tab 2: Phasor diagram
        self.setup_phasor_diagram()

        # Tab 3: Characteristics
        self.setup_characteristics_plot()

        # Tab 4: Alternator Lab
        self.setup_alternator_tab()

        # Tab 5: Advanced Analysis
        self.setup_advanced_tab()

    def setup_realtime_plots(self):
        """Setup real-time plotting tab"""
        plot_frame = ttk.Frame(self.notebook)
        self.notebook.add(plot_frame, text="Real-Time Simulation")

        # Create figure with subplots
        self.fig_realtime = Figure(figsize=(8, 6), dpi=100)
        self.fig_realtime.subplots_adjust(hspace=0.4)

        self.ax_current = self.fig_realtime.add_subplot(221)
        self.ax_torque = self.fig_realtime.add_subplot(222)
        self.ax_speed = self.fig_realtime.add_subplot(223)
        self.ax_power = self.fig_realtime.add_subplot(224)

        self.ax_current.set_title('Armature Current')
        self.ax_current.set_xlabel('Time (s)')
        self.ax_current.set_ylabel('Current (A)')
        self.ax_current.grid(True)

        self.ax_torque.set_title('Electromagnetic Torque')
        self.ax_torque.set_xlabel('Time (s)')
        self.ax_torque.set_ylabel('Torque (N.m)')
        self.ax_torque.grid(True)

        self.ax_speed.set_title('Rotor Speed')
        self.ax_speed.set_xlabel('Time (s)')
        self.ax_speed.set_ylabel('Speed (rad/s)')
        self.ax_speed.grid(True)

        self.ax_power.set_title('Output Power')
        self.ax_power.set_xlabel('Time (s)')
        self.ax_power.set_ylabel('Power (W)')
        self.ax_power.grid(True)

        self.canvas_realtime = FigureCanvasTkAgg(self.fig_realtime, master=plot_frame)
        self.canvas_realtime.draw()
        self.canvas_realtime.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def setup_phasor_diagram(self):
        """Setup phasor diagram tab"""
        phasor_frame = ttk.Frame(self.notebook)
        self.notebook.add(phasor_frame, text="Phasor Diagram")

        self.fig_phasor = Figure(figsize=(8, 6), dpi=100)
        self.ax_phasor = self.fig_phasor.add_subplot(111)

        self.canvas_phasor = FigureCanvasTkAgg(self.fig_phasor, master=phasor_frame)
        self.canvas_phasor.draw()
        self.canvas_phasor.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def setup_characteristics_plot(self):
        """Setup motor characteristics plot tab"""
        char_frame = ttk.Frame(self.notebook)
        self.notebook.add(char_frame, text="Motor Characteristics")

        self.fig_char = Figure(figsize=(8, 6), dpi=100)
        self.fig_char.subplots_adjust(hspace=0.4)

        self.ax_torque_delta = self.fig_char.add_subplot(221)
        self.ax_power_delta = self.fig_char.add_subplot(222)
        self.ax_efficiency_delta = self.fig_char.add_subplot(223)
        self.ax_pf_delta = self.fig_char.add_subplot(224)

        self.canvas_char = FigureCanvasTkAgg(self.fig_char, master=char_frame)
        self.canvas_char.draw()

    def setup_alternator_tab(self):
        """Create alternator-focused calculation tab."""
        alt_frame = ttk.Frame(self.notebook)
        self.notebook.add(alt_frame, text="Alternator Lab")

        # Inputs
        input_frame = ttk.LabelFrame(alt_frame, text="Input Parameters", padding="10")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.alt_entries = {}
        alt_params = [
            ("e_oc", "Open-Circuit Voltage (V)", 1500.0),
            ("i_sc", "Short-Circuit Current (A)", 250.0),
            ("r_a", "Armature Resistance (Ω)", 2.0),
            ("v_line", "Rated Line Voltage (V)", 6600.0),
            ("i_load", "Load Current (A)", 250.0),
            ("pf", "Power Factor (lagging)", 0.8),
        ]

        for idx, (key, label, default) in enumerate(alt_params):
            row = ttk.Frame(input_frame)
            row.pack(fill=tk.X, pady=3)
            ttk.Label(row, text=label, width=28).pack(side=tk.LEFT)
            var = tk.DoubleVar(value=default)
            entry = ttk.Entry(row, textvariable=var, width=18)
            entry.pack(side=tk.LEFT, padx=5)
            self.alt_entries[key] = var

        ttk.Button(input_frame, text="Compute", command=self.compute_alternator).pack(
            fill=tk.X, pady=10
        )

        # Results
        result_frame = ttk.LabelFrame(alt_frame, text="Results", padding="10")
        result_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.alt_results = {}
        result_labels = [
            ("Synchronous Impedance |Zs| (Ω)", "z_s"),
            ("Synchronous Reactance Xs (Ω)", "x_s"),
            ("No-load Line Voltage (V)", "no_load_voltage_line"),
            ("Voltage Regulation (%)", "voltage_regulation_pct"),
            ("Power Factor Angle (deg)", "power_factor_angle_deg"),
        ]

        for label, key in result_labels:
            row = ttk.Frame(result_frame)
            row.pack(fill=tk.X, pady=2)
            ttk.Label(row, text=label + ":", width=28).pack(side=tk.LEFT)
            val_label = ttk.Label(row, text="-")
            val_label.pack(side=tk.LEFT, padx=5)
            self.alt_results[key] = val_label

        info = (
            "The alternator lab models synchronous impedance using open-circuit "
            "and short-circuit tests, then estimates the no-load terminal voltage "
            "after a sudden load rejection."
        )
        ttk.Label(result_frame, text=info, wraplength=360, foreground="gray").pack(
            fill=tk.X, pady=10
        )

        self.compute_alternator()

    def setup_advanced_tab(self):
        """Comprehensive analysis tab with voltage recovery simulation."""
        adv_frame = ttk.Frame(self.notebook)
        self.notebook.add(adv_frame, text="Advanced Analysis")

        control_frame = ttk.LabelFrame(adv_frame, text="Load Rejection Simulation", padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        ttk.Label(control_frame, text="Duration (s)").pack(anchor=tk.W)
        self.adv_duration = tk.DoubleVar(value=0.6)
        ttk.Scale(control_frame, from_=0.1, to=2.0, variable=self.adv_duration,
                  orient=tk.HORIZONTAL).pack(fill=tk.X, pady=3)

        ttk.Label(control_frame, text="Step (s)").pack(anchor=tk.W)
        self.adv_step = tk.DoubleVar(value=0.01)
        ttk.Scale(control_frame, from_=0.001, to=0.05, variable=self.adv_step,
                  orient=tk.HORIZONTAL).pack(fill=tk.X, pady=3)

        ttk.Button(control_frame, text="Run Load-Shedding Study",
                  command=self.run_advanced_study).pack(fill=tk.X, pady=8)

        ttk.Label(control_frame, text="Solver").pack(anchor=tk.W, pady=(8, 0))
        ttk.Radiobutton(control_frame, text="RK45", variable=self.solver_method,
                       value="RK45").pack(anchor=tk.W)
        ttk.Radiobutton(control_frame, text="Euler", variable=self.solver_method,
                       value="Euler").pack(anchor=tk.W)

        plot_frame = ttk.LabelFrame(adv_frame, text="Dynamic Response", padding="10")
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.fig_adv = Figure(figsize=(7, 5), dpi=100)
        self.fig_adv.subplots_adjust(hspace=0.35)
        self.ax_adv_voltage = self.fig_adv.add_subplot(211)
        self.ax_adv_current = self.fig_adv.add_subplot(212)

        self.ax_adv_voltage.set_ylabel("Voltage (V)")
        self.ax_adv_current.set_ylabel("Current (A)")
        self.ax_adv_current.set_xlabel("Time (s)")

        self.canvas_adv = FigureCanvasTkAgg(self.fig_adv, master=plot_frame)
        self.canvas_adv.draw()
        self.canvas_adv.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def compute_alternator(self):
        """Compute alternator parameters and update labels."""
        e_oc = self.alt_entries["e_oc"].get()
        i_sc = self.alt_entries["i_sc"].get()
        r_a = self.alt_entries["r_a"].get()
        v_line = self.alt_entries["v_line"].get()
        i_load = self.alt_entries["i_load"].get()
        pf = self.alt_entries["pf"].get()

        results = self.alternator.compute_parameters(e_oc, i_sc, r_a, v_line, i_load, pf)

        for key, label in self.alt_results.items():
            value = results.get(key, 0)
            if isinstance(value, complex):
                text = f"{abs(value):.2f} ∠ {math.degrees(math.atan2(value.imag, value.real)):.1f}°"
            else:
                text = f"{value:.3f}" if abs(value) < 1e4 else f"{value:,.1f}"
            label.config(text=text)

        self.latest_alternator = results
        self.status_var.set("Alternator parameters updated")

    def run_advanced_study(self):
        """Run voltage recovery simulation for the alternator."""
        if not hasattr(self, "latest_alternator"):
            self.compute_alternator()

        r_a = self.alt_entries["r_a"].get()
        i_load = self.alt_entries["i_load"].get()
        e_phase = self.latest_alternator.get("e_phase", 0)
        x_s = self.latest_alternator.get("x_s", 0)
        phi_deg = self.latest_alternator.get("power_factor_angle_deg", 0)

        duration = self.adv_duration.get()
        step = self.adv_step.get()

        times, currents, voltages = self.alternator.simulate_load_rejection(
            e_phase, r_a, x_s, i_load, phi_deg, duration, step, method=self.solver_method.get()
        )

        self.ax_adv_voltage.clear()
        self.ax_adv_current.clear()
        self.ax_adv_voltage.plot(times, voltages, 'b-', label='Terminal Voltage (line)')
        self.ax_adv_current.plot(times, currents, 'r-', label='Armature Current')

        self.ax_adv_voltage.set_title("Load Rejection Voltage Recovery")
        self.ax_adv_voltage.set_ylabel("Voltage (V)")
        self.ax_adv_voltage.grid(True, alpha=0.3)
        self.ax_adv_voltage.legend()

        self.ax_adv_current.set_title("Current Decay")
        self.ax_adv_current.set_ylabel("Current (A)")
        self.ax_adv_current.set_xlabel("Time (s)")
        self.ax_adv_current.grid(True, alpha=0.3)
        self.ax_adv_current.legend()

        self.fig_adv.tight_layout()
        self.canvas_adv.draw()
        self.status_var.set("Advanced voltage recovery study completed")
        self.canvas_char.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initial plot
        self.plot_characteristics()

    def setup_bottom_panel(self, parent):
        """Setup the bottom panel with results"""
        bottom_frame = ttk.LabelFrame(parent, text="Results & Status", padding="10")
        bottom_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)

        # Create results display
        results_container = ttk.Frame(bottom_frame)
        results_container.pack(fill=tk.BOTH, expand=True)

        # Split into two columns
        left_results = ttk.Frame(results_container)
        left_results.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        right_results = ttk.Frame(results_container)
        right_results.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        # Result labels
        self.result_labels = {}

        results_left = [
            ("Power Angle δ", "delta"),
            ("Armature Current I_a", "I_a"),
            ("d-axis Current I_d", "I_d"),
            ("q-axis Current I_q", "I_q")
        ]

        results_right = [
            ("Output Power P_out", "P_out"),
            ("Shaft Torque T", "T_shaft"),
            ("Efficiency η", "efficiency"),
            ("Power Factor cos φ", "power_factor")
        ]

        for i, (label, key) in enumerate(results_left):
            frame = ttk.Frame(left_results)
            frame.pack(fill=tk.X, pady=2)
            ttk.Label(frame, text=label + ":", width=20).pack(side=tk.LEFT)
            lbl = ttk.Label(frame, text="--", font=('TkDefaultFont', 10, 'bold'))
            lbl.pack(side=tk.LEFT)
            self.result_labels[key] = lbl

        for i, (label, key) in enumerate(results_right):
            frame = ttk.Frame(right_results)
            frame.pack(fill=tk.X, pady=2)
            ttk.Label(frame, text=label + ":", width=20).pack(side=tk.LEFT)
            lbl = ttk.Label(frame, text="--", font=('TkDefaultFont', 10, 'bold'))
            lbl.pack(side=tk.LEFT)
            self.result_labels[key] = lbl

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(bottom_frame, textvariable=self.status_var,
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(10, 0))

    def update_parameter(self, param_name):
        """Update motor parameter"""
        try:
            value = self.param_entries[param_name].get()
            setattr(self.motor, param_name, value)

            # Recalculate derived parameters
            self.motor.V_phase = self.motor.V_line / np.sqrt(3)
            self.motor.omega_s = 2 * np.pi * self.motor.freq
            self.motor.n_s = 120 * self.motor.freq / self.motor.poles

            self.calculate_steady_state()
            self.status_var.set(f"Updated {param_name} = {value}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to update parameter: {str(e)}")

    def on_delta_change(self, value):
        """Handle power angle slider change"""
        delta = float(value)
        self.delta_label.config(text=f"δ = {delta:.1f}°")
        self.calculate_steady_state()

    def calculate_steady_state(self):
        """Calculate and display steady-state results"""
        try:
            delta = self.delta_var.get()
            results = self.motor.calculate_steady_state(delta)

            # Update result labels
            self.result_labels['delta'].config(text=f"{results['delta']:.2f}°")
            self.result_labels['I_a'].config(text=f"{results['I_a']:.3f} A")
            self.result_labels['I_d'].config(text=f"{results['I_d']:.3f} A")
            self.result_labels['I_q'].config(text=f"{results['I_q']:.3f} A")
            self.result_labels['P_out'].config(text=f"{results['P_out']:.2f} W")
            self.result_labels['T_shaft'].config(text=f"{results['T_shaft']:.3f} N.m")
            self.result_labels['efficiency'].config(text=f"{results['efficiency']:.2f} %")
            self.result_labels['power_factor'].config(text=f"{results['power_factor']:.3f}")

            # Update phasor diagram
            self.update_phasor_diagram(results)

            self.status_var.set("Steady-state calculation completed")
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {str(e)}")

    def update_phasor_diagram(self, results):
        """Update the phasor diagram"""
        self.ax_phasor.clear()

        # Get parameters
        delta_rad = np.deg2rad(results['delta'])
        phi_rad = np.deg2rad(results['phi_deg'])

        # Voltage phasor (reference)
        V = self.motor.V_phase
        V_angle = delta_rad

        # Current phasor
        I = results['I_a']
        I_angle = phi_rad

        # Plot phasors
        self.ax_phasor.quiver(0, 0, V * np.cos(V_angle), V * np.sin(V_angle),
                             angles='xy', scale_units='xy', scale=1, color='blue',
                             width=0.005, label=f'V = {V:.1f}V')

        self.ax_phasor.quiver(0, 0, I * 20 * np.cos(I_angle), I * 20 * np.sin(I_angle),
                             angles='xy', scale_units='xy', scale=1, color='red',
                             width=0.005, label=f'I = {I:.2f}A (×20)')

        # d-q axes
        self.ax_phasor.plot([0, 100], [0, 0], 'k--', alpha=0.3, label='q-axis')
        self.ax_phasor.plot([0, 0], [0, 100], 'k--', alpha=0.3, label='d-axis')

        self.ax_phasor.set_xlim(-150, 150)
        self.ax_phasor.set_ylim(-150, 150)
        self.ax_phasor.set_aspect('equal')
        self.ax_phasor.grid(True, alpha=0.3)
        self.ax_phasor.legend()
        self.ax_phasor.set_title('Phasor Diagram (dq reference frame)')
        self.ax_phasor.set_xlabel('q-axis')
        self.ax_phasor.set_ylabel('d-axis')

        self.canvas_phasor.draw()

    def plot_characteristics(self):
        """Plot motor characteristics vs power angle"""
        delta_range = np.linspace(1, 89, 100)

        torques = []
        powers = []
        efficiencies = []
        power_factors = []

        for delta in delta_range:
            try:
                results = self.motor.calculate_steady_state(delta)
                torques.append(results['T_shaft'])
                powers.append(results['P_out'])
                efficiencies.append(results['efficiency'])
                power_factors.append(results['power_factor'])
            except:
                torques.append(0)
                powers.append(0)
                efficiencies.append(0)
                power_factors.append(0)

        # Plot torque vs delta
        self.ax_torque_delta.clear()
        self.ax_torque_delta.plot(delta_range, torques, 'b-', linewidth=2)
        self.ax_torque_delta.set_xlabel('Power Angle δ (degrees)')
        self.ax_torque_delta.set_ylabel('Torque (N.m)')
        self.ax_torque_delta.set_title('Torque vs Power Angle')
        self.ax_torque_delta.grid(True, alpha=0.3)

        # Plot power vs delta
        self.ax_power_delta.clear()
        self.ax_power_delta.plot(delta_range, powers, 'r-', linewidth=2)
        self.ax_power_delta.set_xlabel('Power Angle δ (degrees)')
        self.ax_power_delta.set_ylabel('Power (W)')
        self.ax_power_delta.set_title('Output Power vs Power Angle')
        self.ax_power_delta.grid(True, alpha=0.3)

        # Plot efficiency vs delta
        self.ax_efficiency_delta.clear()
        self.ax_efficiency_delta.plot(delta_range, efficiencies, 'g-', linewidth=2)
        self.ax_efficiency_delta.set_xlabel('Power Angle δ (degrees)')
        self.ax_efficiency_delta.set_ylabel('Efficiency (%)')
        self.ax_efficiency_delta.set_title('Efficiency vs Power Angle')
        self.ax_efficiency_delta.grid(True, alpha=0.3)

        # Plot power factor vs delta
        self.ax_pf_delta.clear()
        self.ax_pf_delta.plot(delta_range, power_factors, 'm-', linewidth=2)
        self.ax_pf_delta.set_xlabel('Power Angle δ (degrees)')
        self.ax_pf_delta.set_ylabel('Power Factor')
        self.ax_pf_delta.set_title('Power Factor vs Power Angle')
        self.ax_pf_delta.grid(True, alpha=0.3)

        self.canvas_char.draw()

    def start_simulation(self):
        """Start dynamic simulation"""
        if not self.is_running:
            self.is_running = True
            self.status_var.set("Simulation running...")
            self.start_btn.config(state='disabled')
            self.run_dynamic_simulation()

    def stop_simulation(self):
        """Stop dynamic simulation"""
        self.is_running = False
        self.status_var.set("Simulation stopped")
        self.start_btn.config(state='normal')

    def reset_simulation(self):
        """Reset simulation data"""
        self.stop_simulation()
        self.simulation_time = 0
        self.simulation_data = {
            'time': [],
            'current': [],
            'torque': [],
            'speed': [],
            'power': []
        }

        # Clear plots
        self.ax_current.clear()
        self.ax_torque.clear()
        self.ax_speed.clear()
        self.ax_power.clear()

        self.ax_current.set_title('Armature Current')
        self.ax_current.set_xlabel('Time (s)')
        self.ax_current.set_ylabel('Current (A)')
        self.ax_current.grid(True)

        self.ax_torque.set_title('Electromagnetic Torque')
        self.ax_torque.set_xlabel('Time (s)')
        self.ax_torque.set_ylabel('Torque (N.m)')
        self.ax_torque.grid(True)

        self.ax_speed.set_title('Rotor Speed')
        self.ax_speed.set_xlabel('Time (s)')
        self.ax_speed.set_ylabel('Speed (rad/s)')
        self.ax_speed.grid(True)

        self.ax_power.set_title('Output Power')
        self.ax_power.set_xlabel('Time (s)')
        self.ax_power.set_ylabel('Power (W)')
        self.ax_power.grid(True)

        self.canvas_realtime.draw()
        self.status_var.set("Simulation reset")

    def run_dynamic_simulation(self):
        """Run dynamic simulation with ODE solver"""
        if not self.is_running:
            return

        T_load = self.torque_var.get()
        self.torque_label.config(text=f"T_load = {T_load:.1f} N.m")

        # Initial conditions
        if len(self.simulation_data['time']) == 0:
            y0 = [0, 0, 0, 0]  # [i_d, i_q, omega_r, theta_r]
            t0 = 0
        else:
            # Continue from last state (simplified)
            y0 = [self.motor.I_d, self.motor.I_q, self.motor.omega_s * 2 / self.motor.poles, 0]
            t0 = self.simulation_data['time'][-1]

        t_span = [t0, t0 + 0.1]  # Simulate 0.1 seconds

        if self.solver_method.get() == "RK45":
            # Use scipy's RK45 solver
            sol = solve_ivp(
                lambda t, y: self.motor.dynamic_model(t, y, T_load, self.motor.V_phase),
                t_span, y0, method='RK45', max_step=self.time_step
            )

            # Extract data
            for i in range(len(sol.t)):
                i_d, i_q, omega_r, theta_r = sol.y[:, i]
                i_a = np.sqrt(i_d**2 + i_q**2)

                Lsd = self.motor.Xsd / self.motor.omega_s
                Lsq = self.motor.Xsq / self.motor.omega_s
                T_em = 1.5 * self.motor.poles / 2 * (Lsd - Lsq) * i_d * i_q
                P_out = T_em * omega_r

                self.simulation_data['time'].append(sol.t[i])
                self.simulation_data['current'].append(i_a)
                self.simulation_data['torque'].append(T_em)
                self.simulation_data['speed'].append(omega_r)
                self.simulation_data['power'].append(P_out)

        else:  # Euler method
            t = t0
            y = np.array(y0)
            dt = self.time_step

            for _ in range(int(0.1 / dt)):
                dydt = self.motor.dynamic_model(t, y, T_load, self.motor.V_phase)
                y = y + np.array(dydt) * dt
                t = t + dt

                i_d, i_q, omega_r, theta_r = y
                i_a = np.sqrt(i_d**2 + i_q**2)

                Lsd = self.motor.Xsd / self.motor.omega_s
                Lsq = self.motor.Xsq / self.motor.omega_s
                T_em = 1.5 * self.motor.poles / 2 * (Lsd - Lsq) * i_d * i_q
                P_out = T_em * omega_r

                self.simulation_data['time'].append(t)
                self.simulation_data['current'].append(i_a)
                self.simulation_data['torque'].append(T_em)
                self.simulation_data['speed'].append(omega_r)
                self.simulation_data['power'].append(P_out)

        # Update plots
        self.update_realtime_plots()

        # Continue simulation
        self.root.after(10, self.run_dynamic_simulation)

    def update_realtime_plots(self):
        """Update real-time plots"""
        if len(self.simulation_data['time']) > 1:
            # Limit data points for better performance
            max_points = 1000
            if len(self.simulation_data['time']) > max_points:
                for key in self.simulation_data:
                    self.simulation_data[key] = self.simulation_data[key][-max_points:]

            time = self.simulation_data['time']

            # Current plot
            self.ax_current.clear()
            self.ax_current.plot(time, self.simulation_data['current'], 'b-')
            self.ax_current.set_title('Armature Current')
            self.ax_current.set_xlabel('Time (s)')
            self.ax_current.set_ylabel('Current (A)')
            self.ax_current.grid(True, alpha=0.3)

            # Torque plot
            self.ax_torque.clear()
            self.ax_torque.plot(time, self.simulation_data['torque'], 'r-')
            self.ax_torque.set_title('Electromagnetic Torque')
            self.ax_torque.set_xlabel('Time (s)')
            self.ax_torque.set_ylabel('Torque (N.m)')
            self.ax_torque.grid(True, alpha=0.3)

            # Speed plot
            self.ax_speed.clear()
            self.ax_speed.plot(time, self.simulation_data['speed'], 'g-')
            self.ax_speed.set_title('Rotor Speed')
            self.ax_speed.set_xlabel('Time (s)')
            self.ax_speed.set_ylabel('Speed (rad/s)')
            self.ax_speed.grid(True, alpha=0.3)

            # Power plot
            self.ax_power.clear()
            self.ax_power.plot(time, self.simulation_data['power'], 'm-')
            self.ax_power.set_title('Output Power')
            self.ax_power.set_xlabel('Time (s)')
            self.ax_power.set_ylabel('Power (W)')
            self.ax_power.grid(True, alpha=0.3)

            self.canvas_realtime.draw()

    def on_window_resize(self, event):
        """Handle window resize event for auto-scaling"""
        figures = [self.fig_realtime, self.fig_phasor, self.fig_char]
        if hasattr(self, "fig_adv"):
            figures.append(self.fig_adv)

        for fig in figures:
            fig.tight_layout()

        if hasattr(self, "canvas_realtime"):
            self.canvas_realtime.draw_idle()
        if hasattr(self, "canvas_phasor"):
            self.canvas_phasor.draw_idle()
        if hasattr(self, "canvas_char"):
            self.canvas_char.draw_idle()
        if hasattr(self, "canvas_adv"):
            self.canvas_adv.draw_idle()

    def save_results(self):
        """Save results to file"""
        try:
            with open('motor_results.txt', 'w') as f:
                f.write("Synchronous Reluctance Motor Analysis Results\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Power Angle: {self.motor.delta:.2f}°\n")
                f.write(f"Armature Current: {self.motor.I_a:.3f} A\n")
                f.write(f"d-axis Current: {self.motor.I_d:.3f} A\n")
                f.write(f"q-axis Current: {self.motor.I_q:.3f} A\n")
                f.write(f"Output Power: {self.motor.P_out:.2f} W\n")
                f.write(f"Shaft Torque: {self.motor.T_shaft:.3f} N.m\n")
                f.write(f"Efficiency: {self.motor.efficiency:.2f} %\n")
                f.write(f"Power Factor: {self.motor.power_factor:.3f}\n")

            messagebox.showinfo("Success", "Results saved to motor_results.txt")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save results: {str(e)}")

    def export_data(self):
        """Export simulation data to CSV"""
        try:
            if len(self.simulation_data['time']) > 0:
                with open('simulation_data.csv', 'w') as f:
                    f.write("Time,Current,Torque,Speed,Power\n")
                    for i in range(len(self.simulation_data['time'])):
                        f.write(f"{self.simulation_data['time'][i]:.6f},")
                        f.write(f"{self.simulation_data['current'][i]:.6f},")
                        f.write(f"{self.simulation_data['torque'][i]:.6f},")
                        f.write(f"{self.simulation_data['speed'][i]:.6f},")
                        f.write(f"{self.simulation_data['power'][i]:.6f}\n")

                messagebox.showinfo("Success", "Data exported to simulation_data.csv")
            else:
                messagebox.showwarning("Warning", "No simulation data to export")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export data: {str(e)}")

    def show_steady_state(self):
        """Show steady-state analysis view"""
        self.notebook.select(1)  # Phasor diagram tab

    def show_dynamic_sim(self):
        """Show dynamic simulation view"""
        self.notebook.select(0)  # Real-time plots tab

    def show_characteristics(self):
        """Show characteristics curves view"""
        self.notebook.select(2)  # Characteristics tab

    def show_about(self):
        """Show about dialog"""
        about_text = """
Synchronous Reluctance Motor Analysis & Simulation
Version 1.0

Features:
• Steady-state analysis
• Dynamic simulation with RK45 and Euler ODE solvers
• Real-time visualization
• Interactive parameter adjustment
• Motor characteristics curves
• Phasor diagrams

Developed for Electrical Engineering Applications
        """
        messagebox.showinfo("About", about_text)

    def show_guide(self):
        """Show user guide"""
        guide_text = """
User Guide:

1. Steady-State Analysis:
   - Adjust power angle using the slider
   - Click "Calculate Steady-State" to update results
   - View phasor diagram in the second tab

2. Dynamic Simulation:
   - Select ODE solver (RK45 or Euler)
   - Adjust load torque
   - Click "Start" to begin simulation
   - View real-time plots in the first tab

3. Characteristics:
   - View motor characteristics vs power angle
   - Third tab shows torque, power, efficiency, and power factor curves

4. Export:
   - Save results using File > Save Results
   - Export simulation data using File > Export Data
        """
        messagebox.showinfo("User Guide", guide_text)


def main():
    """Main application entry point"""
    root = tk.Tk()
    app = MotorSimulationGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
