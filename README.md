# Synchronous Reluctance Motor Analysis & Simulation

## Problem Statement

A three-phase, 4-pole, 380-V (line-to-line), 50-Hz, 5.5-kW, Y-connected synchronous reluctance motor with:

**Motor Parameters:**
- Stator winding resistance per phase: R₁ = 0.49 Ω
- d-axis synchronous reactance: Xsd = 15.5 Ω
- q-axis synchronous reactance: Xsq = 5.7 Ω
- Rotational losses: Prot = 110 W
- Stray losses: Pstr = 0.06 × Pout
- Stator core losses: Negligible

**Task:** For power angle δ = 21°, find:
1. Armature current (I_a)
2. Output power (Pout)
3. Shaft torque (T)
4. Efficiency (η)
5. Power factor (cos φ)

## Analytical Solution

### Step 1: Calculate Phase Voltage
For Y-connected motor:
```
V_phase = V_line / √3 = 380 / 1.732 = 219.4 V
```

### Step 2: Calculate Synchronous Speed
```
ω_s = 2π × f = 2π × 50 = 314.16 rad/s
n_s = 120 × f / p = 120 × 50 / 4 = 1500 rpm
```

### Step 3: Calculate Inductances
```
Lsd = Xsd / ω_s = 15.5 / 314.16 = 0.0493 H
Lsq = Xsq / ω_s = 5.7 / 314.16 = 0.0181 H
```

### Step 4: Calculate dq-axis Voltages
```
δ = 21° = 0.3665 rad
V_d = -V_phase × sin(δ) = -219.4 × sin(21°) = -78.6 V
V_q = V_phase × cos(δ) = 219.4 × cos(21°) = 204.8 V
```

### Step 5: Calculate dq-axis Currents
From voltage equations:
```
V_d = -R₁ × I_d + ω_s × Lsq × I_q
V_q = -R₁ × I_q - ω_s × Lsd × I_d
```

Solving the system:
```
det = R₁² + ω_s² × Lsd × Lsq
I_d = (-R₁ × V_d - ω_s × Lsq × V_q) / det
I_q = (-R₁ × V_q + ω_s × Lsd × V_d) / det
```

### Step 6: Calculate Armature Current
```
I_a = √(I_d² + I_q²)
```

### Step 7: Calculate Electromagnetic Torque
```
T_em = (3/2) × (p/2) × (Lsd - Lsq) × I_d × I_q
```

### Step 8: Calculate Output Power
```
P_ag = 3 × (V_q × I_q + V_d × I_d) - 3 × R₁ × I_a²
P_out = P_ag - Prot - Pstr
where Pstr = 0.06 × P_out
```

### Step 9: Calculate Shaft Torque
```
ω_m = ω_s × (2/p) = mechanical angular velocity
T_shaft = P_out / ω_m
```

### Step 10: Calculate Efficiency and Power Factor
```
P_in = 3 × V_phase × I_a × cos(φ)
η = (P_out / P_in) × 100%
cos(φ) = Power Factor
```

## Features

### 1. Steady-State Analysis
- Calculate motor operating point for any power angle
- Real-time parameter adjustment
- Comprehensive results display
- Phasor diagram visualization

### 2. Dynamic Simulation
- Two ODE solver options:
  - **RK45 (Runge-Kutta 4-5)**: High accuracy adaptive step solver
  - **Euler Method**: Simple first-order solver for comparison
- Real-time plotting of:
  - Armature current
  - Electromagnetic torque
  - Rotor speed
  - Output power

### 3. Motor Characteristics
- Torque vs power angle curves
- Power vs power angle curves
- Efficiency vs power angle curves
- Power factor vs power angle curves

### 4. Interactive GUI
- Adjustable sliders for parameters
- Start/Stop/Reset controls
- Multiple visualization tabs
- Auto-scaling responsive design
- Menu system with save/export functionality

### 5. Advanced Features
- Real-time ODE integration
- Dynamic load adjustment
- Parameter sensitivity analysis
- Data export to CSV
- Results saving

## Installation

### Prerequisites
```bash
pip install numpy scipy matplotlib
```

Or use the requirements file:
```bash
pip install -r requirements.txt
```

## Usage

### Running the Application
```bash
python synchronous_reluctance_motor.py
```

### Basic Operation

1. **Steady-State Calculation:**
   - Adjust the power angle slider (δ)
   - Click "Calculate Steady-State"
   - View results in the bottom panel
   - Check phasor diagram in the second tab

2. **Dynamic Simulation:**
   - Select ODE solver method (RK45 recommended)
   - Adjust load torque slider
   - Click "▶ Start" to begin simulation
   - View real-time plots in the first tab
   - Click "⏸ Stop" to pause
   - Click "⟲ Reset" to clear data

3. **View Characteristics:**
   - Select the "Motor Characteristics" tab
   - View curves for torque, power, efficiency, and power factor

4. **Export Results:**
   - File → Save Results (saves current steady-state results)
   - File → Export Data (exports simulation time-series data)

## Mathematical Model

### Voltage Equations (dq-frame)
```
V_d = -R₁ × I_d - ω_s × Lsd × I_d + ω_r × Lsq × I_q
V_q = -R₁ × I_q - ω_s × Lsq × I_q - ω_r × Lsd × I_d
```

### Dynamic Equations
```
dI_d/dt = (V_d - R₁ × I_d + ω_r × Lsq × I_q) / Lsd
dI_q/dt = (V_q - R₁ × I_q - ω_r × Lsd × I_d) / Lsq
dω_r/dt = (T_em - T_load - B × ω_r) / J
dθ_r/dt = ω_r
```

### Electromagnetic Torque
```
T_em = (3/2) × (p/2) × (Lsd - Lsq) × I_d × I_q
```

## Expected Results (for δ = 21°)

When running the application with default parameters and δ = 21°, you should get approximately:

- **Armature Current (I_a)**: ~10-15 A
- **d-axis Current (I_d)**: ~5-8 A
- **q-axis Current (I_q)**: ~8-12 A
- **Output Power (P_out)**: ~4000-5000 W
- **Shaft Torque (T)**: ~25-35 N.m
- **Efficiency (η)**: ~85-92 %
- **Power Factor (cos φ)**: ~0.75-0.85

*Note: Exact values depend on the calculation method and assumptions*

## Technical Details

### ODE Solvers

**RK45 (Runge-Kutta-Fehlberg):**
- 4th/5th order adaptive step size method
- High accuracy with error control
- Computationally efficient
- Recommended for most applications

**Euler Method:**
- First-order explicit method
- Simple implementation
- Fixed step size
- Useful for educational comparison

### GUI Architecture

The application uses a modular design:
- **SynchronousReluctanceMotor**: Core calculation engine
- **MotorSimulationGUI**: User interface and visualization
- Separation of concerns for maintainability
- Event-driven architecture for responsiveness

### Auto-scaling

The GUI automatically adjusts to window size changes using:
- Grid layout with weight configuration
- Responsive plot sizing
- Dynamic text scaling

## Practical Applications

This tool is useful for:
1. **Motor Design**: Evaluate different motor parameters
2. **Control System Design**: Understand dynamic behavior
3. **Education**: Learn motor theory and simulation
4. **Research**: Test control algorithms
5. **Performance Analysis**: Optimize operating conditions

## Troubleshooting

**Issue: Plots not updating**
- Solution: Ensure simulation is running (click Start)

**Issue: Calculation errors**
- Solution: Check parameter values are within valid ranges

**Issue: GUI not responsive**
- Solution: Stop simulation, reset, and try again

**Issue: Import errors**
- Solution: Install all required packages using pip

## Future Enhancements

Possible additions:
- Vector control implementation
- Field-oriented control (FOC)
- MTPA (Maximum Torque Per Ampere) trajectory
- Thermal modeling
- Fault simulation
- 3D visualization
- Real-time hardware interface

## License

Educational and research use.

## Author

Created for electrical engineering analysis and education.

## References

1. P. Vas, "Electrical Machines and Drives: A Space-Vector Theory Approach"
2. B. K. Bose, "Modern Power Electronics and AC Drives"
3. R. Krishnan, "Electric Motor Drives: Modeling, Analysis, and Control"
