"""
Test script to verify motor calculations for δ = 21°
"""

import sys
import numpy as np

# Import the motor class (without GUI dependencies for testing)
class SynchronousReluctanceMotor:
    """Mathematical model of a synchronous reluctance motor"""

    def __init__(self, V_line=380, poles=4, freq=50, P_rated=5500,
                 R1=0.49, Xsd=15.5, Xsq=5.7, Prot=110, Pstr_factor=0.06):
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
        self.V_phase = V_line / np.sqrt(3)
        self.omega_s = 2 * np.pi * freq
        self.n_s = 120 * freq / poles

    def calculate_steady_state(self, delta_deg):
        delta_rad = np.deg2rad(delta_deg)

        # Reference: d-axis leads q-axis by 90 degrees (rotor reference frame)
        # Power angle: angle by which rotor d-axis leads stator voltage
        # For motor operation: V_q is positive (aligned with voltage), I_q > 0, I_d > 0

        # Voltage components in rotor dq frame
        V_d = -self.V_phase * np.sin(delta_rad)
        V_q = self.V_phase * np.cos(delta_rad)

        # Calculate inductances
        Lsd = self.Xsd / self.omega_s
        Lsq = self.Xsq / self.omega_s

        # Steady-state voltage equations:
        # V_d = -R1*I_d + omega_s*Lsq*I_q  (simplified, phasor balance)
        # V_q = -R1*I_q - omega_s*Lsd*I_d

        # Solving for currents:
        det = self.R1**2 + self.omega_s**2 * Lsd * Lsq

        I_d = (-self.R1 * V_d - self.omega_s * Lsq * V_q) / det
        I_q = (-self.R1 * V_q + self.omega_s * Lsd * V_d) / det

        # Total armature current magnitude
        I_a = np.sqrt(I_d**2 + I_q**2)

        # Electromagnetic torque (3-phase, reluctance torque)
        T_em = 1.5 * (self.poles / 2) * (Lsd - Lsq) * I_d * I_q

        # Input power (3-phase)
        P_in = 3 * (V_d * I_d + V_q * I_q)

        # Copper losses
        P_cu = 3 * self.R1 * I_a**2

        # Air gap power
        P_ag = P_in - P_cu

        # Mechanical angular velocity
        omega_m = self.omega_s * 2 / self.poles

        # Mechanical power
        P_mech = T_em * omega_m

        # Output power (subtract rotational losses)
        P_out_prelim = P_mech - self.Prot

        # Stray losses
        P_str = self.Pstr_factor * abs(P_out_prelim)

        # Final output power
        P_out = P_out_prelim - P_str if P_out_prelim > 0 else P_out_prelim + P_str

        # Shaft torque
        T_shaft = P_out / omega_m if omega_m != 0 else 0

        # Current phase angle relative to voltage
        I_phase = np.arctan2(I_d, I_q)
        V_phase_ang = delta_rad

        # Power factor
        phi_angle = I_phase - V_phase_ang
        power_factor = np.cos(phi_angle)

        # Efficiency
        efficiency = (abs(P_out) / abs(P_in) * 100) if P_in != 0 else 0

        return {
            'delta': delta_deg,
            'I_a': I_a,
            'I_d': I_d,
            'I_q': I_q,
            'P_out': P_out,
            'T_shaft': T_shaft,
            'efficiency': efficiency,
            'power_factor': power_factor,
            'phi_deg': np.rad2deg(phi_angle),
            'T_em': T_em,
            'P_cu': P_cu,
            'P_in': P_in,
            'V_phase': self.V_phase,
            'omega_s': self.omega_s,
            'n_s': self.n_s
        }


def main():
    print("=" * 70)
    print("SYNCHRONOUS RELUCTANCE MOTOR ANALYSIS")
    print("=" * 70)
    print()

    # Create motor with given parameters
    motor = SynchronousReluctanceMotor(
        V_line=380,      # V
        poles=4,
        freq=50,         # Hz
        P_rated=5500,    # W
        R1=0.49,         # Ohm
        Xsd=15.5,        # Ohm
        Xsq=5.7,         # Ohm
        Prot=110,        # W
        Pstr_factor=0.06
    )

    print("MOTOR PARAMETERS:")
    print("-" * 70)
    print(f"Line-to-line Voltage (V_line)      : {motor.V_line} V")
    print(f"Phase Voltage (V_phase)            : {motor.V_phase:.2f} V")
    print(f"Number of Poles                    : {motor.poles}")
    print(f"Frequency                          : {motor.freq} Hz")
    print(f"Synchronous Speed                  : {motor.n_s} rpm")
    print(f"Synchronous Angular Velocity       : {motor.omega_s:.2f} rad/s")
    print(f"Rated Power                        : {motor.P_rated} W")
    print(f"Stator Resistance (R1)             : {motor.R1} Ω")
    print(f"d-axis Synchronous Reactance (Xsd) : {motor.Xsd} Ω")
    print(f"q-axis Synchronous Reactance (Xsq) : {motor.Xsq} Ω")
    print(f"Rotational Losses                  : {motor.Prot} W")
    print(f"Stray Losses Factor                : {motor.Pstr_factor}")
    print()

    # Calculate for δ = 21°
    delta = 21
    print(f"STEADY-STATE ANALYSIS FOR δ = {delta}°")
    print("-" * 70)

    results = motor.calculate_steady_state(delta)

    print(f"\nOPERATING POINT:")
    print(f"  Power Angle (δ)                  : {results['delta']:.2f}°")
    print(f"  Phase Angle (φ)                  : {results['phi_deg']:.2f}°")
    print()

    print(f"CURRENTS:")
    print(f"  Armature Current (I_a)           : {results['I_a']:.4f} A")
    print(f"  d-axis Current (I_d)             : {results['I_d']:.4f} A")
    print(f"  q-axis Current (I_q)             : {results['I_q']:.4f} A")
    print()

    print(f"TORQUE:")
    print(f"  Electromagnetic Torque (T_em)    : {results['T_em']:.4f} N.m")
    print(f"  Shaft Torque (T_shaft)           : {results['T_shaft']:.4f} N.m")
    print()

    print(f"POWER:")
    print(f"  Input Power (P_in)               : {results['P_in']:.2f} W")
    print(f"  Output Power (P_out)             : {results['P_out']:.2f} W")
    print(f"  Copper Losses (P_cu)             : {results['P_cu']:.2f} W")
    print()

    print(f"PERFORMANCE:")
    print(f"  Efficiency (η)                   : {results['efficiency']:.2f} %")
    print(f"  Power Factor (cos φ)             : {results['power_factor']:.4f}")
    print()

    print("=" * 70)
    print("CALCULATION COMPLETE")
    print("=" * 70)
    print()
    print("To run the full GUI application with visualization:")
    print("  python synchronous_reluctance_motor.py")
    print()


if __name__ == "__main__":
    main()
