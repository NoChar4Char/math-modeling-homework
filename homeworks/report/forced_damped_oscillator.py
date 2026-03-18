import matplotlib.pyplot as plt
import numpy as np
from diff_eq_copy import DiffEq

PI = np.pi
DT = 0.01
T_FINAL = 100
FAKE_Y = np.zeros((int(T_FINAL/DT)+1, 2)) 
Y_0 = np.array([50, 0]) 
K = 100
B = 5
M = 200

C = np.array([[0, 1], [-K/M, -B/M]])

def f(t, u):
    return C@u

forced_damped_oscillator = DiffEq("Forced Damped Oscillator", f, T_FINAL, DT, Y_0, FAKE_Y)
forced_damped_oscillator.forward_euler()
forced_damped_oscillator.plot("red", False, False, "Forward Euler")
forced_damped_oscillator.backward_euler_oscillator(C)
forced_damped_oscillator.plot("green", False, False, "Backward Euler")
forced_damped_oscillator.explicit_midpoint()
forced_damped_oscillator.plot("green", False, False, "Explicit Midpoint")
forced_damped_oscillator.rk4()
forced_damped_oscillator.plot("green", False, False, "Low Storage Runge-Kutta 4")

plt.title("Damped Oscillator Plotted with all 4 Methods")
plt.xlabel("Time")
plt.ylabel("Distance from starting point")
plt.legend()
plt.show()