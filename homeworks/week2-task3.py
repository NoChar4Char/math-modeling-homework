import matplotlib.pyplot as plt
import numpy as np
from diff_eq import DiffEq

DT = 0.01
T_FINAL = 10
LIN_ACTUAL_Y = np.zeros((int(T_FINAL/DT)+1, 2))
NONLIN_ACTUAL_Y = np.zeros((int(T_FINAL/DT)+1, 2)) # for now, don't know how to use ellipical stuff
Y_0 = np.array([np.pi/4, 0])

def lin_g(t, u):
    L = 10
    G = 9.81
    A = np.array([[0, 1], [-G / L, 0]])
    return A@u

def lin_actual(lin_g, y_0, t_final, DT):
    L = 10
    G = 9.81
    n_steps = int(t_final / DT)
    t = 0
    for i in range(n_steps+1):
        LIN_ACTUAL_Y[i][0] = np.pi/4 * np.cos(t*np.sqrt(G/L))
        LIN_ACTUAL_Y[i][1] = -np.pi/4 * np.sqrt(G/L) * np.sin(t*np.sqrt(G/L))
        t += DT
    return LIN_ACTUAL_Y[-1]

lin_actual(lin_g, np.array([np.pi/4, 0]), T_FINAL, DT)
lin_pendulum = DiffEq(lin_g, T_FINAL, DT, Y_0, LIN_ACTUAL_Y)

color_list = ["red", "green", "blue", "purple"]
# time_list = [0.1, 0.01, 0.001]
time_list = [0.5, 0.2, 0.1, 0.05, 0.04, 0.02, 0.01, 0.008, 0.005, 0.004]
method_list = ["Forward Euler", "Explicit Midpoint", "Low Storage Runge-Kutta 4", "Backward Euler"]
lin_pendulum.plot_error_with_actual(color_list, time_list, method_list)