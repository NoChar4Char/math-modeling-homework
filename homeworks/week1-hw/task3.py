import matplotlib.pyplot as plt
import numpy as np
from diff_eq import DiffEq


# task 3
'''
1) f(u) = [-g/l * sin theta / theta, omega]
2) plot phase portrait - plot x:theta to y:omega, get ellipse/circle unless omega very high...
'''

DT = 0.01
T_FINAL = 80
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

def nonlin_g(t, u):
    L = 10
    G = 9.81
    R = 10
    A = np.array([u[1], -G * np.sin(u[0]) / R])
    return A


lin_actual(lin_g, np.array([np.pi/4, 0]), T_FINAL, DT)

lin_pendulum = DiffEq(lin_g, T_FINAL, DT, Y_0, LIN_ACTUAL_Y)
'''
lin_pendulum.rk4()
lin_pendulum.plot("red", False, False, "Linear pendulum")
'''

nonlin_pendulum = DiffEq(nonlin_g, T_FINAL, DT, Y_0, lin_pendulum.y_vals)
'''
nonlin_pendulum.rk4()
nonlin_pendulum.plot("blue", False, False, "Nonlinear pendulum")
plt.legend()
plt.show()
'''

'''
color_list = ["red", "green", "blue"]
time_list = [0.1, 0.01, 0.001]
time_list = [1, 0.5, 0.2, 0.1, 0.05, 0.04, 0.02, 0.01]
method_list = ["Forward Euler", "Explicit Midpoint", "Low Storage Runge-Kutta 4"]
nonlin_pendulum.plot_error_with_other(lin_pendulum, color_list, time_list, method_list)
'''

method = "Forward Euler"
nonlin_pendulum.phase_diagram("green", method)