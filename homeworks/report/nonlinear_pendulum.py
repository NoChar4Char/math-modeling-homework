import matplotlib.pyplot as plt
import numpy as np
from diff_eq_copy import DiffEq
from scipy.integrate import solve_ivp


PI = np.pi
DT = 0.01
T_FINAL = 10
NONLIN_ACTUAL_Y = np.zeros((int(T_FINAL/DT)+1, 2)) # for now, don't know how to use ellipical stuff
Y_0 = np.array([PI/4, 0]) # default: np.array([np.pi/4, 0]), which matches the actual for linearized.
G = 9.81
L = 10

def get_actual_nonlinear_trajectory(T_FINAL, DT, L=10, G=9.81):
    t_eval_points = np.arange(0, T_FINAL + DT, DT)
    
    def pendulum_deriv(t, y):
        theta, omega = y
        return [omega, -(G / L) * np.sin(theta)]

    sol = solve_ivp(pendulum_deriv, [0, T_FINAL], [np.pi / 4, 0.0],  method='RK45', t_eval=t_eval_points, rtol=1e-12, atol=1e-12)

    NONLIN_ACTUAL_Y = sol.y.T
    return NONLIN_ACTUAL_Y

def nonlin_g(t, u):
    L = 10
    G = 9.81
    R = 10
    A = np.array([u[1], -G * np.sin(u[0]) / R])
    return A

NONLIN_ACTUAL_Y = get_actual_nonlinear_trajectory(T_FINAL, DT)
nonlin_pendulum = DiffEq("Nonlinear Pendulum", nonlin_g, T_FINAL, DT, Y_0, NONLIN_ACTUAL_Y)
'''
nonlin_pendulum.forward_euler()
nonlin_pendulum.plot("red", False, False, "Forward Euler")
nonlin_pendulum.explicit_midpoint()
nonlin_pendulum.plot("green", False, False, "Explicit Midpoint")
nonlin_pendulum.rk4()
nonlin_pendulum.plot("blue", False, False, "Low Storage Runge-Kutta 4")
nonlin_pendulum.backward_euler_nonlin_pendulum()
nonlin_pendulum.plot("purple", False, False, "Backward Euler")
nonlin_pendulum.plot_actual("black", False)
'''

color_list = ["red", "green", "blue", "purple"]
# time_list = [0.1, 0.01, 0.001]
time_list = [0.5, 0.25, 0.2 ,0.1, 0.05, 0.04, 0.02, 0.025, 0.01]
method_list = ["Forward Euler", "Explicit Midpoint", "Low Storage Runge-Kutta 4", "Backward Euler"]
nonlin_pendulum.plot_error_with_actual(color_list, time_list, method_list)

'''plt.title("Nonlinear Pendulum Plotted with all 4 Methods")
plt.xlabel("Time")
plt.ylabel("Distance from starting point")
plt.legend()
plt.show()'''