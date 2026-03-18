import matplotlib.pyplot as plt
import numpy as np
from diff_eq_copy import DiffEq

PI = np.pi
DT = 0.01
T_FINAL = 4
FAKE_Y = np.zeros((int(T_FINAL/DT)+1, 2)) # for now, don't know how to use ellipical stuff
Y_0 = np.array([0, 0]) # default: np.array([np.pi/4, 0]), which matches the actual for linearized.
G = 10

B = 1.25
M = 20

def f(t, u):
    arr = np.array([u[1], G - B*u[1]*u[1]/M])
    return arr

parachute_with_drag = DiffEq("Parachute with Drag", f, T_FINAL, DT, Y_0, FAKE_Y)
parachute_with_drag.forward_euler()
parachute_with_drag.plot("red", False, False, "Forward Euler")
parachute_with_drag.backward_euler_parachute(B, M, G)
parachute_with_drag.plot("green", False, False, "Backward Euler")
parachute_with_drag.explicit_midpoint()
parachute_with_drag.plot("blue", False, False, "Explicit Midpoint")
parachute_with_drag.rk4()
parachute_with_drag.plot("purple", False, False, "Low Storage Runge-Kutta 4")

plt.title("Parachute With Drag Plotted with all 4 Methods")
plt.xlabel("Time")
plt.ylabel("Distance from starting point")
plt.legend()
plt.show()