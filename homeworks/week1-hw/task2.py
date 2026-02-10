import matplotlib.pyplot as plt
import numpy as np
from diff_eq import DiffEq

''' pointers:
1) preallocate via np.linespace(), np.zeros
2) put DT as argument in functions for more transparent
3) error vs stepsize, use absolute value
4) plot error on log log scale
'''

DT = 0.01
T_FINAL = 2
Y_0 = 400
ACTUAL_Y = np.zeros(int(T_FINAL/DT)+1)

def f(t:float, y:float): # du/dt
   return -np.log(2) * y / 2

def actual(f, Y_0) -> None:
    t = 0
    n_steps = int(T_FINAL / DT)
    k = np.log(2) / 2  # because y' = -(ln2/2)*y
    for i in range(n_steps+1):
        ACTUAL_Y[i] = Y_0 * np.exp(-k*t)
        t += DT

actual(f, Y_0)
drug_decay = DiffEq(f, T_FINAL, DT, Y_0, ACTUAL_Y)
# drug_decay.forward_euler()
# drug_decay.plot("red", True, True, "Forward Euler")
# # print(abs(drug_decay.y_vals-drug_decay.ACTUAL_Y), "\n\n\n")
# drug_decay.explicit_midpoint()
# drug_decay.plot("green", True, True, "Explicit Midpoint")
# # print(abs(drug_decay.y_vals-drug_decay.ACTUAL_Y), "\n\n\n")
# drug_decay.rk4()
# drug_decay.plot("blue", True, True, "Low Storage Runge-Kutta 4")
# # print(abs(drug_decay.y_vals-drug_decay.ACTUAL_Y), "\n\n\n")
# plt.xlabel("Timestep", **drug_decay.labels)
# plt.ylabel("Error", **drug_decay.labels)
# plt.legend()
# plt.show()

color_list = ["red", "green", "blue"]
# time_list = [0.1, 0.01, 0.001]
time_list = [0.2, 0.1, 0.05, 0.04, 0.02, 0.01, 0.008, 0.005, 0.004, 0.001]
method_list = ["Forward Euler", "Explicit Midpoint", "Low Storage Runge-Kutta 4"]
drug_decay.plot_error(color_list, time_list, method_list)