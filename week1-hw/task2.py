import matplotlib.pyplot as plt
import numpy as np
from diff_eq import DiffEq

''' pointers:
1) preallocate via np.linespace(), np.zeros
2) put DT as argument in functions for more transparent
3) error vs stepsize, use absolute value
4) plot error on log log scale
'''

DT = 0.1
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
drug_decay.forward_euler()
drug_decay.plot("red", False, False)