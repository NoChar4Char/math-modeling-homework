import matplotlib.pyplot as plt
import numpy as np
from diff_eq_copy import DiffEq

PI = np.pi
DT = 0.01
T_FINAL = 10
NONLIN_ACTUAL_Y = np.zeros((int(T_FINAL/DT)+1, 2)) # for now, don't know how to use ellipical stuff
Y_0 = np.array([PI/4, 0]) # default: np.array([np.pi/4, 0]), which matches the actual for linearized.
G = 10

