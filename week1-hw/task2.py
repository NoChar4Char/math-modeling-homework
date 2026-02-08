import matplotlib.pyplot as plt
import numpy as np

''' pointers:
1) preallocate via np.linespace(), np.zeros
2) put DT as argument in functions for more transparent
3) error vs stepsize, use absolute value
4) plot error on log log scale
'''

DT = 0.001
T_FINAL = 2
actual_y = np.zeros(int(T_FINAL/DT)+1)
t_vals = np.linspace(0, T_FINAL, int(T_FINAL/DT)+1)
y_vals = np.zeros(int(T_FINAL/DT)+1)

def f(t: float, y: float): # du/dt
   return -np.log(2) * y / 2


def actual(f, y_0, t_final, DT):
    n_steps = int(t_final / DT)
    k = np.log(2) / 2  # because y' = -(ln2/2)*y
    actual_y[0] = y_0
    for i in range(n_steps):
        actual_y[i+1] = y_0 * np.exp(-k * t_vals[i])
    return actual_y[-1]

def forward_euler(f, y_0, t_final, DT):
    n_steps = int(t_final/DT)
    y_vals[0] = y_0
    for i in range(n_steps):
        y_vals[i+1] = y_vals[i]+DT * f(t_vals[i+1], y_vals[i])
    return y_vals[-1]

def explicit_midpoint(f, y_0, t_final, DT):
    n_steps = int(t_final/DT)
    y_vals[0] = y_0
    for i in range(n_steps):
        y_mid = y_vals[i]+DT/2 * f(t_vals[i]+DT/2, y_vals[i])
        y_vals[i+1] = y_vals[i]+DT * f(t_vals[i+1], y_mid)
    return y_vals[-1]

def rk4(f, y_0, t_final, DT):
    n_steps = int(t_final/DT)
    y_vals[0] = y_0
    for i in range(n_steps):
        y1 = y_vals[i]+DT/4 * f(t_vals[i]+DT/4, y_vals[i])
        y2 = y_vals[i]+DT/3 * f(t_vals[i]+DT/3, y1)
        y3 = y_vals[i]+DT/2 * f(t_vals[i]+DT/2, y2)
        y_vals[i+1] = y_vals[i]+DT * f(t_vals[i+1], y3)
    return y_vals[-1]


line_style = dict(
    marker="o",
    ms = 0.1,
    lw=0.5,
)

labels = dict(
    fontsize=10,
    family="Arial",
    fontweight="bold"
)

actual(f, 400, 2, DT)
for i in range(3):
    match i:
        case 0:
            print(forward_euler(f, 400, 2, DT))
            plt.loglog(t_vals, abs(y_vals-actual_y), mfc = "red", **line_style)
        case 1:
            print(explicit_midpoint(f, 400, 2, DT))
            plt.loglog(t_vals, abs(y_vals-actual_y), mfc = "green", **line_style)
        case 2:
            print(rk4(f, 400, 2, DT))
            plt.loglog(t_vals, abs(y_vals-actual_y), mfc = "blue", **line_style)
    plt.xlabel("Time (hours)", **labels)
    plt.ylabel("Deviation (mg)", **labels)

plt.title("All 3", **labels)
plt.show()
