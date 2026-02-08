import matplotlib.pyplot as plt
import numpy as np

# task 3
'''
1) f(u) = [-g/l * sin theta / theta, omega]
2) plot phase portrait - plot x:theta to y:omega, get ellipse/circle unless omega very high...
'''

DT = 0.001
T_FINAL = 10
actual_y = np.zeros((int(T_FINAL/DT)+1, 2))
t_vals = np.linspace(0, T_FINAL, int(T_FINAL/DT)+1)
y_vals = np.zeros((int(T_FINAL/DT)+1, 2))

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

L = 10
G = 9.81
A = np.array([[0, 1], [-G / L, 0]])
I = np.array([[1, 0], [0, 1]])

def du(t, u):
    return A@u

def actual(du, y_0, t_final, DT):
    n_steps = int(t_final / DT)
    actual_y[0] = y_0
    for i in range(n_steps):
        actual_y[i+1] = np.pi/4 * np.cos(t_vals[i]*np.sqrt(G/L))
    return actual_y[-1]

actual(du, np.array([np.pi/4, 0]), 10, DT)
for i in range(4):
    match 0:
        case 0:
            print(forward_euler(du, np.array([np.pi/4, 0]), 10, DT))
            plt.plot(t_vals, abs(y_vals[:,0]-actual_y[:,0]), mfc = "red", **line_style)
        case 1:
            print(explicit_midpoint(du, np.array([np.pi/4, 0]), 10, DT))
            plt.plot(t_vals, abs(y_vals[:,0]-actual_y[:,0]), mfc = "green", **line_style)
        case 2:
            print(rk4(du, np.array([np.pi/4, 0]), 10, DT))
            plt.plot(t_vals, abs(y_vals[:,0]-actual_y[:,0]), mfc = "blue", **line_style)
    plt.xlabel("Time (s)", **labels)
    plt.ylabel("Angle (rad)", **labels)

plt.title("All 3", **labels)
plt.show()