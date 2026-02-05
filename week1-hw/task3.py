import matplotlib.pyplot as plt
import numpy as np

# task 3
'''
1) f(u) = [-g/l * sin theta / theta, omega]
2) plot phase portrait - plot x:theta to y:omega, get ellipse/circle unless omega very high...
'''

actual_t = [0]
actual_y = []
t_vals = []
y_vals = []
DT = 0.5

def forward_euler(f, y_0, t_final):
    n_steps = int(t_final / DT)
    u_k = y_0
    t = 0
    # y_vals.append(y_0)
    for i in range(n_steps):
        t += DT
        u_k += DT * f(t, u_k)
        t_vals.append(t)
        y_vals.append(u_k.copy() - actual_y[i])
        # y_vals.append(u_k.copy())
    return u_k


def explicit_midpoint(f, y_0, t_final):
    n_steps = int(t_final/DT)
    u_k = y_0
    t = 0
    # y_vals.append(y_0)
    for i in range(n_steps):
        t_mid = t + DT/2
        u_mid = u_k + DT / 2 * f(t_mid, u_k)
        t += DT
        u_k += DT * f(t, u_mid)
        t_vals.append(t)
        y_vals.append(u_k.copy() - actual_y[i])
        # y_vals.append(u_k.copy())
    return u_k


def rk4(f, y_0, t_final):
    n_steps = int(t_final/DT)
    u_k = y_0
    t = 0
    # y_vals.append(y_0)
    for i in range(n_steps):
        y1 = u_k + DT/4 * f(t+DT/4, u_k)
        y2 = u_k + DT/3 * f(t+DT/3, y1)
        y3 = u_k + DT/2 * f(t+DT/2, y2)
        t += DT
        u_k += DT * f(t, y3)
        t_vals.append(t)
        y_vals.append(u_k.copy() - actual_y[i])
        # y_vals.append(u_k.copy())
    return u_k


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

def actual(du, y_0, t_final):
    n_steps = int(t_final / DT)
    t = 0
    actual_y.append(y_0)

    for _ in range(n_steps):
        t += DT
        y = np.pi/4 * np.cos(t*np.sqrt(G/L))
        actual_t.append(t)
        actual_y.append(y)

    return actual_y[-1]

for i in range(4):
    t_vals = [0]
    y_vals = [[np.pi/4, 0]]
    DT = 0.5
    match i:
        case 0:
            actual(du, np.array([np.pi/4, 0]), 10)
        case 1:
            print(forward_euler(du, np.array([np.pi/4, 0]), 10))
            Y = np.array(y_vals)
            plt.plot(t_vals, Y[:,0], mfc = "red", **line_style)
            y_vals.clear()
            # plt.title("Forward Euler", **labels)
        case 2:
            print(explicit_midpoint(du, np.array([np.pi/4, 0]), 10))
            Y = np.array(y_vals)
            plt.plot(t_vals, Y[:,0], mfc = "green", **line_style)
            y_vals.clear()
            # plt.title("Explicit Midpoint", **labels)
        case 3:
            print(rk4(du, np.array([np.pi/4, 0]), 10))
            Y = np.array(y_vals)
            plt.plot(t_vals, Y[:,0], mfc = "blue", **line_style)
            y_vals.clear()
            # plt.title("Low Storage Runge-Kutta 4", **labels)
    plt.xlabel("Time (s)", **labels)
    plt.ylabel("Angle (rad)", **labels)

plt.title("All 3", **labels)
plt.show()