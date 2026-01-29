import matplotlib.pyplot as plt
import numpy as np


t_vals = [0]
y_vals = []
DT = 0.5


def f(t, y): # du/dt
   return -np.log(2) * y / 2


# def actual(f, y_0, t_final, DT=0.1):




def forward_euler(f, y_0, t_final):
    n_steps = int(t_final / DT)
    u_k = y_0
    t = 0
    y_vals.append(y_0)
    for _ in range(n_steps):
        t += DT
        u_k += DT * f(t, u_k)
        t_vals.append(t)
        y_vals.append(u_k)
    return u_k


def explicit_midpoint(f, y_0, t_final):
    n_steps = int(t_final/DT)
    u_k = y_0
    t = 0
    y_vals.append(y_0)
    for _ in range(n_steps):
        t_mid = t + int(DT/2)
        u_mid = u_k + DT / 2 * f(t_mid, u_k)
        t += DT
        u_k += DT * f(t, u_mid)
        t_vals.append(t)
        y_vals.append(u_k)
    return u_k


def rk4(f, y_0, t_final):
    n_steps = int(t_final/DT)
    u_k = y_0
    t = 0
    y_vals.append(y_0)
    for _ in range(n_steps):
        y1 = u_k + DT/4 * f(t+DT/4, u_k)
        y2 = u_k + DT/3 * f(t+DT/3, y1)
        y3 = u_k + DT/2 * f(t+DT/2, y2)
        t += DT
        u_k += DT * f(t, y3)
        t_vals.append(t)
        y_vals.append(u_k)
    return u_k


using = 1
using_dict = {
    1: "Forward Euler",
    2: "Explicit Midpoint",
    3: "Low Storage Runge-Kutta 4"
}


line_style = dict(
    marker="o",
    ms = 0.1,
    lw=0.1,
)


labels = dict(
    fontsize=10,
    family="Arial",
    fontweight="bold"
)
for i in range(1,4):
    match i:
        case 1:
            forward_euler(f, 400, 2)
            plt.plot(t_vals, y_vals, mfc = "red", **line_style)
            # plt.title("Forward Euler", **labels)
        case 2:
            explicit_midpoint(f, 400, 2)
            plt.plot(t_vals, y_vals, mfc = "green", **line_style)
            # plt.title("Explicit Midpoint", **labels)
        case 3:
            rk4(f, 400, 2)
            plt.plot(t_vals, y_vals, mfc = "blue", **line_style)
            # plt.title("Low Storage Runge-Kutta 4", **labels)
    t_vals = [0]
    y_vals = []
    plt.xlabel("Time (hours)", **labels)
    plt.ylabel("Amount of NSAID (mg)", **labels)






plt.title("All 3", **labels)
plt.show()