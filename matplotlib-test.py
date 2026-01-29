import matplotlib.pyplot as plt
import numpy as np


# print(plt.__version__)
x = np.array([2020, 2021, 2023, 2024, 2025, 2026])
y1 = np.array([450, 455, 469, 470, 488, 505])
y2 = np.array([455, 460, 474, 475, 500, 530])
y3 = np.array([400, 405, 410, 415, 425, 435])


line_style = dict(
    marker="o",
    ms = 5,
    mfc = "#a0ff00",
    mec = "#0fff00",
    ls="dotted",
    lw=4,
    c="#1b41d9",
)


# plt.plot(x, y, marker="o", ms = 5, mfc = "#ffff00", mec = "#ffff00", ls="dotted", lw=4, c="#1b41d9")
# markersize = ms, markerfacecolor = mfc, markeredgecolor = mec, line style = ls, linewidth = lw, c = color
plt.title("Class size", fontsize=20, family="Arial", fontweight="bold")
plt.xlabel("test x")
plt.ylabel("test y")


plt.plot(x, y1, **line_style)
# plt.plot(x, y2, **line_style)
# plt.plot(x, y3, **line_style)


plt.xticks()


plt.show()