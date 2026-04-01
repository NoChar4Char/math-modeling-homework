import numpy as np
import matplotlib.pyplot as plt

NUM_SAMPLES = 4000
rng = np.random.default_rng()
rand_x = rng.random(size=NUM_SAMPLES)
rand_y = rng.random(size=NUM_SAMPLES)
samples = np.linspace(1, NUM_SAMPLES, NUM_SAMPLES)

x_good_ls = []
y_good_ls = []
x_bad_ls = []
y_bad_ls = []
all_pi_estimates = np.zeros(NUM_SAMPLES)

indicated = 0
for i in range(NUM_SAMPLES):
    if (rand_x[i]*rand_x[i] + rand_y[i]*rand_y[i] <= 1):
        indicated += 1
        x_good_ls.append(rand_x[i])
        y_good_ls.append(rand_y[i])
    else:
        x_bad_ls.append(rand_x[i])
        y_bad_ls.append(rand_y[i])
    all_pi_estimates[i] = 4*(indicated / (i+1))
# 
pi = 4*(indicated / NUM_SAMPLES)
# all_pis = np.linspace(np.pi, np.pi, NUM_SAMPLES)

# x_vals = np.linspace(0, 1, 20000)
# y_vals = np.sqrt(1-x_vals**2)
# plt.plot(x_vals, y_vals, c="black")
# plt.scatter(x_good_ls, y_good_ls, facecolors='none', edgecolors='blue')
# plt.scatter(x_bad_ls, y_bad_ls, facecolors='none', edgecolors='red')
# plt.plot(samples, all_pis, color="red")
# plt.plot(samples, all_pi_estimates, color="blue")

NUM_TIMES = 300
colors = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
    "navy",
    "teal",
    "gold",
    "darkgreen",
    "crimson",
    "indigo",
    "coral",
    "turquoise",
    "magenta",
    "lime",
    "maroon",
    "slateblue",
    "darkorange",
    "forestgreen",
    "deeppink",
    "dodgerblue",
    "chocolate",
    "mediumseagreen",
    "orchid",
    "black",
]
reset_pi_estimates = np.zeros(NUM_TIMES)
for time in range(NUM_TIMES):
    indicated = 0
    # reset_pi_estimates.fill(0)
    for i in range(NUM_SAMPLES):
        rand_x = rng.random(size=NUM_SAMPLES)
        rand_y = rng.random(size=NUM_SAMPLES)
        if (rand_x[i]*rand_x[i] + rand_y[i]*rand_y[i] <= 1):
            indicated += 1
    reset_pi_estimates[time] = 4*(indicated / NUM_SAMPLES)
    # plt.plot(samples, reset_pi_estimates, c=colors[time])
    # print(reset_pi_estimates)

print(reset_pi_estimates)
# bar_types = np.linspace(3.1, 3.18, 20)
bar_types = np.linspace(2.9, 3.38, 20)
# bar_types = np.linspace(1, 5, 20)
plt.hist(reset_pi_estimates, bins=bar_types, color="blue")
plt.title(f'N = {NUM_SAMPLES}, Mean = {np.mean(reset_pi_estimates)}, Stddev = {np.std(reset_pi_estimates)}')
plt.show()