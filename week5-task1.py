import numpy as np

NUM_SAMPLES = 10000000
rng = np.random.default_rng()
rand_x = rng.random(size=NUM_SAMPLES)
rand_y = rng.random(size=NUM_SAMPLES)

indicated = 0
for i in range(NUM_SAMPLES):
    if (rand_x[i]*rand_x[i] + rand_y[i]*rand_y[i] <= 1):
        indicated += 1

pi = 4*(indicated / NUM_SAMPLES)
print(pi)