import numpy as np
import matplotlib.pyplot as plt

data = np.load('landing_data.npy')
# np.set_printoptions(threshold=np.inf) 

# plt.hist(data)
# plt.show()

print(f'Raw Data. Mean: {np.mean(data)}; Standard Deviation: {np.std(data)}.')
print(data)

reshaped_data = np.reshape(data, (500, 100))

def is_danger_zone(data_col):
    is_danger = 0
    for distance in data_col:
        if distance < 190 or distance > 260:
            is_danger += 1
    return is_danger / len(data_col)

def estimate_danger_prob():
    prob_of_cols = np.zeros(len(reshaped_data))
    for i in range(len(reshaped_data)):
        column = reshaped_data[i]
        prob_of_cols[i] = is_danger_zone(column)
    np.set_printoptions(threshold=np.inf) 
    print(prob_of_cols)
    plt.hist(prob_of_cols)
    plt.show()
    print(f'Probabilities. Mean: {np.mean(prob_of_cols)}; Standard Deviation: {np.std(prob_of_cols)}.')

estimate_danger_prob()