import matplotlib.pyplot as plt
import numpy as np
DT = 0.00000001
TAU = 0.0000000001

def f(x):
    try:
        return 20.2*np.exp(-157990000.0/(4157.0*(3416.0-361.0*x))) - (x*x / ((12.0-0.5*x)*(1.0-0.5*x)))
    except ZeroDivisionError:
        return -1000000000 #asymptote - either -infinity or +infinity
def g(x)->float:
    try:
        return x - (20.2*np.exp(-157990000.0/(4157.0*(3416.0-361.0*x))) - (x*x / ((12.0-0.5*x)*(1.0-0.5*x))))
    except ZeroDivisionError:
        return -1000000000 #asymptote - either -infinity or +infinity
    
def df(f:function, x:float):
    return (f(x+DT) - f(x)) / DT
# Task 2
def bisection():
    lo = -80.0
    hi = -20.0
    while hi-lo > TAU:
        mid = lo + (hi-lo)/2
        if f(mid)*f(lo)<0: # diff signs
            hi = mid
        else:
            lo = mid
    return lo

# print(bisection())

#Task 3
def bisection(xa:float, xb:float, TAU:float):
    lo = xa
    hi = xb
    if f(lo) * f(hi) > 0:
        raise ValueError("Bad values")
    while hi-lo > TAU:
        mid = lo + (hi-lo)/2
        if f(mid)*f(lo)<0: # diff signs
            hi = mid
        else:
            lo = mid
    return lo
#Task 4
truth_actual = 0.05687710815030454 #via 100 iterations
truth_tau = 0.0568762004442215 #using tau = 1e-6
def bisection_error(xa:float, xb:float, num_of_iterations:int, truth_value:float):
    lo = xa
    hi = xb
    if f(lo) * f(hi) > 0:
        raise ValueError("Bad values")
    for _ in range(num_of_iterations):
        mid = lo + (hi-lo)/2
        if f(mid)*f(lo)<0: # diff signs
            hi = mid
        else:
            lo = mid
    return abs((hi+lo)/2 - truth_value)
def plot_log_bisection_error(xa:float, xb:float, truth_value:float, num_of_iterator_values:int):
    iterators = np.arange(1, num_of_iterator_values+1)
    errors = np.zeros(num_of_iterator_values)
    labels = dict(fontsize=10, family="Arial", fontweight="bold")
    line_style = dict(marker="o", ms = 0.1, lw=0.5)
    for i in range(num_of_iterator_values):
        errors[i] = bisection_error(xa, xb, iterators[i], truth_value)
    plt.plot(iterators, errors, label="Error Plot of Biesction", **line_style)
    plt.yscale('log', base=2)
    plt.legend()
    plt.title("Error Plot of Biesction", **labels)
    plt.xlabel("Number of iterations", **labels)
    plt.ylabel("Error Value", **labels)
    plt.show()

# plot_log_bisection_error(10**-6, 1.99, truth_actual, 60)

#Task 5
def fixed_point_iteration(x0:float, TAU:float=10**-9):
    list_vals = [g(x0)]
    i = 0
    while i == 0 or abs(list_vals[i] - list_vals[i-1]) > TAU and i != 10000:
        list_vals.append(g(list_vals[i]))
        i += 1
    return list_vals

# print(fixed_point_iteration(-30))

def newton_raphson(x0:float, f:function, df:function, TAU:float=10**-9):
    x_vals = [x0]
    i = 0
    while i == 0 or abs(x_vals[i] - x_vals[i-1]) > TAU and i != 10000:
        prev_x_val = x_vals[i]
        x_vals.append(prev_x_val - f(prev_x_val) / df(f, prev_x_val))
        i += 1
    return x_vals[-1]

print(newton_raphson(-30, f, df))