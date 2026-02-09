import matplotlib.pyplot as plt
import numpy as np

class DiffEq:
    def __init__(self, f:function, T_FINAL:int, DT:float, Y_0:float, ACTUAL_Y):
        self.dimensions = None
        try:
            self.dimensions = ACTUAL_Y.shape[1]
        except IndexError:
            self.dimensions = 1
        self.f = f
        self.T_FINAL = T_FINAL
        self.DT = DT
        self.Y_0 = Y_0
        self.ACTUAL_Y = ACTUAL_Y
        self.t_vals = np.linspace(0, T_FINAL, int(T_FINAL/DT)+1)
        self.y_vals = np.zeros((int(T_FINAL/DT)+1)) if self.dimensions == 1 else np.zeros((int(T_FINAL/DT)+1, self.dimensions))
        # self.y_vals = np.zeros((int(T_FINAL/DT)+1, self.dimensions))
        self.is_filled = False
        self.labels = dict(
            fontsize=10,
            family="Arial",
            fontweight="bold"
        )
        self.line_style = dict(
            marker="o",
            ms = 0.1,
            lw=0.5,
        )

    def forward_euler(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], self.y_vals[i])
        self.is_filled = True
        return self.y_vals[-1]

    def explicit_midpoint(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            y_mid = self.y_vals[i]+self.DT/2 * self.f(self.t_vals[i]+self.DT/2, self.y_vals[i])
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], y_mid)
        self.is_filled = True
        return self.y_vals[-1]

    def rk4(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            y1 = self.y_vals[i]+self.DT/4 * self.f(self.t_vals[i]+self.DT/4, self.y_vals[i])
            y2 = self.y_vals[i]+self.DT/3 * self.f(self.t_vals[i]+self.DT/3, y1)
            y3 = self.y_vals[i]+self.DT/2 * self.f(self.t_vals[i]+self.DT/2, y2)
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], y3)
        self.is_filled = True
        return self.y_vals[-1]
    
    def plot(self, color:str, is_loglog:bool, is_diff:bool, label_val:str) -> None:
        if (not self.is_filled):
            raise ValueError("y_vals is empty")
        y_arr = self.y_vals
        y_actual = self.ACTUAL_Y
        if (self.dimensions > 1):
            y_arr = y_arr[:,0]
            y_actual = y_actual[:,0]

        if(is_diff): 
            y_arr = np.abs(y_arr - y_actual)

        if (is_loglog):
            plt.loglog(self.t_vals, y_arr, c=color, label=label_val, **self.line_style)
        else:
            plt.plot(self.t_vals, y_arr, c=color, label=label_val, **self.line_style)
        print(y_arr)
        print("\n\n\n")
        # print(self.y_actual)
        # print(y_arr, "\n\n")
    
    def plot_actual(self, color:str, is_loglog:bool) -> None:
        if (is_loglog):
            plt.loglog(self.t_vals, self.ACTUAL_Y, mfc=color, **self.line_style)
        else:
            plt.plot(self.t_vals, self.ACTUAL_Y, mfc=color, **self.line_style)
        plt.show()

    def phase_diagram(self, color:str) -> None:
        if (self.dimensions <= 1 or not self.is_filled):
            raise ValueError("Bad conditions")
        # only works for theta and omega for now

        plt.plot(self.y_vals[:, 0], self.y_vals[:, 1], c=color, **self.line_style)
        plt.show()