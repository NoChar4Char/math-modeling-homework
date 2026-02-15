import matplotlib.pyplot as plt
import numpy as np

class DiffEq:
    def __init__(self, name:str, f:function, T_FINAL:int, DT:float, Y_0:float, ACTUAL_Y): # broken if DT not divisible by T_FINAL
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
        if T_FINAL/DT != int(T_FINAL/DT):
            raise ValueError(f"DT={DT} doesn't divide into T_FINAL={T_FINAL}, remainder is {int(T_FINAL/DT)}")
        self.t_vals = np.linspace(0, T_FINAL, int(T_FINAL/DT) + 1)
        self.y_vals = np.zeros((int(T_FINAL/DT)+1)) if self.dimensions == 1 else np.zeros((int(T_FINAL/DT)+1, self.dimensions))
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
        self.name = name

    def forward_euler(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], self.y_vals[i])
        self.is_filled = True
        # return self.y_vals[-1]

    def explicit_midpoint(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            y_mid = self.y_vals[i]+self.DT/2 * self.f(self.t_vals[i]+self.DT/2, self.y_vals[i])
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], y_mid)
        self.is_filled = True
        # return self.y_vals[-1]

    def rk4(self):
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            y1 = self.y_vals[i]+self.DT/4 * self.f(self.t_vals[i]+self.DT/4, self.y_vals[i])
            y2 = self.y_vals[i]+self.DT/3 * self.f(self.t_vals[i]+self.DT/3, y1)
            y3 = self.y_vals[i]+self.DT/2 * self.f(self.t_vals[i]+self.DT/2, y2)
            self.y_vals[i+1] = self.y_vals[i]+self.DT * self.f(self.t_vals[i+1], y3)
        self.is_filled = True
        # return self.y_vals[-1]
    
    def backward_euler(self):
        # only for linear pendulum
        if self.name.lower() != "linear pendulum":
            raise NameError("Bad differential equation")
        I = np.array([[1, 0], [0, 1]])
        A = np.array([[0, 1], [-9.81/10, 0]])
        n_steps = int(self.T_FINAL/self.DT)
        self.y_vals[0] = self.Y_0
        for i in range(n_steps):
            self.y_vals[i+1] = np.linalg.solve(I - self.DT*A, self.y_vals[i])
        self.is_filled = True

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
    
    def plot_actual(self, color:str, is_loglog:bool) -> None:
        y_actual = self.ACTUAL_Y
        if (self.dimensions > 1):
            y_actual = y_actual[:,0]

        if (is_loglog):
            plt.loglog(self.t_vals, y_actual, mfc=color, label = "Actual Graph", **self.line_style)
        else:
            plt.plot(self.t_vals, y_actual, mfc=color, label = "Actual Graph", **self.line_style)

    def phase_portrait(self, color:str, method:str) -> None:
        match method.lower():
            case "forward euler" | "1":
                self.forward_euler()
            case "explicit midpoint" | "2":
                self.explicit_midpoint()
            case "low storage runge-kutta 4" | "3":
                self.rk4()
            case _:
                print("???")
                return None
        if not self.dimensions == 2:
            raise ValueError("Bad dimensions")

        plt.plot(self.y_vals[:,0], self.y_vals[:,1], c=color, label=method, **self.line_style)
        print(self.y_vals[:,0], "\n", self.y_vals[:,1])
        plt.title(f"Phase Portrait for {self.name}", **self.labels)
        plt.xlabel("Theta (Angular Displacement)", **self.labels)
        plt.ylabel("Omega (Angular Velocity)", **self.labels)
        plt.legend()
        plt.show()

    def plot_error_with_actual(self, color_list:list, time_list:list, method_list:list) -> None: # actual task1
        original_DT = self.DT
        for i in range(len(method_list)):
            points = np.zeros(len(time_list))
            for j in range(len(time_list)):
                self.DT = time_list[j]
                self.t_vals = np.linspace(0, self.T_FINAL, int(self.T_FINAL/self.DT) + 1)
                self.y_vals = np.zeros((int(self.T_FINAL/self.DT)+1)) if self.dimensions == 1 else np.zeros((int(self.T_FINAL/self.DT)+1, self.dimensions))
                if self.T_FINAL/self.DT != int(self.T_FINAL/self.DT):
                    raise ValueError(f"DT={self.DT} doesn't divide into T_FINAL={self.T_FINAL}")
                match method_list[i].lower():
                    case "forward euler" | "1":
                        self.forward_euler()
                    case "explicit midpoint" | "2":
                        self.explicit_midpoint()
                    case "low storage runge-kutta 4" | "3":
                        self.rk4()
                    case "backward euler" | "4":
                        self.backward_euler()
                    case _:
                        raise ValueError("Bad methods")
                
                points[j] = np.abs(self.y_vals[-1] - self.ACTUAL_Y[-1]) if self.dimensions == 1 else np.abs(self.y_vals[-1][0] - self.ACTUAL_Y[-1][0])
            plt.loglog(time_list, points, c=color_list[i], label=method_list[i], **self.line_style)
        plt.legend()
        plt.title(f"Plot {self.name} Error with Actual", **self.labels)
        plt.xlabel("Timestep", **self.labels)
        plt.ylabel("Error", **self.labels)
        plt.show()
        self.DT = original_DT

    def plot_error_with_other(self, other_dq:DiffEq, color_list:list, time_list:list, method_list:list) -> None: # actual task1
        original_DT = self.DT # should be equal to other_dq.DT
        for i in range(len(method_list)):
            points = np.zeros(len(time_list))
            for j in range(len(time_list)):
                self.DT = time_list[j]
                self.t_vals = np.linspace(0, self.T_FINAL, int(self.T_FINAL/self.DT) + 1)
                self.y_vals = np.zeros((int(self.T_FINAL/self.DT)+1)) if self.dimensions == 1 else np.zeros((int(self.T_FINAL/self.DT)+1, self.dimensions))
                other_dq.t_vals = np.linspace(0, other_dq.T_FINAL, int(other_dq.T_FINAL/other_dq.DT) + 1)
                other_dq.y_vals = np.zeros((int(other_dq.T_FINAL/other_dq.DT)+1)) if other_dq.dimensions == 1 else np.zeros((int(other_dq.T_FINAL/other_dq.DT)+1, other_dq.dimensions))
                if self.T_FINAL/self.DT != int(self.T_FINAL/self.DT):
                    raise ValueError(f"DT={self.DT} doesn't divide into T_FINAL={self.T_FINAL}")
                match method_list[i].lower():
                    case "forward euler" | "1":
                        self.forward_euler()
                        other_dq.forward_euler()
                    case "explicit midpoint" | "2":
                        self.explicit_midpoint()
                        other_dq.explicit_midpoint()
                    case "low storage runge-kutta 4" | "3":
                        self.rk4()
                        other_dq.rk4()
                    case _:
                        print("???")
                        continue
                
                points[j] = np.abs(self.y_vals[-1] - other_dq.y_vals[-1]) if self.dimensions == 1 else np.abs(self.y_vals[-1][0] - other_dq.y_vals[-1][0])
            plt.loglog(time_list, points, c=color_list[i], label=method_list[i], **self.line_style)
        plt.legend()
        plt.title(f"Plot {self.name} Error with Other {other_dq.name}", **self.labels)
        plt.xlabel("Timestep", **self.labels)
        plt.ylabel("Error", **self.labels)
        plt.show()
        self.DT = original_DT
        other_dq.DT = original_DT