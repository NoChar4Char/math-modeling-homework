import numpy as np

def f_drag(prev_state, var, C_D=1.8, AREA=25*np.pi, RHO=1.2, MASS=500):
    '''
    prev_state: [x, y, dx, dy]
    var: True for x, False for y
    C_d: coefficient of drag
    area: surface area
    rho: air density
    '''

    velocity = prev_state[2] if var else prev_state[3]
    velocity_mag = np.sqrt(prev_state[2] * prev_state[2] + prev_state[3] * prev_state[3])
    return -0.5 * C_D * AREA * RHO * velocity * velocity_mag/MASS

def ds(prev_state, t, W, T_DEPLOY, G=9.81):
    '''
    prev_state: [x, y, dx, dy]
    (new_state: [dx, dy, ddx, ddy])
    t: current time
    W: random horizontal wind acceleration
    T_DEPLOY: time at which the parachute deploys and the drag force becomes nonzero
    '''
    dx = prev_state[2]
    dy = prev_state[3]
    ddx = f_drag(prev_state, True) + W if t > T_DEPLOY else W
    ddy = -G + f_drag(prev_state, False) if t > T_DEPLOY else -G
    return np.array([dx, dy, ddx, ddy])

def landing_site(H_0, VX_0 = 50, DT=0.005, T_FINAL=5000):
    '''
    NOTE THAT THIS WILL ONLY OUTPUT THE FINAL RESULT
    uses forward euler
    H_0: initial height
    W: random horizontal wind acceleration
    T_DEPLOY: time at which the parachute deploys and the drag force becomes nonzero
    VX_0: (initial) horizontal velocity
    DT: how small the step is
    T_FINAL: where the simulation should end. NOTE THAT THIS SHOULD NEVER HAPPEN
    '''
    W = np.random.normal(0, 1)
    T_DEPLOY = np.random.normal(4, 1)
    state = np.array([0, H_0, VX_0, 0])
    for t in np.arange(0, T_FINAL, DT):
        state = DT * ds(state, t, W, T_DEPLOY) + state
        if state[1] <= 0:
            break
    
    return state[0]

def get_dataset(ARRAY_LEN: int):
    landing_x_locs = np.zeros(ARRAY_LEN)
    for i in range(ARRAY_LEN):
        landing_x_locs[i] = landing_site(100)
    return landing_x_locs

np.save('landing_data', get_dataset(50000))