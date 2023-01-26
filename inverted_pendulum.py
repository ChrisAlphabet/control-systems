import control
import math
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

m = 1 
M = 5 
L = 2 
g = -9.8 
d = 1

def odes(x, t, controller):
    # constants

    position = x[0] # x, note that this is a free variable
    velocity = x[1] # x_dot
    angle = x[2] # theta
    angular_velocity = x[3] # theta_dot

    u = controller(x)

    Sx = math.sin(angle)
    Cx = math.cos(angle)
    D = m*L*L*(M+m*(1-Cx**2))

    dx = [velocity, 
         (1/D)*(-m**2*L**2*g*Cx*Sx + m*L**2*(m*L*angular_velocity**2*Sx - d*velocity)) + m*L*L*(1/D)*u,
         angular_velocity,
         (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*angular_velocity**2*Sx - d*velocity)) - m*L*Cx*(1/D)*u]

    return dx

# initial conditions
x0 = np.array([0, 0, 1, 0])

controller = lambda x: 0
t = np.linspace(0, 100, 500)
x = odeint(odes, x0, t, args=(controller,)) # note that I could pass args to my function using args=(S,) for example

position = x[:,0] # x, note that this is a free variable
velocity = x[:,1] # x_dot
angle = x[:,2] # angle
angular_velocity = x[:,3] # angle_dot

fig,ax = plt.subplots()
ax.plot(t, position,label='position')
ax.plot(t, velocity,label='velocity')
ax.plot(t, angle,label='angle')
ax.plot(t, angular_velocity,label='angular_velocity')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('x')
plt.show()


# controller
b = 1 # Pendulum up (b=1)
A = np.array([
    [0, 1, 0, 0],
    [0, -d/M, b*m*g/M, 0],
    [0, 0, 0, 1],
    [0, -b*d/(M*L), -b*(m+M)*g/(M*L), 0]
])

B = np.array([[0], [1/M], [0], [b*1/(M*L)]])

eigen_values, eigen_vectors = np.linalg.eig(A)
# print(eigen_values)

controllability = np.linalg.matrix_rank(control.ctrb(A, B))
# print(controllability) # expect 4 if controllable

# lqr controller
Q = np.identity(4)
R = 0.0001
K, S, E = control.lqr(A, B, Q, R)
print(K)

# reference position
wr = np.array([1, 0, math.pi, 0])
controller = lambda x: -K * (x - wr) # not sure if this will work
print(controller(x0))
x = odeint(odes, x0, t, args=(controller,)) # hmm, unhappy here!!

position = x[:,0] # x, note that this is a free variable
velocity = x[:,1] # x_dot
angle = x[:,2] # angle
angular_velocity = x[:,3] # angle_dot

fig,ax = plt.subplots()
ax.plot(t, position,label='position')
ax.plot(t, velocity,label='velocity')
ax.plot(t, angle,label='angle')
ax.plot(t, angular_velocity,label='angular_velocity')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('x')
plt.show()
