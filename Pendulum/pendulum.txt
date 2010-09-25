[States]
# The state order and the default initial conditions
omega = 0.0
theta = 0.0
# This gets mapped to an ordered dictionary of the states, description and
# intial conditions
# x = {}
# x['0'] = ['omega', 0.0]
# x['1'] = ['theta', 0.0]

[Parameters]
m = 1.0
g = 9.81
l = 1.0

[Equations of Motion]
z[1] = sin(theta)
omegap = -g/l*z[1] + torque/(m*l*l)
thetap = omega
# this gets mapped to the following in the integration function
# f[0] = -g/l*sin(theta)
# f[1] = omega

[Model Inputs]
torque = 10*sin(2*pi*t + pi/6)

[Model Outputs]
ke = 1/2*m*(l*thetap)**2
