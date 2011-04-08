% This is a basic pendulum with a torque input and extra outputs
beepsound off
overwrite all
autoz on

unitsystem kg, meter, sec

% define the frames and bodies
newtonian n

bodies a

% define the constants
% g : gravity
% l : pendulum length
% m : pendulum mass
% i : pendulum inertia
constants g, l, m, i

% states are [theta, omega, theta', omega']

variables theta'

motionvariables' omega'

% mass
mass a=m

% only has inertia about the 3 axis
inertia a,0,0,i

simprot(n,a,3,theta)

% the kinematic differential equation
theta' = omega

% the pendulum cg is at half the length
p_no_ao> = l/2 * a2>

% the kinematics
angvel(n,a)

v_ao_n>=dt(p_no_ao>, n)

alf_a_n>=dt(w_a_n>, n)

a_ao_n>=dt(v_ao_n>, n)

% the kinetics
gravity(g * n2>,a)

% input forces and torques
specified torque, force

torque = 10 * theta * sin(2 * pi * t + pi/6)
torque(n/a,torque * n3>)

force = 5 * cos(theta) + 0.5 * sin(theta')
force_ao> += force * n1>

% form the equations of motion
zero = fr() + frstar()

solve(zero, omega')

% add some extra things to output
k = ke(a)
p = -m * g * dot(p_no_ao>, n2>)
th2 = 2 * theta

% linear model
A[1,1] = d(theta', theta)
A[1,2] = d(theta', omega)
A[2,1] = d(omega', theta)
A[2,2] = d(omega', omega)

B[1,1] = d(theta', torque)
B[2,1] = d(omega', torque)
B[1,2] = d(theta', force)
B[2,2] = d(omega', force)

C[1,1] = d(theta, theta)
C[1,2] = d(theta, omega)
C[2,1] = d(omega, theta)
C[2,2] = d(omega, omega)
C[3,1] = d(k, theta)
C[3,2] = d(k, omega)
C[4,1] = d(p, theta)
C[4,2] = d(p, omega)
C[5,1] = d(th2, theta)
C[5,2] = d(th2, omega)

D[1,1] = d(theta, torque)
D[2,1] = d(omega, torque)
D[3,1] = d(k, torque)
D[4,1] = d(p, torque)
D[5,1] = d(th2, torque)
D[1,1] = d(theta, force)
D[2,1] = d(omega, force)
D[3,1] = d(k, force)
D[4,1] = d(p, force)
D[5,1] = d(th2, force)

% make sure the linear model shows up in the code
encode A,B,C,D

% give it some values for the constants
input g=9.81 meter/sec^2, l=2 meter, m=4 kg, i=0.5 kg*meter^2

% tell it what to output
output theta radian, omega radian/sec, k joules, p joules, th2 radian

code dynamics() PendulumDynamics.c

save Pendulum.all
