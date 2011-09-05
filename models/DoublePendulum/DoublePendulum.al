% This is a basic double pendulum with inputs and outputs
beepsound off
overwrite all
autoz on

unitsystem kg, meter, sec

% define the frames and bodies
newtonian n

bodies a, b

% define the constants
% g : gravity
% l1 : first pendulum length
% l2 : second pendulum length
% m1 : pendulum mass
% m2 : pendulum mass
% i1 : pendulum inertia
% i2 : pendulum inertia
constants g, l{2}, m{2}, i{2}

% states are [omega1, omega2, theta1, theta2]

variables theta{2}'

motionvariables' omega{2}'

% mass
mass a=m1, b=m2

% only has inertia about the 3 axis
inertia a, 0, 0, i1
inertia b, 0, 0, i2

simprot(n, a, 3, theta1)
simprot(a, b, 3, theta2)

% the kinematic differential equation
theta1' = omega1
theta2' = omega2

% the cg is at the middle of each rod
p_no_ao> = l1 / 2 * a2>
p_no_bo> = l1 * a2> + l2 / 2 * b2>

% the kinematics
angvel(n, a)
angvel(n, b)

v_ao_n>=dt(p_no_ao>, n)
v_bo_n>=dt(p_no_bo>, n)

alf_a_n>=dt(w_a_n>, n)
alf_b_n>=dt(w_b_n>, n)

a_ao_n>=dt(v_ao_n>, n)
a_bo_n>=dt(v_bo_n>, n)

% the kinetics
gravity(g * n2>, a, b)

% input forces and torques
% these should not depend on states
% give a time dependent expression for torque, and nothing for force
specified torque, force

torque = 10 * sin(2 * pi * t + pi/6)
torque(n/a, torque * n3>)

force_bo> += force * (b1> + b2>)

% form the equations of motion
zero = fr() + frstar()

solve(zero, omega1', omega2')

% add some extra things to output
kinetic = ke(a, b)
potential = -m1 * g * dot(p_no_ao>, n2>) - m2 * g * dot(p_no_bo>, n2>)
energy = kinetic + potential

% linear model
A[1,1] = d(omega1', omega1)
A[1,2] = d(omega1', omega2)
A[1,3] = d(omega1', theta1)
A[1,4] = d(omega1', theta2)
A[2,1] = d(omega2', omega1)
A[2,2] = d(omega2', omega2)
A[2,3] = d(omega2', theta1)
A[2,4] = d(omega2', theta2)
A[3,1] = d(theta1', omega1)
A[3,2] = d(theta1', omega2)
A[3,3] = d(theta1', theta1)
A[3,4] = d(theta1', theta2)
A[4,1] = d(theta2', omega1)
A[4,2] = d(theta2', omega2)
A[4,3] = d(theta2', theta1)
A[4,4] = d(theta2', theta2)

B[1,1] = d(omega1', torque)
B[1,2] = d(omega1', force)
B[2,1] = d(omega2', torque)
B[2,2] = d(omega2', force)
B[3,1] = d(theta1', torque)
B[3,2] = d(theta1', force)
B[4,1] = d(theta2', torque)
B[4,2] = d(theta2', force)

C[1,1] = d(omega2, omega1)
C[1,2] = d(omega2, omega2)
C[1,3] = d(omega2, theta1)
C[1,4] = d(omega2, theta2)
C[2,1] = d(theta2, omega1)
C[2,2] = d(theta2, omega2)
C[2,3] = d(theta2, theta1)
C[2,4] = d(theta2, theta2)
C[3,1] = d(kinetic, omega1)
C[3,2] = d(kinetic, omega2)
C[3,3] = d(kinetic, theta1)
C[3,4] = d(kinetic, theta2)
C[4,1] = d(potential, omega1)
C[4,2] = d(potential, omega2)
C[4,3] = d(potential, theta1)
C[4,4] = d(potential, theta2)
C[5,1] = d(energy, omega1)
C[5,2] = d(energy, omega2)
C[5,3] = d(energy, theta1)
C[5,4] = d(energy, theta2)

D[1,1] = d(omega2, torque)
D[2,1] = d(theta2, torque)
D[3,1] = d(kinetic, torque)
D[4,1] = d(potential, torque)
D[5,1] = d(energy, torque)
D[1,2] = d(omega2, force)
D[2,2] = d(theta2, force)
D[3,2] = d(kinetic, force)
D[4,2] = d(potential, force)
D[5,1] = d(energy, force)

% make sure the linear model shows up in the code
encode A,B,C,D

% give it some values for the constants
input g=9.81 meter/sec^2
input l1=2.0 meter, m1=4.0 kg, i1=0.5 kg*meter^2
input l2=2.0 meter, m2=4.0 kg, i2=0.5 kg*meter^2

% tell it what to output, these should be in order of the outputs in the C
% matrix
output omega2 radian/sec, theta2 radian, kinetic joules, potential joules, energy joules

code dynamics() DoublePendulumDynamics.c

save DoublePendulum.all
