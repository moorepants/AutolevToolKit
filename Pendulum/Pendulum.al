% This is a basic pendulum with a torque input and extra outputs
beepsound off
overwrite all
autoz on

newtonian n

bodies a

% g : gravity
% l : pendulum length
% m : pendulum mass
% i : pendulum inertia
constants g, l, m, i

% states are [theta, omega]

variables theta'

motionvariables' omega'

mass a=m

% only has inertia about the 3 axis
inertia a,0,0,i

simprot(n,a,3,theta)

% the pendulum cg is at half the length
p_no_ao> = l/2*a2>

theta' = omega

angvel(n,a)

v_ao_n>=dt(p_no_ao>,n)

alf_a_n>=dt(w_a_n>,n)

a_ao_n>=dt(v_ao_n>,n)

gravity(g*n2>,a)

specified torque, force

torque = 10*theta*sin(2*pi*t + pi/6)
force = 5*cos(theta) + 0.5*sin(theta')

torque(n/a,torque*n3>)

force_ao> = force*n1>

zero=fr()+frstar()

solve(zero, omega')

% kinetic and potential energy
k = ke(a)
p = -m*g*dot(p_no_ao>,n2>)
% a random output
th2 = 2*theta


A[1,1] = d(theta', theta)
A[1,2] = d(theta', omega)
A[2,1] = d(omega', theta)
A[2,2] = d(omega', omega)

B[1,1] = d(theta', torque)
B[2,1] = d(omega', torque)

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

input g=9.81, l=2, m=4, i=0.5
encode A,B,C,D
output t,theta,omega,k,p,th2

code dynamics() PendulumDynamics.c

save Pendulum.all
