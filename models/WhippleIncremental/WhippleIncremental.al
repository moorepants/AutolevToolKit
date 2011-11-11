%---------------------------------------------------------------------%
% File: WhippleIncremental.al
% Creation Date: November 11, 2011
% Author: Jason Moore
% Description: Generates the nonlinear and linear equations of motion for the
% Whipple bicycle model with three inputs.
%---------------------------------------------------------------------%
%         Default Settings
%---------------------------------------------------------------------%

autoz on
autorhs off
overwrite all
beepsound off

%---------------------------------------------------------------------%
%         newtonian, bodies, frames, particles, points
%---------------------------------------------------------------------%

% declare the inertial reference frame

newtonian n

% declare two intermediate frames
% a: yaw frame
% b: roll frame

frames a,b

% declare the four bodies
% c: bicycle frame
% d: rear wheel
% e: fork/handlebar
% f: front wheel

bodies c,d,e,f

% declare four points
% nd: rear contact point on ground
% dn: rear contact point on wheel
% nf: front contact point on ground
% fn: front contact point on wheel

points nd,dn,nf,fn

points ce

%---------------------------------------------------------------------%
%         constants and variables
%---------------------------------------------------------------------%

constants rF, rR, d{1:3}, l{1:4}

% gravity

constants g

% masses

constants mc, md, me, mf

% inertia

constants ic11, ic22, ic33, ic31
constants id11, id22
constants ie11, ie22, ie33, ie31
constants if11, if22

% input torques
% T4: roll torque
% T6: rear wheel torque
% T7: steer torque

specified T4, T6, T7

%---------------------------------------------------------------------%
% declare the generalized coordinates
%---------------------------------------------------------------------%

% q1:  perpendicular distance from the n2> axis to the rear contact
%      point in the ground plane
% q2:  perpendicular distance from the n1> axis to the rear contact
%      point in the ground plane
% q3:  frame yaw angle
% q4:  frame roll angle
% q5:  frame pitch angle
% q6:  rear wheel rotation angle
% q7:  steering rotation angle
% q8:  front wheel rotation angle

variables q{8}'

%---------------------------------------------------------------------%
%         generalized speeds
%---------------------------------------------------------------------%

motionvariables' u{8}'

%---------------------------------------------------------------------%
%         mass and inertia properties
%---------------------------------------------------------------------%

mass c=mc,d=md,e=me,f=mf
inertia c,ic11,ic22,ic33,ic12,ic23,ic31
inertia d,id11,id22,id33
inertia e,ie11,ie22,ie33,ie12,ie23,ie31
inertia f,if11,if22,if33

%---------------------------------------------------------------------%
%         angular relationships                                       %
%---------------------------------------------------------------------%

% frame yaw
simprot(n,a,3,q3)

% frame roll
simprot(a,b,1,q4)

% frame pitch
simprot(b,c,2,q5)

% rear wheel rotation
simprot(c,d,2,q6)

% steering angle
simprot(c,e,3,q7)

% front wheel rotation
simprot(e,f,2,q8)

%---------------------------------------------------------------------%
%         position vectors
%---------------------------------------------------------------------%

% locate the center of mass for each body
p_no_do>=q1*n1>+q2*n2>-rR*b3> % newtonian origin to rear wheel center
p_do_co>=l1*c1>+l2*c3> % rear wheel center to bicycle frame center
p_do_ce> = d1 * c1>
p_ce_fo>= d2*e3> + d3*e1> % rear wheel center to the front wheel center
p_fo_eo>=l3*e1>+l4*e3> % rear wheel center to fork/handlebar center

% locate the ground contact points
p_do_dn>=rr*b3>
p_dn_nd>=0>
p_fo_fn>=rF*unitvec(n3>-dot(e2>,n3>)*e2>)
p_fn_nf>=0>

%---------------------------------------------------------------------%
%         define the kinematical differential equations
%---------------------------------------------------------------------%

q1'=u1
q2'=u2
q3'=u3
q4'=u4
q5'=u5
q6'=u6
q7'=u7
q8'=u8

%---------------------------------------------------------------------%
%         angular velocities
%---------------------------------------------------------------------%

angvel(n,a)
angvel(n,b)
angvel(n,c)
angvel(n,d)
angvel(n,e)
angvel(n,f)

%---------------------------------------------------------------------%
%         velocities
%---------------------------------------------------------------------%

v_do_n> = dt(p_no_do>, n)
v2pts(n, c, do, co)
v2pts(n, c, do, ce)
v2pts(n, e, ce, fo)
v2pts(n, e, fo, eo)

v2pts(n,d,do,dn)
v2pts(n,f,fo,fn)

%---------------------------------------------------------------------%
%         define the pitch configuration constraint
%---------------------------------------------------------------------%

% set the n3> component of p_nd_nf> equal to zero
pzero=dot(p_nd_nf>,n3>)

%---------------------------------------------------------------------%
%         motion constraints
%---------------------------------------------------------------------%

% due to the assumptions of no side slip and no slip rolling the
% velocities of the front and rear wheel contact points, cn and gn,
% cannot have components of velocity in the ground plane

dependent[1]=dot(v_dn_n>,a1>)
dependent[2]=dot(v_dn_n>,a2>)
dependent[3]=dot(v_fn_n>,a1>)
dependent[4]=dot(v_fn_n>,a2>)
dependent[5]=dt(pzero)

% the rear wheel angular speed, u6, the roll rate, u4,the steering rate, u7,
% are taken to be the independent generalized speeds

constrain(dependent[u1,u2,u3,u5,u8])

%---------------------------------------------------------------------%
%         angular accelerations
%---------------------------------------------------------------------%

alf_c_n>=dt(w_c_n>,n)
alf_d_n>=dt(w_d_n>,n)
alf_e_n>=dt(w_e_n>,n)
alf_f_n>=dt(w_f_n>,n)

%---------------------------------------------------------------------%
%         accelerations
%---------------------------------------------------------------------%

a_co_n>=dt(v_co_n>,n)
a_do_n>=dt(v_do_n>,n)
a_eo_n>=dt(v_eo_n>,n)
a_fo_n>=dt(v_fo_n>,n)

%---------------------------------------------------------------------%
%         forces and torques
%---------------------------------------------------------------------%

gravity(g*n3>,c,d,e,f)
torque(a/b,T4*a1>) % roll torque
torque(c/d,T6*c2>) % rear wheel torque
torque(c/e,T7*c3>) % steer torque

%---------------------------------------------------------------------%
%         equations of motion
%---------------------------------------------------------------------%

zero=fr()+frstar()
solve(zero,u4',u6',u7')

%---------------------------------------------------------------------%
%         linearization
%---------------------------------------------------------------------%
% linearizes the equations of motion

aMat[1,1]=d(q1',q1)
aMat[1,2]=d(q1',q2)
aMat[1,3]=d(q1',q3)
aMat[1,4]=d(q1',q4)
aMat[1,5]=d(q1',q5)
aMat[1,6]=d(q1',q6)
aMat[1,7]=d(q1',q7)
aMat[1,8]=d(q1',q8)
aMat[1,9]=d(q1',u4)
aMat[1,10]=d(q1',u6)
aMat[1,11]=d(q1',u7)

aMat[2,1]=d(q2',q1)
aMat[2,2]=d(q2',q2)
aMat[2,3]=d(q2',q3)
aMat[2,4]=d(q2',q4)
aMat[2,5]=d(q2',q5)
aMat[2,6]=d(q2',q6)
aMat[2,7]=d(q2',q7)
aMat[2,8]=d(q2',q8)
aMat[2,9]=d(q2',u4)
aMat[2,10]=d(q2',u6)
aMat[2,11]=d(q2',u7)

aMat[3,1]=d(q3',q1)
aMat[3,2]=d(q3',q2)
aMat[3,3]=d(q3',q3)
aMat[3,4]=d(q3',q4)
aMat[3,5]=d(q3',q5)
aMat[3,6]=d(q3',q6)
aMat[3,7]=d(q3',q7)
aMat[3,8]=d(q3',q8)
aMat[3,9]=d(q3',u4)
aMat[3,10]=d(q3',u6)
aMat[3,11]=d(q3',u7)

aMat[4,1]=d(q4',q1)
aMat[4,2]=d(q4',q2)
aMat[4,3]=d(q4',q3)
aMat[4,4]=d(q4',q4)
aMat[4,5]=d(q4',q5)
aMat[4,6]=d(q4',q6)
aMat[4,7]=d(q4',q7)
aMat[4,8]=d(q4',q8)
aMat[4,9]=d(q4',u4)
aMat[4,10]=d(q4',u6)
aMat[4,11]=d(q4',u7)

aMat[5,1]=d(q5',q1)
aMat[5,2]=d(q5',q2)
aMat[5,3]=d(q5',q3)
aMat[5,4]=d(q5',q4)
aMat[5,5]=d(q5',q5)
aMat[5,6]=d(q5',q6)
aMat[5,7]=d(q5',q7)
aMat[5,8]=d(q5',q8)
aMat[5,9]=d(q5',u4)
aMat[5,10]=d(q5',u6)
aMat[5,11]=d(q5',u7)

aMat[6,1]=d(q6',q1)
aMat[6,2]=d(q6',q2)
aMat[6,3]=d(q6',q3)
aMat[6,4]=d(q6',q4)
aMat[6,5]=d(q6',q5)
aMat[6,6]=d(q6',q6)
aMat[6,7]=d(q6',q7)
aMat[6,8]=d(q6',q8)
aMat[6,9]=d(q6',u4)
aMat[6,10]=d(q6',u6)
aMat[6,11]=d(q6',u7)

aMat[7,1]=d(q7',q1)
aMat[7,2]=d(q7',q2)
aMat[7,3]=d(q7',q3)
aMat[7,4]=d(q7',q4)
aMat[7,5]=d(q7',q5)
aMat[7,6]=d(q7',q6)
aMat[7,7]=d(q7',q7)
aMat[7,8]=d(q7',q8)
aMat[7,9]=d(q7',u4)
aMat[7,10]=d(q7',u6)
aMat[7,11]=d(q7',u7)

aMat[8,1]=d(q8',q1)
aMat[8,2]=d(q8',q2)
aMat[8,3]=d(q8',q3)
aMat[8,4]=d(q8',q4)
aMat[8,5]=d(q8',q5)
aMat[8,6]=d(q8',q6)
aMat[8,7]=d(q8',q7)
aMat[8,8]=d(q8',q8)
aMat[8,9]=d(q8',u4)
aMat[8,10]=d(q8',u6)
aMat[8,11]=d(q8',u7)

aMat[9,1]=d(u4',q1)
aMat[9,2]=d(u4',q2)
aMat[9,3]=d(u4',q3)
aMat[9,4]=d(u4',q4)
aMat[9,5]=d(u4',q5)
aMat[9,6]=d(u4',q6)
aMat[9,7]=d(u4',q7)
aMat[9,8]=d(u4',q8)
aMat[9,9]=d(u4',u4)
aMat[9,10]=d(u4',u6)
aMat[9,11]=d(u4',u7)

aMat[10,1]=d(u6',q1)
aMat[10,2]=d(u6',q2)
aMat[10,3]=d(u6',q3)
aMat[10,4]=d(u6',q4)
aMat[10,5]=d(u6',q5)
aMat[10,6]=d(u6',q6)
aMat[10,7]=d(u6',q7)
aMat[10,8]=d(u6',q8)
aMat[10,9]=d(u6',u4)
aMat[10,10]=d(u6',u6)
aMat[10,11]=d(u6',u7)

aMat[11,1]=d(u7',q1)
aMat[11,2]=d(u7',q2)
aMat[11,3]=d(u7',q3)
aMat[11,4]=d(u7',q4)
aMat[11,5]=d(u7',q5)
aMat[11,6]=d(u7',q6)
aMat[11,7]=d(u7',q7)
aMat[11,8]=d(u7',q8)
aMat[11,9]=d(u7',u4)
aMat[11,10]=d(u7',u6)
aMat[11,11]=d(u7',u7)

bMat[1,1]=d(q1',T4)
bMat[1,2]=d(q1',T7)

bMat[2,1]=d(q2',T4)
bMat[2,2]=d(q2',T7)

bMat[3,1]=d(q3',T4)
bMat[3,2]=d(q3',T7)

bMat[4,1]=d(q4',T4)
bMat[4,2]=d(q4',T7)

bMat[5,1]=d(q5',T4)
bMat[5,2]=d(q5',T7)

bMat[6,1]=d(q6',T4)
bMat[6,2]=d(q6',T7)

bMat[7,1]=d(q7',T4)
bMat[7,2]=d(q7',T7)

bMat[8,1]=d(q8',T4)
bMat[8,2]=d(q8',T7)

bMat[9,1]=d(u4',T4)
bMat[9,2]=d(u4',T7)

bMat[10,1]=d(u6',T4)
bMat[10,2]=d(u6',T7)

bMat[11,1]=d(u7',T4)
bMat[11,2]=d(u7',T7)

cMat[1,1]=d(q1,q1)
cMat[1,2]=d(q1,q2)
cMat[1,3]=d(q1,q3)
cMat[1,4]=d(q1,q4)
cMat[1,5]=d(q1,q5)
cMat[1,6]=d(q1,q6)
cMat[1,7]=d(q1,q7)
cMat[1,8]=d(q1,q8)
cMat[1,9]=d(q1,u4)
cMat[1,10]=d(q1,u6)
cMat[1,11]=d(q1,u7)

cMat[2,1]=d(q2,q1)
cMat[2,2]=d(q2,q2)
cMat[2,3]=d(q2,q3)
cMat[2,4]=d(q2,q4)
cMat[2,5]=d(q2,q5)
cMat[2,6]=d(q2,q6)
cMat[2,7]=d(q2,q7)
cMat[2,8]=d(q2,q8)
cMat[2,9]=d(q2,u4)
cMat[2,10]=d(q2,u6)
cMat[2,11]=d(q2,u7)

cMat[3,1]=d(q3,q1)
cMat[3,2]=d(q3,q2)
cMat[3,3]=d(q3,q3)
cMat[3,4]=d(q3,q4)
cMat[3,5]=d(q3,q5)
cMat[3,6]=d(q3,q6)
cMat[3,7]=d(q3,q7)
cMat[3,8]=d(q3,q8)
cMat[3,9]=d(q3,u4)
cMat[3,10]=d(q3,u6)
cMat[3,11]=d(q3,u7)

cMat[4,1]=d(q4,q1)
cMat[4,2]=d(q4,q2)
cMat[4,3]=d(q4,q3)
cMat[4,4]=d(q4,q4)
cMat[4,5]=d(q4,q5)
cMat[4,6]=d(q4,q6)
cMat[4,7]=d(q4,q7)
cMat[4,8]=d(q4,q8)
cMat[4,9]=d(q4,u4)
cMat[4,10]=d(q4,u6)
cMat[4,11]=d(q4,u7)

cMat[5,1]=d(q5,q1)
cMat[5,2]=d(q5,q2)
cMat[5,3]=d(q5,q3)
cMat[5,4]=d(q5,q4)
cMat[5,5]=d(q5,q5)
cMat[5,6]=d(q5,q6)
cMat[5,7]=d(q5,q7)
cMat[5,8]=d(q5,q8)
cMat[5,9]=d(q5,u4)
cMat[5,10]=d(q5,u6)
cMat[5,11]=d(q5,u7)

cMat[6,1]=d(q6,q1)
cMat[6,2]=d(q6,q2)
cMat[6,3]=d(q6,q3)
cMat[6,4]=d(q6,q4)
cMat[6,5]=d(q6,q5)
cMat[6,6]=d(q6,q6)
cMat[6,7]=d(q6,q7)
cMat[6,8]=d(q6,q8)
cMat[6,9]=d(q6,u4)
cMat[6,10]=d(q6,u6)
cMat[6,11]=d(q6,u7)

cMat[7,1]=d(q7,q1)
cMat[7,2]=d(q7,q2)
cMat[7,3]=d(q7,q3)
cMat[7,4]=d(q7,q4)
cMat[7,5]=d(q7,q5)
cMat[7,6]=d(q7,q6)
cMat[7,7]=d(q7,q7)
cMat[7,8]=d(q7,q8)
cMat[7,9]=d(q7,u4)
cMat[7,10]=d(q7,u6)
cMat[7,11]=d(q7,u7)

cMat[8,1]=d(q8,q1)
cMat[8,2]=d(q8,q2)
cMat[8,3]=d(q8,q3)
cMat[8,4]=d(q8,q4)
cMat[8,5]=d(q8,q5)
cMat[8,6]=d(q8,q6)
cMat[8,7]=d(q8,q7)
cMat[8,8]=d(q8,q8)
cMat[8,9]=d(q8,u4)
cMat[8,10]=d(q8,u6)
cMat[8,11]=d(q8,u7)

cMat[9,1]=d(u1,q1)
cMat[9,2]=d(u1,q2)
cMat[9,3]=d(u1,q3)
cMat[9,4]=d(u1,q4)
cMat[9,5]=d(u1,q5)
cMat[9,6]=d(u1,q6)
cMat[9,7]=d(u1,q7)
cMat[9,8]=d(u1,q8)
cMat[9,9]=d(u1,u4)
cMat[9,10]=d(u1,u6)
cMat[9,11]=d(u1,u7)

cMat[10,1]=d(u2,q1)
cMat[10,2]=d(u2,q2)
cMat[10,3]=d(u2,q3)
cMat[10,4]=d(u2,q4)
cMat[10,5]=d(u2,q5)
cMat[10,6]=d(u2,q6)
cMat[10,7]=d(u2,q7)
cMat[10,8]=d(u2,q8)
cMat[10,9]=d(u2,u4)
cMat[10,10]=d(u2,u6)
cMat[10,11]=d(u2,u7)

cMat[11,1]=d(u3,q1)
cMat[11,2]=d(u3,q2)
cMat[11,3]=d(u3,q3)
cMat[11,4]=d(u3,q4)
cMat[11,5]=d(u3,q5)
cMat[11,6]=d(u3,q6)
cMat[11,7]=d(u3,q7)
cMat[11,8]=d(u3,q8)
cMat[11,9]=d(u3,u4)
cMat[11,10]=d(u3,u6)
cMat[11,11]=d(u3,u7)

cMat[12,1]=d(u4,q1)
cMat[12,2]=d(u4,q2)
cMat[12,3]=d(u4,q3)
cMat[12,4]=d(u4,q4)
cMat[12,5]=d(u4,q5)
cMat[12,6]=d(u4,q6)
cMat[12,7]=d(u4,q7)
cMat[12,8]=d(u4,q8)
cMat[12,9]=d(u4,u4)
cMat[12,10]=d(u4,u6)
cMat[12,11]=d(u4,u7)

cMat[13,1]=d(u5,q1)
cMat[13,2]=d(u5,q2)
cMat[13,3]=d(u5,q3)
cMat[13,4]=d(u5,q4)
cMat[13,5]=d(u5,q5)
cMat[13,6]=d(u5,q6)
cMat[13,7]=d(u5,q7)
cMat[13,8]=d(u5,q8)
cMat[13,9]=d(u5,u4)
cMat[13,10]=d(u5,u6)
cMat[13,11]=d(u5,u7)

cMat[14,1]=d(u6,q1)
cMat[14,2]=d(u6,q2)
cMat[14,3]=d(u6,q3)
cMat[14,4]=d(u6,q4)
cMat[14,5]=d(u6,q5)
cMat[14,6]=d(u6,q6)
cMat[14,7]=d(u6,q7)
cMat[14,8]=d(u6,q8)
cMat[14,9]=d(u6,u4)
cMat[14,10]=d(u6,u6)
cMat[14,11]=d(u6,u7)

cMat[15,1]=d(u7,q1)
cMat[15,2]=d(u7,q2)
cMat[15,3]=d(u7,q3)
cMat[15,4]=d(u7,q4)
cMat[15,5]=d(u7,q5)
cMat[15,6]=d(u7,q6)
cMat[15,7]=d(u7,q7)
cMat[15,8]=d(u7,q8)
cMat[15,9]=d(u7,u4)
cMat[15,10]=d(u7,u6)
cMat[15,11]=d(u7,u7)

cMat[16,1]=d(u8,q1)
cMat[16,2]=d(u8,q2)
cMat[16,3]=d(u8,q3)
cMat[16,4]=d(u8,q4)
cMat[16,5]=d(u8,q5)
cMat[16,6]=d(u8,q6)
cMat[16,7]=d(u8,q7)
cMat[16,8]=d(u8,q8)
cMat[16,9]=d(u8,u4)
cMat[16,10]=d(u8,u6)
cMat[16,11]=d(u8,u7)

dMat[1,1]=d(q1,T4)
dMat[1,2]=d(q1,T7)

dMat[2,1]=d(q2,T4)
dMat[2,2]=d(q2,T7)

dMat[3,1]=d(q3,T4)
dMat[3,2]=d(q3,T7)

dMat[4,1]=d(q4,T4)
dMat[4,2]=d(q4,T7)

dMat[5,1]=d(q5,T4)
dMat[5,2]=d(q5,T7)

dMat[6,1]=d(q6,T4)
dMat[6,2]=d(q6,T7)

dMat[7,1]=d(q7,T4)
dMat[7,2]=d(q7,T7)

dMat[8,1]=d(q8,T4)
dMat[8,2]=d(q8,T7)

dMat[9,1]=d(u1,T4)
dMat[9,2]=d(u1,T7)

dMat[10,1]=d(u2,T4)
dMat[10,2]=d(u2,T7)

dMat[11,1]=d(u3,T4)
dMat[11,2]=d(u3,T7)

dMat[12,1]=d(u4,T4)
dMat[12,2]=d(u4,T7)

dMat[13,1]=d(u5,T4)
dMat[13,2]=d(u5,T7)

dMat[14,1]=d(u6,T4)
dMat[14,2]=d(u6,T7)

dMat[15,1]=d(u7,T4)
dMat[15,2]=d(u7,T7)

dMat[16,1]=d(u8,T4)
dMat[16,2]=d(u8,T7)


unitsystem kg, meter, sec

input q1=0.0 m, q2=0.0 m, q3=0.0 rad, q4=0.0 rad, q5=0.0 rad, q6=0.0 rad, q7=0.0 rad, q8=0.0 rad
input u4=0.0 rad/s, u6=0.0 rad/s, u7=0.0 rad/s

output q1 m, q2 m, q3 rad, q4 rad, q5 rad, q6 rad, q7 rad, q8 rad
output u1 m/s, u2 m/s, u3 rad/s, u4 rad/s, u5 rad/s, u6 rad/s, u7 rad/s, u8 rad/s

code dynamics() WhippleIncrementalDynamics.m

encode aMat,bMat,cMat,dMat

code dynamics() WhippleIncrementalDynamics.c
code algebraic() WhippleIncrementalAlgebraic.m

%---------------------------------------------------------------------%
%         save output
%---------------------------------------------------------------------%

save WhippleIncremental.all

%---------------------------------------------------------------------%
