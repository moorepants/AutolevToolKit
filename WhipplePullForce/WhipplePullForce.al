%---------------------------------------------------------------------%
% File: WhipplePullForce.al
% Creation Date: March 31, 2011
% Author: Jason Moore
% Description: Generates the nonlinear and linear equations of motion for the
% Whipple bicycle model.
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

% declare six bodies
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

%---------------------------------------------------------------------%
%         constants and variables
%---------------------------------------------------------------------%
% define the benchmark parameters from meijaard et. al, 2007

% w:        wheelbase                          [m]
% c:        trail                              [m]
% lambda:   steer axis tilt                    [rad]
% g:        gravity                            [n/kg]
% v:        forward speed                      [m/s]
% rR:       rear wheel radius                  [m]
% mR:      rear wheel mass                    [kg]
% IRxx:     rear wheel mass moment of inertia  [kg*m^2]
% IRyy:     rear wheel mass moment of inertia  [kg*m^2]
% xB:       rear body center of mass location  [m]
% zB:       rear body center of mass location  [m]
% mB:      rear body mass                     [kg]
% IBxx:     rear body mass moment of inertia   [kg*m^2]
% IByy:     rear body mass moment of inertia   [kg*m^2]
% IBzz:     rear body mass moment of inertia   [kg*m^2]
% IBxz:     rear body mass moment of inertia   [kg*m^2]
% xH:       fork center of mass location       [m]
% zH:       fork center of mass location       [m]
% mH:      fork mass                          [kg]
% IHxx:     fork mass moment of inertia        [kg*m^2]
% IHyy:     fork mass moment of inertia        [kg*m^2]
% IHzz:     fork mass moment of inertia        [kg*m^2]
% IHxz:     fork mass moment of inertia        [kg*m^2]
% rF:       front wheel radius                 [m]
% mF:      front wheel mass                   [kg]
% IFxx:     front wheel mass moment of inertia [kg*m^2]
% IFyy:     front wheel mass moment of inertia [kg*m^2]
% Tphi:    roll torque                        [n*m]
% Tdelta:  steer torque                       [n*m]
% TthetaR: rear wheel torque                  [n*m]

constants w,c,lambda,g,v
constants rR,mR,IRxx,IRyy
constants xB,zB,mB,IBxx,IByy,IBzz,IBxz
constants xH,zH,mH,IHxx,IHyy,IHzz,IHxz
constants rF,mF,IFxx,IFyy
specified Tphi,Tdelta,TthetaR

% convert the benchmark constants to this model's constants
% rF: radius of front wheel
% rR: radius of rear wheel
% d1: the perpendicular distance from the head tube axis to the center
%     of the rear wheel
% d3: the perpendicular distance from the head tube axis to the center
%     of the front wheel (fork offset)
% l1: the distance in the d1> direction from the center of the rear
%     wheel to the frame center of mass
% l2: the distance in the d3> direction from the center of the rear
%     wheel to the frame center of mass
% l3: the distance in the f1> direction from the steer point to the
%     center of mass of the fork
% l4: the distance in the f3> direction from the steer point to the
%     center of mass of the fork

d1    =  cos(lambda)*(c+w-rR*tan(lambda))
d3    = -cos(lambda)*(c-rF*tan(lambda))
d2    = (rR+d1*sin(lambda)-rF+d3*sin(lambda))/cos(lambda)
% rear wheel inertia
id11  =  IRxx
id22  =  IRyy
id33  =  IRxx
cf    =  [cos(lambda),0,-sin(lambda);0,1,0;sin(lambda),0,cos(lambda)]
% rotate bicycle frame inertia through lambda
IB    =  [IBxx,0,IBxz;0,IByy,0;IBxz,0,IBzz]
IBrot =  cf*IB*transpose(cf)
% bicycle frame inertia
ic11  =  IBrot[1,1]
ic12  =  IBrot[1,2]
ic22  =  IBrot[2,2]
ic23  =  IBrot[2,3]
ic31  =  IBrot[3,1]
ic33  =  IBrot[3,3]
% rotate fork inertia matrix through lambda
IH    =  [IHxx,0,IHxz;0,IHyy,0;IHxz,0,IHzz]
IHrot =  cf*IH*transpose(cf)
% fork/handlebar inertia
ie11  =  IHrot[1,1]
ie12  =  IHrot[1,2]
ie22  =  IHrot[2,2]
ie23  =  IHrot[2,3]
ie31  =  IHrot[3,1]
ie33  =  IHrot[3,3]
% front wheel inertia
if11  =  IFxx
if22  =  IFyy
if33  =  IFxx
% mass center locations
l1    =  xB*cos(lambda)-zB*sin(lambda)-rR*sin(lambda)
l2    =  xB*sin(lambda)+zB*cos(lambda)+rR*cos(lambda)
l3    =  cos(lambda)*xH-sin(lambda)*zH-c*cos(lambda)-w*cos(lambda)
l4    =  rR*cos(lambda)+xH*sin(lambda)+zH*cos(lambda)
% masses
massc    =  mB
massd    =  mR
masse    =  mH
massf    =  mF
% input torques
T4    =  Tphi
T6    =  TthetaR
T7    =  Tdelta

% add stuff for the lateral roll disturbance force
% pf: point at which the lateral force is applied
points pf
% location of the point
% xpf: the distance from the rear wheel contact point to the pull force point
% zpf: the distance from the rear wheel contact point to the pull force point
constants xpf, zpf
l5    =  xpf*cos(lambda)-zpf*sin(lambda)-rR*sin(lambda)
l6    =  xpf*sin(lambda)+zpf*cos(lambda)+rR*cos(lambda)
% Fphi:    lateral roll disturbance force     [n*m]
specified Fphi

% declare the generalized coordinates
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

mass c=massc,d=massd,e=masse,f=massf
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
simprot(b,c,2,q5+lambda)

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
p_do_eo>=d1*c1>+l3*e1>+l4*e3> % rear wheel center to fork/handlebar center

% rear wheel center to the front wheel center
p_do_fo>=d1*c1>+d2*e3>+d3*e1>

% locate the ground contact points
p_do_dn>=rr*b3>
p_dn_nd>=0>
p_fo_fn>=rF*unitvec(n3>-dot(e2>,n3>)*e2>)
p_fn_nf>=0>

% locate the pull force point
p_do_pf>=l5*c1>+l6*c3> % rear wheel center to pull force point

%---------------------------------------------------------------------%
%         define the generalized speeds
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

v_co_n>=dt(p_no_co>,n)
v_do_n>=dt(p_no_do>,n)
v_eo_n>=dt(p_no_eo>,n)
v_fo_n>=dt(p_no_fo>,n)

v2pts(n,d,do,dn)
v2pts(n,f,fo,fn)

v_pf_n>=dt(p_no_pf>,n)

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

% the rear wheel angular speed, u6, the roll rate, u4,the
% steering rate, u7, and the rider lean rate, u9 are taken to be the independent
% generalized speeds

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
force_pf>+=Fphi*(0*n1>+1*n2>)

%---------------------------------------------------------------------%
%         equations of motion
%---------------------------------------------------------------------%

zero=fr()+frstar()
solve(zero,u4',u6',u7')

%---------------------------------------------------------------------%
%       some extra outputs
%---------------------------------------------------------------------%

% front wheel contact location
q9 = dot(p_no_nf>, n1>)
q10 = dot(p_no_nf>, n2>)

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
bMat[1,3]=d(q1',Fphi)

bMat[2,1]=d(q2',T4)
bMat[2,2]=d(q2',T7)
bMat[2,3]=d(q2',Fphi)

bMat[3,1]=d(q3',T4)
bMat[3,2]=d(q3',T7)
bMat[3,3]=d(q3',Fphi)

bMat[4,1]=d(q4',T4)
bMat[4,2]=d(q4',T7)
bMat[4,3]=d(q4',Fphi)

bMat[5,1]=d(q5',T4)
bMat[5,2]=d(q5',T7)
bMat[5,3]=d(q5',Fphi)

bMat[6,1]=d(q6',T4)
bMat[6,2]=d(q6',T7)
bMat[6,3]=d(q6',Fphi)

bMat[7,1]=d(q7',T4)
bMat[7,2]=d(q7',T7)
bMat[7,3]=d(q7',Fphi)

bMat[8,1]=d(q8',T4)
bMat[8,2]=d(q8',T7)
bMat[8,3]=d(q8',Fphi)

bMat[9,1]=d(u4',T4)
bMat[9,2]=d(u4',T7)
bMat[9,3]=d(u4',Fphi)

bMat[10,1]=d(u6',T4)
bMat[10,2]=d(u6',T7)
bMat[10,3]=d(u6',Fphi)

bMat[11,1]=d(u7',T4)
bMat[11,2]=d(u7',T7)
bMat[11,3]=d(u7',Fphi)

cMat[1,1]=d(q9,q1)
cMat[1,2]=d(q9,q2)
cMat[1,3]=d(q9,q3)
cMat[1,4]=d(q9,q4)
cMat[1,5]=d(q9,q5)
cMat[1,6]=d(q9,q6)
cMat[1,7]=d(q9,q7)
cMat[1,8]=d(q9,q8)
cMat[1,9]=d(q9,u4)
cMat[1,10]=d(q9,u6)
cMat[1,11]=d(q9,u7)

cMat[2,1]=d(q10,q1)
cMat[2,2]=d(q10,q2)
cMat[2,3]=d(q10,q3)
cMat[2,4]=d(q10,q4)
cMat[2,5]=d(q10,q5)
cMat[2,6]=d(q10,q6)
cMat[2,7]=d(q10,q7)
cMat[2,8]=d(q10,q8)
cMat[2,9]=d(q10,u4)
cMat[2,10]=d(q10,u6)
cMat[2,11]=d(q10,u7)

dMat[1,1]=d(q9,T4)
dMat[1,2]=d(q9,T7)
dMat[1,3]=d(q9,Fphi)

dMat[2,1]=d(q10,T4)
dMat[2,2]=d(q10,T7)
dMat[2,3]=d(q10,Fphi)

encode aMat,bMat,cMat,dMat

code dynamics() WhipplePullForce.c
code algebraic() WhipplePullForce.m

%---------------------------------------------------------------------%
%         save output
%---------------------------------------------------------------------%

save WhipplePullForce.all

%---------------------------------------------------------------------%
