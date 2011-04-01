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
% rr:       rear wheel radius                  [m]
% m_r:      rear wheel mass                    [kg]
% irxx:     rear wheel mass moment of inertia  [kg*m^2]
% iryy:     rear wheel mass moment of inertia  [kg*m^2]
% xb:       rear body center of mass location  [m]
% zb:       rear body center of mass location  [m]
% m_b:      rear body mass                     [kg]
% ibxx:     rear body mass moment of inertia   [kg*m^2]
% ibyy:     rear body mass moment of inertia   [kg*m^2]
% ibzz:     rear body mass moment of inertia   [kg*m^2]
% ibxz:     rear body mass moment of inertia   [kg*m^2]
% xh:       fork center of mass location       [m]
% zh:       fork center of mass location       [m]
% m_h:      fork mass                          [kg]
% ihxx:     fork mass moment of inertia        [kg*m^2]
% ihyy:     fork mass moment of inertia        [kg*m^2]
% ihzz:     fork mass moment of inertia        [kg*m^2]
% ihxz:     fork mass moment of inertia        [kg*m^2]
% rf:       front wheel radius                 [m]
% m_f:      front wheel mass                   [kg]
% ifxx:     front wheel mass moment of inertia [kg*m^2]
% ifyy:     front wheel mass moment of inertia [kg*m^2]
% t_phi:    roll torque                        [n*m]
% t_delta:  steer torque                       [n*m]
% t_thetar: rear wheel torque                  [n*m]

constants w,c,lambda,g,v
constants rr,m_r,irxx,iryy
constants xb,zb,m_b,ibxx,ibyy,ibzz,ibxz
constants xh,zh,m_h,ihxx,ihyy,ihzz,ihxz
constants rf,m_f,ifxx,ifyy
specified t_phi,t_delta,t_thetar

% convert the benchmark constants to this model's constants
% rf: radius of front wheel
% rr: radius of rear wheel
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

d1    =  cos(lambda)*(c+w-rr*tan(lambda))
d3    = -cos(lambda)*(c-rf*tan(lambda))
d2    = (rr+d1*sin(lambda)-rf+d3*sin(lambda))/cos(lambda)
% rear wheel inertia
id11  =  irxx
id22  =  iryy
id33  =  irxx
cf    =  [cos(lambda),0,-sin(lambda);0,1,0;sin(lambda),0,cos(lambda)]
% rotate bicycle frame inertia through lambda
ib    =  [ibxx,0,ibxz;0,ibyy,0;ibxz,0,ibzz]
ibrot =  cf*ib*transpose(cf)
% bicycle frame inertia
ic11  =  ibrot[1,1]
ic12  =  ibrot[1,2]
ic22  =  ibrot[2,2]
ic23  =  ibrot[2,3]
ic31  =  ibrot[3,1]
ic33  =  ibrot[3,3]
% rotate fork inertia matrix through lambda
ih    =  [ihxx,0,ihxz;0,ihyy,0;ihxz,0,ihzz]
ihrot =  cf*ih*transpose(cf)
% fork/handlebar inertia
ie11  =  ihrot[1,1]
ie12  =  ihrot[1,2]
ie22  =  ihrot[2,2]
ie23  =  ihrot[2,3]
ie31  =  ihrot[3,1]
ie33  =  ihrot[3,3]
% front wheel inertia
if11  =  ifxx
if22  =  ifyy
if33  =  ifxx
% mass center locations
l1    =  xb*cos(lambda)-zb*sin(lambda)-rr*sin(lambda)
l2    =  xb*sin(lambda)+zb*cos(lambda)+rr*cos(lambda)
l3    =  cos(lambda)*xh-sin(lambda)*zh-c*cos(lambda)-w*cos(lambda)
l4    =  rr*cos(lambda)+xh*sin(lambda)+zh*cos(lambda)
% masses
mc    =  m_b
md    =  m_r
me    =  m_h
mf    =  m_f
% input torques
t4    =  t_phi
t6    =  t_thetar
t7    =  t_delta

% add stuff for the lateral roll disturbance force
% pf: point at which the lateral force is applied
points pf
% location of the point
% xpf: the distance from the rear wheel contact point to the pull force point
% zpf: the distance from the rear wheel contact point to the pull force point
constants xpf, zpf
l5    =  xpf*cos(lambda)-zpf*sin(lambda)-rr*sin(lambda)
l6    =  xpf*sin(lambda)+zpf*cos(lambda)+rr*cos(lambda)
% f_phi:    lateral roll disturbance force     [n*m]
specified f_phi

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
p_no_do>=q1*n1>+q2*n2>-rr*b3> % newtonian origin to rear wheel center
p_do_co>=l1*c1>+l2*c3> % rear wheel center to bicycle frame center
p_do_eo>=d1*c1>+l3*e1>+l4*e3> % rear wheel center to fork/handlebar center

% rear wheel center to the front wheel center
p_do_fo>=d1*c1>+d2*e3>+d3*e1>

% locate the ground contact points
p_do_dn>=rr*b3>
p_dn_nd>=0>
p_fo_fn>=rf*unitvec(n3>-dot(e2>,n3>)*e2>)
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
torque(a/b,t4*a1>) % roll torque
torque(c/d,t6*c2>) % rear wheel torque
torque(c/e,t7*c3>) % steer torque
force_pf>+=f_phi*(0*n1>+1*n2>)

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
% linearizes the equations of motion about the upright configuration
% and constant forward speed

a[1,1]=d(q1',q1)
a[1,2]=d(q1',q2)
a[1,3]=d(q1',q3)
a[1,4]=d(q1',q4)
a[1,5]=d(q1',q5)
a[1,6]=d(q1',q6)
a[1,7]=d(q1',q7)
a[1,8]=d(q1',q8)
a[1,9]=d(q1',u4)
a[1,10]=d(q1',u6)
a[1,11]=d(q1',u7)

a[1,1]=d(q4',q4)
a[1,2]=d(q4',q7)
a[1,3]=d(q4',u4)
a[1,4]=d(q4',u7)

a[2,1]=d(q7',q4)
a[2,2]=d(q7',q7)
a[2,3]=d(q7',u4)
a[2,4]=d(q7',u7)

a[3,1]=d(u4',q4)
a[3,2]=d(u4',q7)
a[3,3]=d(u4',u4)
a[3,4]=d(u4',u7)

a[4,1]=d(u7',q4)
a[4,2]=d(u7',q7)
a[4,3]=d(u7',u4)
a[4,4]=d(u7',u7)

b[1,1]=d(q4',t4)
b[1,2]=d(q4',t7)
b[1,3]=d(q4',f_phi)

b[2,1]=d(q7',t4)
b[2,2]=d(q7',t7)
b[2,3]=d(q7',f_phi)

b[3,1]=d(u4',t4)
b[3,2]=d(u4',t7)
b[3,3]=d(u4',f_phi)

b[4,1]=d(u7',t4)
b[4,2]=d(u7',t7)
b[4,3]=d(u7',f_phi)

cmat[1,1]=d(q4',q1)
cmat[1,2]=d(q4',q2)
cmat[1,3]=d(q4',q3)
cmat[1,4]=d(q4',q4)
cmat[1,5]=d(q4',q5)
cmat[1,6]=d(q4',q6)
cmat[1,7]=d(q4',q7)
cmat[1,8]=d(q4',q8)
cmat[1,9]=d(q4',q9)
cmat[1,10]=d(q4',q10)

cmat[2,1]=d(q7',q1)
cmat[2,2]=d(q7',q2)
cmat[2,3]=d(q7',q3)
cmat[2,4]=d(q7',q4)
cmat[2,5]=d(q7',q5)
cmat[2,6]=d(q7',q6)
cmat[2,7]=d(q7',q7)
cmat[2,8]=d(q7',q8)
cmat[2,9]=d(q7',q9)
cmat[2,10]=d(q7',q10)

cmat[3,1]=d(u4',q1)
cmat[3,2]=d(u4',q2)
cmat[3,3]=d(u4',q3)
cmat[3,4]=d(u4',q4)
cmat[3,5]=d(u4',q5)
cmat[3,6]=d(u4',q6)
cmat[3,7]=d(u4',q7)
cmat[3,8]=d(u4',q8)
cmat[3,9]=d(u4',q9)
cmat[3,10]=d(u4',q10)

cmat[4,1]=d(u7',q1)
cmat[4,2]=d(u7',q2)
cmat[4,3]=d(u7',q3)
cmat[4,4]=d(u7',q4)
cmat[4,5]=d(u7',q5)
cmat[4,6]=d(u7',q6)
cmat[4,7]=d(u7',q7)
cmat[4,8]=d(u7',q8)
cmat[4,9]=d(u7',q9)
cmat[4,10]=d(u7',q10)

encode a,b,cmat

code dynamics() WhipplePullForce.c
code algebraic() WhipplePullForce.m

%---------------------------------------------------------------------%
%         save output
%---------------------------------------------------------------------%

save WhipplePullForce.all

%---------------------------------------------------------------------%
