%---------------------------------------------------------------------%
% File: WhippleArmsLateralForce.al
% Creation Date: September 15, 2011
% Author: Jason K. Moore
% Description: Generates the nonlinear and linear equations of motion
% for the Whipple bicycle model with four inputs: roll torque, rear
% wheel torque, steer torque and a lateral force at a point on the
% bicycle frame. This model also includes four rigid bodies that models
% the upper and lower arms of the rider and their motion with respect
% to the steer angle.

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

frames a, b

% declare the four bicycle bodies
% c: bicycle frame
% d: rear wheel
% e: fork/handlebar
% f: front wheel

bodies c, d, e, f

% declare the four rider bodies
% g: right upper arm
% h: right lower arm
% i: left upper arm
% j: left lower arm

bodies g, h, i, j

% declare four points to locate the wheel contact points
% nd: rear contact point on ground
% dn: rear contact point on wheel
% nf: front contact point on ground
% fn: front contact point on wheel

points nd, dn, nf, fn

% pf: point at which the lateral force is applied

points pf

% declare points to locate the rider's arms
% sr: right shoulder
% sl: left shoulder
% er : right elbow
% el : left elbow
% hr: right hand
% hl: left hand
% gr: right handlebar grip
% gl: left handlebar grip

points sr, sl, er, el, hr, hl, gr, gl

%---------------------------------------------------------------------%
%         constants and variables
%---------------------------------------------------------------------%

% gravity

constants g

% input torques and forces
% T4: roll torque [n*m]
% T6: rear wheel torque [n*m]
% T7: steer torque [n*m]
% F: lateral roll disturbance force [n]

specified T4, T6, T7, F

% rF: radius of front wheel
% rR: radius of rear wheel
% d1: the perpendicular distance from the steer axis to the center
%     of the rear wheel
% d2: the distance along the steer axis
% d3: the perpendicular distance from the steer axis to the center
%     of the front wheel (fork offset)
% d4: the distance from the rear wheel center to the lateral force point
% d5: the distance from the rear wheel center to the lateral force point
% d6: locates the shoulders relative to the rear wheel center
% d7: locates the shoulders relative to the rear wheel center
% d8: locates the shoulders relative to the rear wheel center
% d9: locates the grips relative to the rear wheel center
% d10: locates the grips relative to the rear wheel center
% d11: locates the grips relative to the rear wheel center
% d12: length of right upper arm
% d13: lenth of right lower arm
% d14: length of left upper arm
% d15: length of left lower arm

constants rF, rR, d{1:15}

% l1: the distance in the c1> direction from the center of the rear
%     wheel to the frame center of mass
% l2: the distance in the c3> direction from the center of the rear
%     wheel to the frame center of mass
% l3: the distance in the e1> direction from the steer point to the
%     center of mass of the fork
% l4: the distance in the e3> direction from the steer point to the
%     center of mass of the fork
% l5: right upper arm com distance
% l6: right lower arm com distance
% l7: left upper arm com distance
% l8: left lower arm com distance

constants l{1:8}

% bicycle masses
constants mc, md, me, mf

% arm masses
constants mg, mh, mi, mj

% bicycle inertia
constants ic11, ic22, ic33, ic12, ic23, ic31
constants id11, id22, id33
constants ie11, ie22, ie33, ie12, ie23, ie31
constants if11, if22, if33

% arm inertia
constants ig11, ig22, ig33, ig12, ig23, ig31
constants ih11, ih22, ih33, ih12, ih23, ih31
constants ii11, ii22, ii33, ii12, ii23, ii31
constants ij11, ij22, ij33, ij12, ij23, ij31

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
% q9:  right shoulder abduction
% q10: right shoulder elevation
% q11: right shoulder rotation
% q12: right elbow elevation
% q13: left shoulder abduction
% q14: left shoulder elevation
% q15: left shoulder rotation
% q16: left elbow elevation

variables q{16}'

%---------------------------------------------------------------------%
%         generalized speeds
%---------------------------------------------------------------------%

motionvariables' u{16}'

%---------------------------------------------------------------------%
%         mass and inertia properties
%---------------------------------------------------------------------%

mass c=mc, d=md, e=me, f=mf, g=mg, h=mh, i=mi, j=mj
inertia c, ic11, ic22, ic33, ic12, ic23, ic31
inertia d, id11, id22, id33
inertia e, ie11, ie22, ie33, ie12, ie23, ie31
inertia f, if11, if22, if33
inertia g, ig11, ig22, ig33, ig12, ig23, ig31
inertia h, ih11, ih22, ih33, ih12, ih23, ih31
inertia i, ii11, ii22, ii33, ii12, ii23, ii31
inertia j, ij11, ij22, ij33, ij12, ij23, ij31

%---------------------------------------------------------------------%
%         angular relationships                                       %
%---------------------------------------------------------------------%

% yaw frame
simprot(n, a, 3, q3)

% roll frame
simprot(a, b, 1, q4)

% pitch frame
simprot(b, c, 2, q5)

% rear wheel rotation
simprot(c, d, 2, q6)

% steering angle
simprot(c, e, 3, q7)

% front wheel rotation
simprot(e, f, 2, q8)

% right shoulder
dircos(c, g, body123, q9, q10, q11)

% right elbow
simprot(g, h, 2, q12)

% left shoulder
dircos(c, i, body123, q13, q14, q15)

% left elbow
simprot(i, j, 2, q16)

%---------------------------------------------------------------------%
%         position vectors
%---------------------------------------------------------------------%

% locate the center of mass for each bicycle body
p_no_do> = q1*n1> + q2*n2> - rR*b3> % newtonian origin to rear wheel center
p_do_co> = l1*c1> + l2*c3> % rear wheel center to bicycle frame center
p_do_eo> = d1*c1> + l3*e1> + l4*e3> % rear wheel center to fork/handlebar center
p_do_fo> = d1*c1> + d2*e3> + d3*e1> % rear wheel center to the front wheel center

% locate the ground contact points
p_do_dn> = rr * b3>
p_dn_nd> = 0>
p_fo_fn> = rF * unitvec(n3> - dot(e2>, n3>) * e2>)
p_fn_nf> = 0>

% locate the lateral force point
p_do_pf> = d4 * c1> + d5 * c3> % rear wheel center to lateral force point

% locate the shoulders
p_do_sr> = d6 * c1> + d7 * c2> + d8 * c3>
p_do_sl> = d6 * c1> - d7 * c2> + d8 * c3>

% locate the grips
p_fo_gr> = d9 * c1> + d10 * c2> + d11 * c3>
p_fo_gl> = d9 * c1> - d10 * c2> + d11 * c3>

% locate the elbows
p_sr_er> = d12 * g3>
p_sl_el> = d14 * i3>

% locate the hands
p_er_hr> = d13 * h3>
p_el_hl> = d15 * j3>

% locate the centers of mass of the arms
p_sr_go> = l5 * g3>
p_er_ho> = l6 * h3>
p_sl_io> = l7 * i3>
p_el_jo> = l8 * j3>

%---------------------------------------------------------------------%
%         define the kinematical differential equations
%---------------------------------------------------------------------%

q1' = u1
q2' = u2
q3' = u3
q4' = u4
q5' = u5
q6' = u6
q7' = u7
q8' = u8
q9' = u9
q10' = u10
q11' = u11
q12' = u12
q13' = u13
q14' = u14
q15' = u15
q16' = u16

%---------------------------------------------------------------------%
%         angular velocities
%---------------------------------------------------------------------%

angvel(n, a)
angvel(n, b)
angvel(n, c)
angvel(n, d)
angvel(n, e)
angvel(n, f)
angvel(n, g)
angvel(n, h)
angvel(n, i)
angvel(n, j)

%---------------------------------------------------------------------%
%         velocities
%---------------------------------------------------------------------%

% center of mass velocities
v_co_n> = dt(p_no_co>, n)
v_do_n> = dt(p_no_do>, n)
v_eo_n> = dt(p_no_eo>, n)
v_fo_n> = dt(p_no_fo>, n)
v_go_n> = dt(p_no_go>, n)
v_ho_n> = dt(p_no_ho>, n)
v_io_n> = dt(p_no_io>, n)
v_jo_n> = dt(p_no_jo>, n)

% velocity of the wheel contacts
v2pts(n, d, do, dn)
v2pts(n, f, fo, fn)

% velocity of the lateral force point
v_pf_n> = dt(p_no_pf>, n)

%---------------------------------------------------------------------%
%         define the configuration constraints
%---------------------------------------------------------------------%

% bicycle frame pitch, both wheels must touch the ground
% set the a3> component of p_nd_nf> equal to zero
pzero = dot(p_nd_nf>, a3>)

% hands must touch the grips
hrzero> = p_sr_hr> - p_sr_gr>
hlzero> = p_sl_hl> - p_sl_gl>

% arms must always hang down
rzero = dot(g2>, b3>)
lzero = dot(i2>, b3>)

%---------------------------------------------------------------------%
%         motion constraints
%---------------------------------------------------------------------%

% non-holonomic
% due to the assumptions of no side slip and no slip rolling the
% velocities of the front and rear wheel contact points, cn and gn,
% cannot have components of velocity in the ground plane
dependent[1] = dot(v_dn_n>, a1>)
dependent[2] = dot(v_dn_n>, a2>)
dependent[3] = dot(v_fn_n>, a1>)
dependent[4] = dot(v_fn_n>, a2>)

% declare the derivative of each configuration constraint (holonomic)
dependent[5] = dt(pzero)
dependent[6] = dt(dot(hrzero>, c1>))
dependent[7] = dt(dot(hrzero>, c2>))
dependent[8] = dt(dot(hrzero>, c3>))
dependent[9] = dt(dot(hlzero>, c1>))
dependent[10] = dt(dot(hlzero>, c2>))
dependent[11] = dt(dot(hlzero>, c3>))
dependent[12] = dt(rzero)
dependent[13] = dt(lzero)

% the rear wheel angular speed, u6, the roll rate, u4, and the steering rate,
% u7, are taken to be the independent generalized speeds

constrain(dependent[u1, u2, u3, u5, u8, u9, u10, u11, u12, u13, u14, u15, u16])

%---------------------------------------------------------------------%
%         angular accelerations
%---------------------------------------------------------------------%

alf_c_n> = dt(w_c_n>, n)
alf_d_n> = dt(w_d_n>, n)
alf_e_n> = dt(w_e_n>, n)
alf_f_n> = dt(w_f_n>, n)
alf_g_n> = dt(w_g_n>, n)
alf_h_n> = dt(w_h_n>, n)
alf_i_n> = dt(w_i_n>, n)
alf_j_n> = dt(w_j_n>, n)

%---------------------------------------------------------------------%
%         accelerations
%---------------------------------------------------------------------%

a_co_n> = dt(v_co_n>, n)
a_do_n> = dt(v_do_n>, n)
a_eo_n> = dt(v_eo_n>, n)
a_fo_n> = dt(v_fo_n>, n)
a_go_n> = dt(v_go_n>, n)
a_ho_n> = dt(v_ho_n>, n)
a_io_n> = dt(v_io_n>, n)
a_jo_n> = dt(v_jo_n>, n)

%---------------------------------------------------------------------%
%         forces and torques
%---------------------------------------------------------------------%

gravity(g * a3>, c, d, e, f, g, h, i, j)
torque(a/b, T4 * a1>) % roll torque
torque(c/d, T6 * c2>) % rear wheel torque
torque(c/e, T7 * c3>) % steer torque
force_pf> += F * n2> % lateral force

%---------------------------------------------------------------------%
%         equations of motion
%---------------------------------------------------------------------%

zero = fr() + frstar()
solve(zero, u4', u6', u7')

%---------------------------------------------------------------------%
%       some extra outputs
%---------------------------------------------------------------------%

% front wheel contact location
q17 = dot(p_no_nf>, n1>)
q18 = dot(p_no_nf>, n2>)

%---------------------------------------------------------------------%
%         linearization
%---------------------------------------------------------------------%
% linearizes the equations of motion

%aMat[1, 1] = d(q1', q1)
%aMat[1, 2] = d(q1', q2)
%aMat[1, 3] = d(q1', q3)
%aMat[1, 4] = d(q1', q4)
%aMat[1, 5] = d(q1', q5)
%aMat[1, 6] = d(q1', q6)
%aMat[1, 7] = d(q1', q7)
%aMat[1, 8] = d(q1', q8)
%aMat[1, 9] = d(q1', q9)
%aMat[1, 10] = d(q1', q10)
%aMat[1, 11] = d(q1', q11)
%aMat[1, 12] = d(q1', q12)
%aMat[1, 13] = d(q1', q13)
%aMat[1, 14] = d(q1', q14)
%aMat[1, 15] = d(q1', q15)
%aMat[1, 16] = d(q1', q16)
%aMat[1, 17] = d(q1', u4)
%aMat[1, 18] = d(q1', u6)
%aMat[1, 19] = d(q1', u7)
%
%aMat[2, 1] = d(q2', q1)
%aMat[2, 2] = d(q2', q2)
%aMat[2, 3] = d(q2', q3)
%aMat[2, 4] = d(q2', q4)
%aMat[2, 5] = d(q2', q5)
%aMat[2, 6] = d(q2', q6)
%aMat[2, 7] = d(q2', q7)
%aMat[2, 8] = d(q2', q8)
%aMat[2, 9] = d(q2', q9)
%aMat[2, 10] = d(q2', q10)
%aMat[2, 11] = d(q2', q11)
%aMat[2, 12] = d(q2', q12)
%aMat[2, 13] = d(q2', q13)
%aMat[2, 14] = d(q2', q14)
%aMat[2, 15] = d(q2', q15)
%aMat[2, 16] = d(q2', q16)
%aMat[2, 17] = d(q2', u4)
%aMat[2, 18] = d(q2', u6)
%aMat[2, 19] = d(q2', u7)
%
%aMat[3, 1] = d(q3', q1)
%aMat[3, 2] = d(q3', q2)
%aMat[3, 3] = d(q3', q3)
%aMat[3, 4] = d(q3', q4)
%aMat[3, 5] = d(q3', q5)
%aMat[3, 6] = d(q3', q6)
%aMat[3, 7] = d(q3', q7)
%aMat[3, 8] = d(q3', q8)
%aMat[3, 9] = d(q3', q9)
%aMat[3, 10] = d(q3', q10)
%aMat[3, 11] = d(q3', q11)
%aMat[3, 12] = d(q3', q12)
%aMat[3, 13] = d(q3', q13)
%aMat[3, 14] = d(q3', q14)
%aMat[3, 15] = d(q3', q15)
%aMat[3, 16] = d(q3', q16)
%aMat[3, 17] = d(q3', u4)
%aMat[3, 18] = d(q3', u6)
%aMat[3, 19] = d(q3', u7)
%
%aMat[4, 1] = d(q4', q1)
%aMat[4, 2] = d(q4', q2)
%aMat[4, 3] = d(q4', q3)
%aMat[4, 4] = d(q4', q4)
%aMat[4, 5] = d(q4', q5)
%aMat[4, 6] = d(q4', q6)
%aMat[4, 7] = d(q4', q7)
%aMat[4, 8] = d(q4', q8)
%aMat[4, 9] = d(q4', q9)
%aMat[4, 10] = d(q4', q10)
%aMat[4, 11] = d(q4', q11)
%aMat[4, 12] = d(q4', q12)
%aMat[4, 13] = d(q4', q13)
%aMat[4, 14] = d(q4', q14)
%aMat[4, 15] = d(q4', q15)
%aMat[4, 16] = d(q4', q16)
%aMat[4, 17] = d(q4', u4)
%aMat[4, 18] = d(q4', u6)
%aMat[4, 19] = d(q4', u7)
%
%aMat[5, 1] = d(q5', q1)
%aMat[5, 2] = d(q5', q2)
%aMat[5, 3] = d(q5', q3)
%aMat[5, 4] = d(q5', q4)
%aMat[5, 5] = d(q5', q5)
%aMat[5, 6] = d(q5', q6)
%aMat[5, 7] = d(q5', q7)
%aMat[5, 8] = d(q5', q8)
%aMat[5, 9] = d(q5', q9)
%aMat[5, 10] = d(q5', q10)
%aMat[5, 11] = d(q5', q11)
%aMat[5, 12] = d(q5', q12)
%aMat[5, 13] = d(q5', q13)
%aMat[5, 14] = d(q5', q14)
%aMat[5, 15] = d(q5', q15)
%aMat[5, 16] = d(q5', q16)
%aMat[5, 17] = d(q5', u4)
%aMat[5, 18] = d(q5', u6)
%aMat[5, 19] = d(q5', u7)
%
%aMat[6, 1] = d(q6', q1)
%aMat[6, 2] = d(q6', q2)
%aMat[6, 3] = d(q6', q3)
%aMat[6, 4] = d(q6', q4)
%aMat[6, 5] = d(q6', q5)
%aMat[6, 6] = d(q6', q6)
%aMat[6, 7] = d(q6', q7)
%aMat[6, 8] = d(q6', q8)
%aMat[6, 9] = d(q6', q9)
%aMat[6, 10] = d(q6', q10)
%aMat[6, 11] = d(q6', q11)
%aMat[6, 12] = d(q6', q12)
%aMat[6, 13] = d(q6', q13)
%aMat[6, 14] = d(q6', q14)
%aMat[6, 15] = d(q6', q15)
%aMat[6, 16] = d(q6', q16)
%aMat[6, 17] = d(q6', u4)
%aMat[6, 18] = d(q6', u6)
%aMat[6, 19] = d(q6', u7)
%
%aMat[7, 1] = d(q7', q1)
%aMat[7, 2] = d(q7', q2)
%aMat[7, 3] = d(q7', q3)
%aMat[7, 4] = d(q7', q4)
%aMat[7, 5] = d(q7', q5)
%aMat[7, 6] = d(q7', q6)
%aMat[7, 7] = d(q7', q7)
%aMat[7, 8] = d(q7', q8)
%aMat[7, 9] = d(q7', q9)
%aMat[7, 10] = d(q7', q10)
%aMat[7, 11] = d(q7', q11)
%aMat[7, 12] = d(q7', q12)
%aMat[7, 13] = d(q7', q13)
%aMat[7, 14] = d(q7', q14)
%aMat[7, 15] = d(q7', q15)
%aMat[7, 16] = d(q7', q16)
%aMat[7, 17] = d(q7', u4)
%aMat[7, 18] = d(q7', u6)
%aMat[7, 19] = d(q7', u7)
%
%aMat[8, 1] = d(q8', q1)
%aMat[8, 2] = d(q8', q2)
%aMat[8, 3] = d(q8', q3)
%aMat[8, 4] = d(q8', q4)
%aMat[8, 5] = d(q8', q5)
%aMat[8, 6] = d(q8', q6)
%aMat[8, 7] = d(q8', q7)
%aMat[8, 8] = d(q8', q8)
%aMat[8, 9] = d(q8', q9)
%aMat[8, 10] = d(q8', q10)
%aMat[8, 11] = d(q8', q11)
%aMat[8, 12] = d(q8', q12)
%aMat[8, 13] = d(q8', q13)
%aMat[8, 14] = d(q8', q14)
%aMat[8, 15] = d(q8', q15)
%aMat[8, 16] = d(q8', q16)
%aMat[8, 17] = d(q8', u4)
%aMat[8, 18] = d(q8', u6)
%aMat[8, 19] = d(q8', u7)
%
%aMat[9, 1] = d(q9', q1)
%aMat[9, 2] = d(q9', q2)
%aMat[9, 3] = d(q9', q3)
%aMat[9, 4] = d(q9', q4)
%aMat[9, 5] = d(q9', q5)
%aMat[9, 6] = d(q9', q6)
%aMat[9, 7] = d(q9', q7)
%aMat[9, 8] = d(q9', q8)
%aMat[9, 9] = d(q9', q9)
%aMat[9, 10] = d(q9', q10)
%aMat[9, 11] = d(q9', q11)
%aMat[9, 12] = d(q9', q12)
%aMat[9, 13] = d(q9', q13)
%aMat[9, 14] = d(q9', q14)
%aMat[9, 15] = d(q9', q15)
%aMat[9, 16] = d(q9', q16)
%aMat[9, 17] = d(q9', u4)
%aMat[9, 18] = d(q9', u6)
%aMat[9, 19] = d(q9', u7)
%
%aMat[10, 1] = d(q10', q1)
%aMat[10, 2] = d(q10', q2)
%aMat[10, 3] = d(q10', q3)
%aMat[10, 4] = d(q10', q4)
%aMat[10, 5] = d(q10', q5)
%aMat[10, 6] = d(q10', q6)
%aMat[10, 7] = d(q10', q7)
%aMat[10, 8] = d(q10', q8)
%aMat[10, 9] = d(q10', q9)
%aMat[10, 10] = d(q10', q10)
%aMat[10, 11] = d(q10', q11)
%aMat[10, 12] = d(q10', q12)
%aMat[10, 13] = d(q10', q13)
%aMat[10, 14] = d(q10', q14)
%aMat[10, 15] = d(q10', q15)
%aMat[10, 16] = d(q10', q16)
%aMat[10, 17] = d(q10', u4)
%aMat[10, 18] = d(q10', u6)
%aMat[10, 19] = d(q10', u7)
%
%aMat[11, 1] = d(q11', q1)
%aMat[11, 2] = d(q11', q2)
%aMat[11, 3] = d(q11', q3)
%aMat[11, 4] = d(q11', q4)
%aMat[11, 5] = d(q11', q5)
%aMat[11, 6] = d(q11', q6)
%aMat[11, 7] = d(q11', q7)
%aMat[11, 8] = d(q11', q8)
%aMat[11, 9] = d(q11', q9)
%aMat[11, 10] = d(q11', q10)
%aMat[11, 11] = d(q11', q11)
%aMat[11, 12] = d(q11', q12)
%aMat[11, 13] = d(q11', q13)
%aMat[11, 14] = d(q11', q14)
%aMat[11, 15] = d(q11', q15)
%aMat[11, 16] = d(q11', q16)
%aMat[11, 17] = d(q11', u4)
%aMat[11, 18] = d(q11', u6)
%aMat[11, 19] = d(q11', u7)
%
%aMat[12, 1] = d(q12', q1)
%aMat[12, 2] = d(q12', q2)
%aMat[12, 3] = d(q12', q3)
%aMat[12, 4] = d(q12', q4)
%aMat[12, 5] = d(q12', q5)
%aMat[12, 6] = d(q12', q6)
%aMat[12, 7] = d(q12', q7)
%aMat[12, 8] = d(q12', q8)
%aMat[12, 9] = d(q12', q9)
%aMat[12, 10] = d(q12', q10)
%aMat[12, 11] = d(q12', q11)
%aMat[12, 12] = d(q12', q12)
%aMat[12, 13] = d(q12', q13)
%aMat[12, 14] = d(q12', q14)
%aMat[12, 15] = d(q12', q15)
%aMat[12, 16] = d(q12', q16)
%aMat[12, 17] = d(q12', u4)
%aMat[12, 18] = d(q12', u6)
%aMat[12, 19] = d(q12', u7)
%
%aMat[13, 1] = d(q13', q1)
%aMat[13, 2] = d(q13', q2)
%aMat[13, 3] = d(q13', q3)
%aMat[13, 4] = d(q13', q4)
%aMat[13, 5] = d(q13', q5)
%aMat[13, 6] = d(q13', q6)
%aMat[13, 7] = d(q13', q7)
%aMat[13, 8] = d(q13', q8)
%aMat[13, 9] = d(q13', q9)
%aMat[13, 10] = d(q13', q10)
%aMat[13, 11] = d(q13', q11)
%aMat[13, 12] = d(q13', q12)
%aMat[13, 13] = d(q13', q13)
%aMat[13, 14] = d(q13', q14)
%aMat[13, 15] = d(q13', q15)
%aMat[13, 16] = d(q13', q16)
%aMat[13, 17] = d(q13', u4)
%aMat[13, 18] = d(q13', u6)
%aMat[13, 19] = d(q13', u7)
%
%aMat[14, 1] = d(q14', q1)
%aMat[14, 2] = d(q14', q2)
%aMat[14, 3] = d(q14', q3)
%aMat[14, 4] = d(q14', q4)
%aMat[14, 5] = d(q14', q5)
%aMat[14, 6] = d(q14', q6)
%aMat[14, 7] = d(q14', q7)
%aMat[14, 8] = d(q14', q8)
%aMat[14, 9] = d(q14', q9)
%aMat[14, 10] = d(q14', q10)
%aMat[14, 11] = d(q14', q11)
%aMat[14, 12] = d(q14', q12)
%aMat[14, 13] = d(q14', q13)
%aMat[14, 14] = d(q14', q14)
%aMat[14, 15] = d(q14', q15)
%aMat[14, 16] = d(q14', q16)
%aMat[14, 17] = d(q14', u4)
%aMat[14, 18] = d(q14', u6)
%aMat[14, 19] = d(q14', u7)
%
%aMat[15, 1] = d(q15', q1)
%aMat[15, 2] = d(q15', q2)
%aMat[15, 3] = d(q15', q3)
%aMat[15, 4] = d(q15', q4)
%aMat[15, 5] = d(q15', q5)
%aMat[15, 6] = d(q15', q6)
%aMat[15, 7] = d(q15', q7)
%aMat[15, 8] = d(q15', q8)
%aMat[15, 9] = d(q15', q9)
%aMat[15, 10] = d(q15', q10)
%aMat[15, 11] = d(q15', q11)
%aMat[15, 12] = d(q15', q12)
%aMat[15, 13] = d(q15', q13)
%aMat[15, 14] = d(q15', q14)
%aMat[15, 15] = d(q15', q15)
%aMat[15, 16] = d(q15', q16)
%aMat[15, 17] = d(q15', u4)
%aMat[15, 18] = d(q15', u6)
%aMat[15, 19] = d(q15', u7)
%
%aMat[16, 1] = d(q16', q1)
%aMat[16, 2] = d(q16', q2)
%aMat[16, 3] = d(q16', q3)
%aMat[16, 4] = d(q16', q4)
%aMat[16, 5] = d(q16', q5)
%aMat[16, 6] = d(q16', q6)
%aMat[16, 7] = d(q16', q7)
%aMat[16, 8] = d(q16', q8)
%aMat[16, 9] = d(q16', q9)
%aMat[16, 10] = d(q16', q10)
%aMat[16, 11] = d(q16', q11)
%aMat[16, 12] = d(q16', q12)
%aMat[16, 13] = d(q16', q13)
%aMat[16, 14] = d(q16', q14)
%aMat[16, 15] = d(q16', q15)
%aMat[16, 16] = d(q16', q16)
%aMat[16, 17] = d(q16', u4)
%aMat[16, 18] = d(q16', u6)
%aMat[16, 19] = d(q16', u7)
%
%aMat[17, 1] = d(u4', q1)
%aMat[17, 2] = d(u4', q2)
%aMat[17, 3] = d(u4', q3)
%aMat[17, 4] = d(u4', q4)
%aMat[17, 5] = d(u4', q5)
%aMat[17, 6] = d(u4', q6)
%aMat[17, 7] = d(u4', q7)
%aMat[17, 8] = d(u4', q8)
%aMat[17, 9] = d(u4', q9)
%aMat[17, 10] = d(u4', q10)
%aMat[17, 11] = d(u4', q11)
%aMat[17, 12] = d(u4', q12)
%aMat[17, 13] = d(u4', q13)
%aMat[17, 14] = d(u4', q14)
%aMat[17, 15] = d(u4', q15)
%aMat[17, 16] = d(u4', q16)
%aMat[17, 17] = d(u4', u4)
%aMat[17, 18] = d(u4', u6)
%aMat[17, 19] = d(u4', u7)
%
%aMat[18, 1] = d(u6', q1)
%aMat[18, 2] = d(u6', q2)
%aMat[18, 3] = d(u6', q3)
%aMat[18, 4] = d(u6', q4)
%aMat[18, 5] = d(u6', q5)
%aMat[18, 6] = d(u6', q6)
%aMat[18, 7] = d(u6', q7)
%aMat[18, 8] = d(u6', q8)
%aMat[18, 9] = d(u6', q9)
%aMat[18, 10] = d(u6', q10)
%aMat[18, 11] = d(u6', q11)
%aMat[18, 12] = d(u6', q12)
%aMat[18, 13] = d(u6', q13)
%aMat[18, 14] = d(u6', q14)
%aMat[18, 15] = d(u6', q15)
%aMat[18, 16] = d(u6', q16)
%aMat[18, 17] = d(u6', u4)
%aMat[18, 18] = d(u6', u6)
%aMat[18, 19] = d(u6', u7)
%
%aMat[19, 1] = d(u7', q1)
%aMat[19, 2] = d(u7', q2)
%aMat[19, 3] = d(u7', q3)
%aMat[19, 4] = d(u7', q4)
%aMat[19, 5] = d(u7', q5)
%aMat[19, 6] = d(u7', q6)
%aMat[19, 7] = d(u7', q7)
%aMat[19, 8] = d(u7', q8)
%aMat[19, 9] = d(u7', q9)
%aMat[19, 10] = d(u7', q10)
%aMat[19, 11] = d(u7', q11)
%aMat[19, 12] = d(u7', q12)
%aMat[19, 13] = d(u7', q13)
%aMat[19, 14] = d(u7', q14)
%aMat[19, 15] = d(u7', q15)
%aMat[19, 16] = d(u7', q16)
%aMat[19, 17] = d(u7', u4)
%aMat[19, 18] = d(u7', u6)
%aMat[19, 19] = d(u7', u7)
%
%bMat[1, 1] = d(q1', T4)
%bMat[1, 2] = d(q1', T6)
%bMat[1, 3] = d(q1', T7)
%bMat[1, 4] = d(q1', F)
%
%bMat[2, 1] = d(q2', T4)
%bMat[2, 2] = d(q2', T6)
%bMat[2, 3] = d(q2', T7)
%bMat[2, 4] = d(q2', F)
%
%bMat[3, 1] = d(q3', T4)
%bMat[3, 2] = d(q3', T6)
%bMat[3, 3] = d(q3', T7)
%bMat[3, 4] = d(q3', F)
%
%bMat[4, 1] = d(q4', T4)
%bMat[4, 2] = d(q4', T6)
%bMat[4, 3] = d(q4', T7)
%bMat[4, 4] = d(q4', F)
%
%bMat[5, 1] = d(q5', T4)
%bMat[5, 2] = d(q5', T6)
%bMat[5, 2] = d(q5', T7)
%bMat[5, 4] = d(q5', F)
%
%bMat[6, 1] = d(q6', T4)
%bMat[6, 2] = d(q6', T6)
%bMat[6, 3] = d(q6', T7)
%bMat[6, 4] = d(q6', F)
%
%bMat[7, 1] = d(q7', T4)
%bMat[7, 2] = d(q7', T6)
%bMat[7, 3] = d(q7', T7)
%bMat[7, 4] = d(q7', F)
%
%bMat[8, 1] = d(q8', T4)
%bMat[8, 2] = d(q8', T6)
%bMat[8, 3] = d(q8', T7)
%bMat[8, 4] = d(q8', F)
%
%bMat[9, 1] = d(q9', T4)
%bMat[9, 2] = d(q9', T6)
%bMat[9, 3] = d(q9', T7)
%bMat[9, 4] = d(q9', F)
%
%bMat[10, 1] = d(q10', T4)
%bMat[10, 2] = d(q10', T6)
%bMat[10, 3] = d(q10', T7)
%bMat[10, 4] = d(q10', F)
%
%bMat[11, 1] = d(q11', T4)
%bMat[11, 2] = d(q11', T6)
%bMat[11, 3] = d(q11', T7)
%bMat[11, 4] = d(q11', F)
%
%bMat[12, 1] = d(q12', T4)
%bMat[12, 2] = d(q12', T6)
%bMat[12, 3] = d(q12', T7)
%bMat[12, 4] = d(q12', F)
%
%bMat[13, 1] = d(q13', T4)
%bMat[13, 2] = d(q13', T6)
%bMat[13, 3] = d(q13', T7)
%bMat[13, 4] = d(q13', F)
%
%bMat[14, 1] = d(q14', T4)
%bMat[14, 2] = d(q14', T6)
%bMat[14, 3] = d(q14', T7)
%bMat[14, 4] = d(q14', F)
%
%bMat[15, 1] = d(q15', T4)
%bMat[15, 2] = d(q15', T6)
%bMat[15, 3] = d(q15', T7)
%bMat[15, 4] = d(q15', F)
%
%bMat[16, 1] = d(q16', T4)
%bMat[16, 2] = d(q16', T6)
%bMat[16, 3] = d(q16', T7)
%bMat[16, 4] = d(q16', F)
%
%bMat[17, 1] = d(u4', T4)
%bMat[17, 2] = d(u4', T6)
%bMat[17, 3] = d(u4', T7)
%bMat[17, 4] = d(u4', F)
%
%bMat[18, 1] = d(u6', T4)
%bMat[18, 2] = d(u6', T6)
%bMat[18, 3] = d(u6', T7)
%bMat[18, 4] = d(u6', F)
%
%bMat[19, 1] = d(u7', T4)
%bMat[19, 2] = d(u7', T6)
%bMat[19, 3] = d(u7', T7)
%bMat[19, 4] = d(u7', F)
%
%cMat[1,1]=d(q1,q1)
%cMat[1,2]=d(q1,q2)
%cMat[1,3]=d(q1,q3)
%cMat[1,4]=d(q1,q4)
%cMat[1,5]=d(q1,q5)
%cMat[1,6]=d(q1,q6)
%cMat[1,7]=d(q1,q7)
%cMat[1,8]=d(q1,q8)
%cMat[1,9]=d(q1,u4)
%cMat[1,10]=d(q1,u6)
%cMat[1,11]=d(q1,u7)
%
%cMat[2,1]=d(q2,q1)
%cMat[2,2]=d(q2,q2)
%cMat[2,3]=d(q2,q3)
%cMat[2,4]=d(q2,q4)
%cMat[2,5]=d(q2,q5)
%cMat[2,6]=d(q2,q6)
%cMat[2,7]=d(q2,q7)
%cMat[2,8]=d(q2,q8)
%cMat[2,9]=d(q2,u4)
%cMat[2,10]=d(q2,u6)
%cMat[2,11]=d(q2,u7)
%
%cMat[3,1]=d(q3,q1)
%cMat[3,2]=d(q3,q2)
%cMat[3,3]=d(q3,q3)
%cMat[3,4]=d(q3,q4)
%cMat[3,5]=d(q3,q5)
%cMat[3,6]=d(q3,q6)
%cMat[3,7]=d(q3,q7)
%cMat[3,8]=d(q3,q8)
%cMat[3,9]=d(q3,u4)
%cMat[3,10]=d(q3,u6)
%cMat[3,11]=d(q3,u7)
%
%cMat[4,1]=d(q4,q1)
%cMat[4,2]=d(q4,q2)
%cMat[4,3]=d(q4,q3)
%cMat[4,4]=d(q4,q4)
%cMat[4,5]=d(q4,q5)
%cMat[4,6]=d(q4,q6)
%cMat[4,7]=d(q4,q7)
%cMat[4,8]=d(q4,q8)
%cMat[4,9]=d(q4,u4)
%cMat[4,10]=d(q4,u6)
%cMat[4,11]=d(q4,u7)
%
%cMat[5,1]=d(q5,q1)
%cMat[5,2]=d(q5,q2)
%cMat[5,3]=d(q5,q3)
%cMat[5,4]=d(q5,q4)
%cMat[5,5]=d(q5,q5)
%cMat[5,6]=d(q5,q6)
%cMat[5,7]=d(q5,q7)
%cMat[5,8]=d(q5,q8)
%cMat[5,9]=d(q5,u4)
%cMat[5,10]=d(q5,u6)
%cMat[5,11]=d(q5,u7)
%
%cMat[6,1]=d(q6,q1)
%cMat[6,2]=d(q6,q2)
%cMat[6,3]=d(q6,q3)
%cMat[6,4]=d(q6,q4)
%cMat[6,5]=d(q6,q5)
%cMat[6,6]=d(q6,q6)
%cMat[6,7]=d(q6,q7)
%cMat[6,8]=d(q6,q8)
%cMat[6,9]=d(q6,u4)
%cMat[6,10]=d(q6,u6)
%cMat[6,11]=d(q6,u7)
%
%cMat[7,1]=d(q7,q1)
%cMat[7,2]=d(q7,q2)
%cMat[7,3]=d(q7,q3)
%cMat[7,4]=d(q7,q4)
%cMat[7,5]=d(q7,q5)
%cMat[7,6]=d(q7,q6)
%cMat[7,7]=d(q7,q7)
%cMat[7,8]=d(q7,q8)
%cMat[7,9]=d(q7,u4)
%cMat[7,10]=d(q7,u6)
%cMat[7,11]=d(q7,u7)
%
%cMat[8,1]=d(q8,q1)
%cMat[8,2]=d(q8,q2)
%cMat[8,3]=d(q8,q3)
%cMat[8,4]=d(q8,q4)
%cMat[8,5]=d(q8,q5)
%cMat[8,6]=d(q8,q6)
%cMat[8,7]=d(q8,q7)
%cMat[8,8]=d(q8,q8)
%cMat[8,9]=d(q8,u4)
%cMat[8,10]=d(q8,u6)
%cMat[8,11]=d(q8,u7)
%
%cMat[9,1]=d(u1,q1)
%cMat[9,2]=d(u1,q2)
%cMat[9,3]=d(u1,q3)
%cMat[9,4]=d(u1,q4)
%cMat[9,5]=d(u1,q5)
%cMat[9,6]=d(u1,q6)
%cMat[9,7]=d(u1,q7)
%cMat[9,8]=d(u1,q8)
%cMat[9,9]=d(u1,u4)
%cMat[9,10]=d(u1,u6)
%cMat[9,11]=d(u1,u7)
%
%cMat[10,1]=d(u2,q1)
%cMat[10,2]=d(u2,q2)
%cMat[10,3]=d(u2,q3)
%cMat[10,4]=d(u2,q4)
%cMat[10,5]=d(u2,q5)
%cMat[10,6]=d(u2,q6)
%cMat[10,7]=d(u2,q7)
%cMat[10,8]=d(u2,q8)
%cMat[10,9]=d(u2,u4)
%cMat[10,10]=d(u2,u6)
%cMat[10,11]=d(u2,u7)
%
%cMat[11,1]=d(u3,q1)
%cMat[11,2]=d(u3,q2)
%cMat[11,3]=d(u3,q3)
%cMat[11,4]=d(u3,q4)
%cMat[11,5]=d(u3,q5)
%cMat[11,6]=d(u3,q6)
%cMat[11,7]=d(u3,q7)
%cMat[11,8]=d(u3,q8)
%cMat[11,9]=d(u3,u4)
%cMat[11,10]=d(u3,u6)
%cMat[11,11]=d(u3,u7)
%
%cMat[12,1]=d(u4,q1)
%cMat[12,2]=d(u4,q2)
%cMat[12,3]=d(u4,q3)
%cMat[12,4]=d(u4,q4)
%cMat[12,5]=d(u4,q5)
%cMat[12,6]=d(u4,q6)
%cMat[12,7]=d(u4,q7)
%cMat[12,8]=d(u4,q8)
%cMat[12,9]=d(u4,u4)
%cMat[12,10]=d(u4,u6)
%cMat[12,11]=d(u4,u7)
%
%cMat[13,1]=d(u5,q1)
%cMat[13,2]=d(u5,q2)
%cMat[13,3]=d(u5,q3)
%cMat[13,4]=d(u5,q4)
%cMat[13,5]=d(u5,q5)
%cMat[13,6]=d(u5,q6)
%cMat[13,7]=d(u5,q7)
%cMat[13,8]=d(u5,q8)
%cMat[13,9]=d(u5,u4)
%cMat[13,10]=d(u5,u6)
%cMat[13,11]=d(u5,u7)
%
%cMat[14,1]=d(u6,q1)
%cMat[14,2]=d(u6,q2)
%cMat[14,3]=d(u6,q3)
%cMat[14,4]=d(u6,q4)
%cMat[14,5]=d(u6,q5)
%cMat[14,6]=d(u6,q6)
%cMat[14,7]=d(u6,q7)
%cMat[14,8]=d(u6,q8)
%cMat[14,9]=d(u6,u4)
%cMat[14,10]=d(u6,u6)
%cMat[14,11]=d(u6,u7)
%
%cMat[15,1]=d(u7,q1)
%cMat[15,2]=d(u7,q2)
%cMat[15,3]=d(u7,q3)
%cMat[15,4]=d(u7,q4)
%cMat[15,5]=d(u7,q5)
%cMat[15,6]=d(u7,q6)
%cMat[15,7]=d(u7,q7)
%cMat[15,8]=d(u7,q8)
%cMat[15,9]=d(u7,u4)
%cMat[15,10]=d(u7,u6)
%cMat[15,11]=d(u7,u7)
%
%cMat[16,1]=d(u8,q1)
%cMat[16,2]=d(u8,q2)
%cMat[16,3]=d(u8,q3)
%cMat[16,4]=d(u8,q4)
%cMat[16,5]=d(u8,q5)
%cMat[16,6]=d(u8,q6)
%cMat[16,7]=d(u8,q7)
%cMat[16,8]=d(u8,q8)
%cMat[16,9]=d(u8,u4)
%cMat[16,10]=d(u8,u6)
%cMat[16,11]=d(u8,u7)
%
%cMat[17,1]=d(q9,q1)
%cMat[17,2]=d(q9,q2)
%cMat[17,3]=d(q9,q3)
%cMat[17,4]=d(q9,q4)
%cMat[17,5]=d(q9,q5)
%cMat[17,6]=d(q9,q6)
%cMat[17,7]=d(q9,q7)
%cMat[17,8]=d(q9,q8)
%cMat[17,9]=d(q9,u4)
%cMat[17,10]=d(q9,u6)
%cMat[17,11]=d(q9,u7)
%
%cMat[18,1]=d(q10,q1)
%cMat[18,2]=d(q10,q2)
%cMat[18,3]=d(q10,q3)
%cMat[18,4]=d(q10,q4)
%cMat[18,5]=d(q10,q5)
%cMat[18,6]=d(q10,q6)
%cMat[18,7]=d(q10,q7)
%cMat[18,8]=d(q10,q8)
%cMat[18,9]=d(q10,u4)
%cMat[18,10]=d(q10,u6)
%cMat[18,11]=d(q10,u7)
%
%dMat[1,1]=d(q1,T4)
%dMat[1,2]=d(q1,T7)
%dMat[1,3]=d(q1,F)
%
%dMat[2,1]=d(q2,T4)
%dMat[2,2]=d(q2,T7)
%dMat[2,3]=d(q2,F)
%
%dMat[3,1]=d(q3,T4)
%dMat[3,2]=d(q3,T7)
%dMat[3,3]=d(q3,F)
%
%dMat[4,1]=d(q4,T4)
%dMat[4,2]=d(q4,T7)
%dMat[4,3]=d(q4,F)
%
%dMat[5,1]=d(q5,T4)
%dMat[5,2]=d(q5,T7)
%dMat[5,3]=d(q5,F)
%
%dMat[6,1]=d(q6,T4)
%dMat[6,2]=d(q6,T7)
%dMat[6,3]=d(q6,F)
%
%dMat[7,1]=d(q7,T4)
%dMat[7,2]=d(q7,T7)
%dMat[7,3]=d(q7,F)
%
%dMat[8,1]=d(q8,T4)
%dMat[8,2]=d(q8,T7)
%dMat[8,3]=d(q8,F)
%
%dMat[9,1]=d(u1,T4)
%dMat[9,2]=d(u1,T7)
%dMat[9,3]=d(u1,F)
%
%dMat[10,1]=d(u2,T4)
%dMat[10,2]=d(u2,T7)
%dMat[10,3]=d(u2,F)
%
%dMat[11,1]=d(u3,T4)
%dMat[11,2]=d(u3,T7)
%dMat[11,3]=d(u3,F)
%
%dMat[12,1]=d(u4,T4)
%dMat[12,2]=d(u4,T7)
%dMat[12,3]=d(u4,F)
%
%dMat[13,1]=d(u5,T4)
%dMat[13,2]=d(u5,T7)
%dMat[13,3]=d(u5,F)
%
%dMat[14,1]=d(u6,T4)
%dMat[14,2]=d(u6,T7)
%dMat[14,3]=d(u6,F)
%
%dMat[15,1]=d(u7,T4)
%dMat[15,2]=d(u7,T7)
%dMat[15,3]=d(u7,F)
%
%dMat[16,1]=d(u8,T4)
%dMat[16,2]=d(u8,T7)
%dMat[16,3]=d(u8,F)
%
%dMat[17,1]=d(q9,T4)
%dMat[17,2]=d(q9,T7)
%dMat[17,3]=d(q9,F)
%
%dMat[18,1]=d(q10,T4)
%dMat[18,2]=d(q10,T7)
%dMat[18,3]=d(q10,F)
%
%encode aMat,bMat%,cMat,dMat

unitsystem kg, meter, sec

output q1 m, q2 m, q3 rad, q4 rad, q5 rad, q6 rad, q7 rad, q8 rad
output q9 rad, q10 rad, q11 rad, q12 rad, q13 rad, q14 rad, q15 rad, q16 rad
output q17 m, q18 m
output u1 m/s, u2 m/s, u3 rad/s, u4 rad/s, u5 rad/s, u6 rad/s, u7 rad/s, u8 rad/s
output u9 rad/s, u10 rad/s, u11 rad/s, u12 rad/s, u13 rad/s, u14 rad/s, u15 rad/s, u16 rad/s

code dynamics() WhippleArmsLateralForceDynamics.c

%---------------------------------------------------------------------%
%         save output
%---------------------------------------------------------------------%

save WhippleArmsLateralForce.all

%---------------------------------------------------------------------%
