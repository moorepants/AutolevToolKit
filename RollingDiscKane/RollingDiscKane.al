autoz off
autorhs off
factoring off
newtonian n
frames a, b
points cn
variables x', y', psi', phi', theta'
motionvariables' wx',wy', wz'
bodies c
constants m, g, r, I, J
mass c=m
inertia c(b), I, J, I, 0, 0, 0

simprot(n, a, 3, psi)
simprot(a, b, 1, phi)
simprot(b, c, 2, theta)

w_c_n> = wx*b1> + wy*b2> + wz*b3>

wcn_coords> = psi'*a3> + phi'*b1> + theta'*b2>
kd_sys = [dot(wcn_coords> - w_c_n>, b1>), &
          dot(wcn_coords> - w_c_n>, b2>), &
          dot(wcn_coords> - w_c_n>, b3>)]
solve(kd_sys, [psi', phi', theta'])

w_a_n> = rhs(psi')*a3>
w_b_a> = rhs(phi')*b1>
w_c_b> = rhs(theta')*b2>

p_cn_co> = -r*b3>
v_co_n> = cross(w_c_n>, p_cn_co>)

vco1> = x'*n1> + y'*n2> + cross(psi'*a3> + phi'*b1>, -r*b3>)
vco2> = cross(psi'*a3> + phi'*b1> + theta'*b2>, -r*b3>)
eq1 = dot(vco1> - vco2>, n1>)
eq2 = dot(vco1> - vco2>, n2>)
solve([eq1, eq2], [x', y'])

alf_c_n> = dt(w_c_n>, n)
a_co_n> = dt(v_co_n>, n)

gravity(g*a3>)
zero = fr() + frstar()

% Energy
ke = dot(v_co_n>, v_co_n>)/2.0 + dot(w_c_n>, dot(I_C_CO>>, w_c_n>))/2.0
pe = m*g*dot(p_cn_co>, a3>)

% Momentum
p> = m*v_co_n>
H_C_CO> = dot(I_C_CO>>, w_c_n>)

unitsystem kg,m,s
output zero[1], zero[2], zero[3]

code algebraic() rollingdisc_kane.c

