clear;
import sym.*

syms P [6 1]
r1 = P(1:3, 1);
r2 = P(4:6, 1);

syms P_dot [6 1]
r1_dot = P_dot(1:3, 1);
r2_dot = P_dot(4:6, 1);

theta1 = norm(r1);
theta2 = norm(r2);

u11_local = sym([-1; 1; 1]);
u12_local = sym([1; -1; -1]);
u22_local = u11_local;

%compute a b and c
a1 = sin(theta1)/theta1;
b1 = (1-cos(theta1))/(theta1*theta1);
c1 = (1-a1)/(theta1*theta1);
a1_dot = (c1-b1)*dot(r1, r1_dot);
b1_dot = (a1-2*b1)/(theta1*theta1)*dot(r1, r1_dot);
c1_dot = (b1-3*c1)/(theta1*theta1)*dot(r1, r1_dot);

a2 = sin(theta2)/theta2;
b2 = (1-cos(theta2))/(theta2*theta2);
c2 = (1-a2)/(theta2*theta2);
a2_dot = (c2-b2)*dot(r2, r2_dot);
b2_dot = (a2-2*b2)/(theta2*theta2)*dot(r2, r2_dot);
c2_dot = (b2-3*c2)/(theta2*theta2)*dot(r2, r2_dot);

%compute angular velocity
w1_local = a1*r1_dot-cross(b1*r1_dot, r1)+c1*dot(r1, r1_dot)*r1;
w2_local = a2*r2_dot-cross(b2*r2_dot, r2)+c2*dot(r2, r2_dot)*r2;
w1_global = w1_local;
w2_global = w1_local+w2_local;

%compute rotation matrix
skew_r1 = get_skew_symmetric(r1);
R1_local = cos(theta1)*eye(3)+a1*skew_r1+b1*(r1*r1.');
R1_global = R1_local;

skew_r2 = get_skew_symmetric(r2);
R2_local = cos(theta2)*eye(3)+a2*skew_r2+b2*(r2*r2.');
R2_global = R1_local * R2_local;

u11_global = R1_global * u11_local;
u12_global = R1_global * u12_local;
u22_global = R2_global * u22_local;

energy_old = compute_total_energy(P, P_dot);

%update state--------------------------------------------------
%compute H
J2 = a2*eye(3)+b2*get_skew_symmetric(r2)+c2*(r2*r2.');
H2_top = get_skew_symmetric(u22_global)*J2;
H2 = vertcat(H2_top, J2);

J1 = a1*eye(3)+b1*get_skew_symmetric(r1)+c1*(r1*r1.');
H1_top = get_skew_symmetric(u11_global)*J1;
H1 = vertcat(H1_top, J1);

%compute D
D2_left = vertcat(eye(3), sym(zeros(3)));
D2_topRight = get_skew_symmetric(u22_global-u12_global);
D2_right = vertcat(D2_topRight, eye(3));
D2 = horzcat(D2_left, D2_right);

%compute gamma
gamma_theta2 = (c2*dot(r2, r2_dot)+a2_dot)*r2_dot-cross(b2_dot*r2_dot, r2) +(c2_dot*dot(r2,r2_dot)+c2*dot(r2_dot, r2_dot))*r2;
gamma2_top = get_skew_symmetric(u22_global)*gamma_theta2 ...
-get_skew_symmetric(w2_global)*get_skew_symmetric(w2_global)*u22_global ...
+get_skew_symmetric(w1_global)*get_skew_symmetric(w1_global)*u12_global;
gamma2 = vertcat(gamma2_top, gamma_theta2);

gamma_theta1 = (c1*dot(r1, r1_dot)+a1_dot)*r1_dot-cross(b1_dot*r1_dot, r1) +(c1_dot*dot(r1, r1_dot)+c1*dot(r1_dot, r1_dot))*r1;
gamma1_top = get_skew_symmetric(u11_global)*gamma_theta1 ...
-get_skew_symmetric(w1_global)*get_skew_symmetric(w1_global)*u11_global;
gamma1 = vertcat(gamma1_top, gamma_theta1);

%compute gamma_t
gamma_t2 = D2*gamma1 + gamma2;
gamma_t1 = gamma1;

%compute Ht
H21 = D2*H1;
H22 = H2;
Ht2 = horzcat(H21,H22);

H11 = H1;
Ht1 = horzcat(H11, sym(zeros(6, 3)));

%compute Mr
M2 = eye(6);
M1 = eye(6);
Mr = Ht1.' * M1 * Ht1 + Ht2.' * M2 * Ht2;

%compute Qr
Fe = sym([0;-9.81;0;0;0;0]);
Qr = Ht1.'*(Fe+M1*gamma_t1)+Ht2.'*(Fe+M2*gamma_t2);

%update state
h = sym(0.001);
P_ddot = inv(Mr) * Qr;
P_dot_new = P_dot+P_ddot*h;
P_new = P+P_dot*h;

%energy_new = compute_total_energy(P_new, P_dot_new);

function energy = compute_total_energy(q, q_dot)
    r1 = q(1:3, 1);
    r2 = q(4:6, 1);
    r1_dot = q_dot(1:3, 1);
    r2_dot = q_dot(4:6, 1);

    theta1 = norm(r1);
    theta2 = norm(r2);

    u11_local = sym([-1; 1; 1]);
    u12_local = sym([1; -1; -1]);
    u22_local = u11_local;

    %compute a b and c
    a1 = sin(theta1)/theta1;    
    b1 = (1-cos(theta1))/(theta1*theta1);
    c1 = (1-a1)/(theta1*theta1);

    a2 = sin(theta2)/theta2;
    b2 = (1-cos(theta2))/(theta2*theta2);
    c2 = (1-a2)/(theta2*theta2);

    %compute angular velocity
    w1_local = a1*r1_dot-cross(b1*r1_dot, r1)+c1*dot(r1, r1_dot)*r1;
    w2_local = a2*r2_dot-cross(b2*r2_dot, r2)+c2*dot(r2, r2_dot)*r2;
    w1_global = w1_local;
    w2_global = w1_local+w2_local;

    %compute rotation matrix
    skew_r1 = get_skew_symmetric(r1);
    R1_local = cos(theta1)*eye(3)+a1*skew_r1+b1*(r1*r1.');  
    R1_global = R1_local;

    skew_r2 = get_skew_symmetric(r2);
    R2_local = cos(theta2)*eye(3)+a2*skew_r2+b2*(r2*r2.');
    R2_global = R1_local * R2_local;
    
    %compute position
    x1 = R1_global*-u11_local;
    x2 = R1_global*2*u12_local+R2_global*-u22_local;

    %compute velocity
    vel1 = cross(w1_local, 2*u12_local);
    vel2 = cross(w1_local, x2)+cross(w2_local, -u22_local);

    %kinatic energy
    T1 = 0.5*(w1_global).'*w1_global+0.5*(vel1).'*vel1;
    T2 = 0.5*(w2_global).'*w2_global+0.5*(vel2).'*vel2;

    %potential energy
    g = sym([0;-9.81;0]);
    V1 = g.'*x1;
    V2 = g.'*x2;

    %total energy before update
    energy = T1 + V1 + T2 + V2;
end

% symmatrix2sym
function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end