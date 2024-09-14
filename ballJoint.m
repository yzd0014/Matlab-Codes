clear;
import sym.*

syms P [6 1] matrix
x1 = P(1:3, 1);
x2 = P(4:6, 1);

syms P_dot [6 1] matrix
x1_dot = P_dot(1:3, 1);
x2_dot = P_dot(4:6, 1);

theta1 = norm(x1);
theta2 = norm(x2);

%compute a b and c
a1 = sin(theta1)/theta1;
b1 = (1-cos(theta1))/(theta1*theta1);
c1 = (1-a1)/(theta1*theta1);
a1_dot = (c1-b1)*dot(symmatrix2sym(x1), symmatrix2sym(x1_dot));
b1_dot = (a1-2*b1)/(theta1*theta1)*dot(symmatrix2sym(x1), symmatrix2sym(x1_dot));
c1_dot = (b1-3*c1)/(theta1*theta1)*dot(symmatrix2sym(x1), symmatrix2sym(x1_dot));

a2 = sin(theta2)/theta2;
b2 = (1-cos(theta2))/(theta2*theta2);
c2 = (1-a2)/(theta2*theta2);
a2_dot = (c2-b2)*dot(symmatrix2sym(x2), symmatrix2sym(x2_dot));
b2_dot = (a2-2*b2)/(theta2*theta2)*dot(symmatrix2sym(x2), symmatrix2sym(x2_dot));
c2_dot = (b2-3*c2)/(theta2*theta2)*dot(symmatrix2sym(x2), symmatrix2sym(x2_dot));

%get glocal rotation matrices for each body
skew_x1 = get_skew_symmetric(x1);
R1_local = cos(theta1)*eye(3)+a1*skew_x1+b1*(x1*x1.');
R1_global = R1_local;

skew_x2 = get_skew_symmetric(x2);
R2_local = cos(theta2)*eye(3)+a2*skew_x2+b2*(x2*x2.');
R2_global = R1_local * R2_local;

%update angular velocity for each body 
w1_local = a1*x1_dot-cross(symmatrix2sym(b1*x1_dot), symmatrix2sym(x1))+c1*dot(symmatrix2sym(x1), symmatrix2sym(x1_dot))*x1;
w2_local = a2*x2_dot-cross(symmatrix2sym(b2*x2_dot), symmatrix2sym(x2))+c2*dot(symmatrix2sym(x2), symmatrix2sym(x2_dot))*x2;
w1_global = w1_local;
w2_global = w1_local + w2_local;

%compute H
u11_local = sym([-1; 1; 1]);
u12_local = sym([1; -1; -1]);
u22_local = u11_local;

u11_global = R1_global * u11_local;
u12_global = R1_global * u12_local;
u22_global = R2_global * u22_local;

J2 = a2*eye(3)+b2*get_skew_symmetric(x2)+c2*(x2*x2.');
H2_top = get_skew_symmetric(symmatrix2sym(u22_global))*J2;
H2 = vertcat(H2_top, J2);

J1 = a1*eye(3)+b1*get_skew_symmetric(x1)+c1*(x1*x1.');
H1_top = get_skew_symmetric(symmatrix2sym(u11_global))*J1;
H1 = vertcat(H1_top, J1);

%compute D
D2_left = vertcat(eye(3), sym(zeros(3)));
D2_topRight = get_skew_symmetric(symmatrix2sym(u22_global)-symmatrix2sym(u12_global));
D2_right = vertcat(D2_topRight, eye(3));
D2 = horzcat(D2_left, D2_right);

%compute gamma
gamma_theta2 = (c2*dot(symmatrix2sym(x2), symmatrix2sym(x2_dot))+a2_dot)*x2_dot-cross(symmatrix2sym(b2_dot*x2_dot), symmatrix2sym(x2)) ...
+(c2_dot*dot(symmatrix2sym(x2),symmatrix2sym(x2_dot))+c2*dot(symmatrix2sym(x2_dot), symmatrix2sym(x2_dot)))*x2;
gamma2_top = get_skew_symmetric(symmatrix2sym(u22_global))*gamma_theta2 ...
-get_skew_symmetric(symmatrix2sym(w2_global))*get_skew_symmetric(symmatrix2sym(w2_global))*u22_global ...
+get_skew_symmetric(symmatrix2sym(w1_global))*get_skew_symmetric(symmatrix2sym(w1_global))*u12_global;
gamma2 = vertcat(gamma2_top, gamma_theta2);

gamma_theta1 = (c1*dot(symmatrix2sym(x1), symmatrix2sym(x1_dot))+a1_dot)*x1_dot-cross(symmatrix2sym(b1_dot*x1_dot), symmatrix2sym(x1)) ...
+(c1_dot*dot(symmatrix2sym(x1),symmatrix2sym(x1_dot))+c1*dot(symmatrix2sym(x1_dot), symmatrix2sym(x1_dot)))*x1;
gamma1_top = get_skew_symmetric(symmatrix2sym(u11_global))*gamma_theta1 ...
-get_skew_symmetric(symmatrix2sym(w1_global))*get_skew_symmetric(symmatrix2sym(w1_global))*u11_global;
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
I2_local = [8/12 0 0;
      0 8/12 0;
      0 0 8/12];     
I1_local = I2_local;
I2_global = R2_global*sym(I2_local)*R2_global.';
I1_global = R1_global*sym(I1_local)*R1_global.';
M2_topLeft = eye(3);
M2_bottomLeft = sym(zeros(3));
M2_topRight = sym(zeros(3));
M2_bottomRight = I2_global;
M2_top = horzcat(M2_topLeft, M2_bottomRight);
M2_bottom = horzcat(M2_bottomLeft, M2_bottomRight);
M2 = vertcat(M2_top, M2_bottom); 

M1_topLeft = eye(3);
M1_bottomLeft = sym(zeros(3));
M1_topRight = sym(zeros(3));
M1_bottomRight = I1_global;
M1_top = horzcat(M1_topLeft, M1_bottomRight);
M1_bottom = horzcat(M1_bottomLeft, M1_bottomRight);
M1 = vertcat(M1_top, M1_bottom);

Mr = Ht1.' * M1 * Ht1 + Ht2.' * M2 * Ht2;

%compute Qr
Fe = sym([0;-9.81;0;0;0;0]);
Fv_top = sym(zeros(3,1));
Fv1_bottom = cross(symmatrix2sym(-w1_global), symmatrix2sym(I1_global*w1_global));
Fv1 = vertcat(Fv_top, Fv1_bottom);
Fv2_bottom = cross(symmatrix2sym(-w2_global), symmatrix2sym(I2_global*w2_global));
Fv2 = vertcat(Fv_top, Fv2_bottom);

Qr = Ht1.' * (Fe + Fv1 + M1 * gamma_t1) + Ht2.' * (Fe + Fv2 + M2 * gamma_t2);

%compute first order derivative
h = sym(0.001);
F = inv(Mr) * Qr * h;
Dp = diff(symmatrix2sym(F), P(1,1))

% symmatrix2sym
function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end