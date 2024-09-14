clear;
syms q [3 1]

q_norm = norm(q);
a = sin(q_norm)/q_norm;
b = (1 - cos(q_norm))/(q_norm * q_norm);
c = (1 - a)/(q_norm * q_norm);

J = eye(3) + b * get_skew_symmetric(q) + c * get_skew_symmetric(q) * get_skew_symmetric(q);
det(J)
eqn = det(J) == 0;
S = solve(eqn, q);
real_q = [double(S.q1);double(S.q2);double(S.q3)];
norm(real_q/(2*pi))

%% 
clear;
syms q [3 1]

q_norm = norm(q);
a = sin(q_norm)/q_norm;
b = (1 - cos(q_norm))/(q_norm * q_norm);
c = (1 - a)/(q_norm * q_norm);

u11_local = [0; 1; 0];
Rotation = cos(q_norm)*eye(3)+a*get_skew_symmetric(q)+b*(q*q.'); 
u11_global = Rotation * u11_local;
Md = 

function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end


