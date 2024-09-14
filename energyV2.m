clear;
dt = linspace(0.0000001,0.01,500);

q0 = linspace(0,pi,100);
q1 = linspace(0,pi,100);
q2 = linspace(0,pi,100);
qs = zeros(100*100*100, 3);
n = 1;
for i=1:1:100
    for j=1:1:100
        for k = 1:1:100
            qs(n, 1) = q0(i);
            qs(n, 2) = q1(j);
            qs(n, 3) = q2(k);
            n = n + 1;
        end
    end
end
q = qs(20,:)';
q_norm = norm(q);

qdot0 = linspace(0,pi,100);
qdot1 = linspace(0,pi,100);
qdot2 = linspace(0,pi,100);
qdots = zeros(100*100*100, 3);
n = 1;
for i=1:1:100
    for j=1:1:100
        for k = 1:1:100
            qdots(n, 1) = q0(i);
            qdots(n, 2) = q1(j);
            qdots(n, 3) = q2(k);
            n = n + 1;
        end
    end
end
% qdot = qdots(1,:);
qdot = [-2;10;0];

%compute a, b, c
if q_norm < 0.0001
    a = 1 - q_norm * q_norm/6;
    b = 0.5 - q_norm * q_norm/24;
    c = 1/6 - q_norm * q_norm/120;
else
    a = sin(q_norm)/q_norm;
    b = (1 - cos(q_norm))/(q_norm * q_norm);
    c = (1 - a)/(q_norm * q_norm);
end

%compute derivative of a, b, c
a_dot = (c-b)*dot(q,qdot);
if q_norm < 0.0001
    b_dot = (-1/12+1/180*q_norm*q_norm)*dot(q, qdot);
    c_dot = (-1/60+1/1260*q_norm*q_norm)*dot(q, qdot);
else
    b_dot = (a-2*b)/(q_norm*q_norm)*dot(q, qdot);
    c_dot = (b-3*c)/(q_norm*q_norm)*dot(q, qdot);
end

%compute rotation matrix
u11_local = [0; 1; 0];
Rotation = cos(q_norm)*eye(3)+a*get_skew_symmetric(q)+b*(q*q.'); 
u11_global = Rotation * u11_local;

%compute Jacobian matrix
J_bottom = eye(3) + b * get_skew_symmetric(q) + c * get_skew_symmetric(q) * get_skew_symmetric(q);
J_top = get_skew_symmetric(u11_global)*J_bottom;
J = vertcat(J_top, J_bottom);

%compute angular velocity
w = J_bottom*qdot;

%compute Jdot_qdot
Jdot_qdot_bottom = (c*dot(q, qdot)+a_dot)*qdot-cross(b_dot*qdot, q)+(c_dot*dot(q, qdot)+c*dot(qdot, qdot))*q;
Jdot_qdot_top = get_skew_symmetric(u11_global)*Jdot_qdot_bottom-get_skew_symmetric(w)*get_skew_symmetric(w)*u11_global;
Jdot_qdot = vertcat(Jdot_qdot_top, Jdot_qdot_bottom);

%compute Mr
Mr = J.'*J;

%compute Qr
Qr = -J.'*Jdot_qdot;

%update
qddot = Mr\Qr;
qdot_new = qdot+qddot*dt;
q_new = q+qdot_new.*dt;

%get energy
E = compute_total_energy_1(q, qdot);
sz = size(dt);
sz = sz(2);
E_new = zeros(1,sz);
for i = 1:sz
    E_new(i) = compute_total_energy_1(q_new(:,i), qdot_new(:,i));
end

%plot
energy_ratio = E_new/E;
plot(dt, energy_ratio);
xlabel('dt');
ylabel('energy ratio');
%% 
clear;


function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end

function E = compute_total_energy_1(q, qdot)
    q_norm = norm(q);    

    %compute a, b, c
    if q_norm < 0.0001
        a = 1 - q_norm * q_norm/6;
        b = 0.5 - q_norm * q_norm/24;
        c = 1/6 - q_norm * q_norm/120;
    else
        a = sin(q_norm)/q_norm;
        b = (1 - cos(q_norm))/(q_norm * q_norm);
        c = (1 - a)/(q_norm * q_norm);
    end
    
    u11_local = [-1; 1; 1];
    R_loccal = cos(q_norm)*eye(3)+a*get_skew_symmetric(q)+b*(q*q.'); 
    R_global = R_loccal;
    u11_global = R_global * u11_local;

    J_bottom = eye(3) + b * get_skew_symmetric(q) + c * get_skew_symmetric(q) * get_skew_symmetric(q);
    J_top = get_skew_symmetric(u11_global)*J_bottom;
    J = vertcat(J_top, J_bottom);

    V = J*qdot;
    E = V.'*V;
end



