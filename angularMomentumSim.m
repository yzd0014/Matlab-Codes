clear;
g = 0;
m11 = @(m1,m2,l1,l2,q) (m1+m2)*l1*l1+m2*l2*l2+2*m2*l1*l2*cos(q(2));
m12 = @(m1,m2,l1,l2,q) m2*l2*l2+m2*l1*l2*cos(q(2));
m21 = @(m1,m2,l1,l2,q) m2*l2*l2+m2*l1*l2*cos(q(2));
m22 = @(m1,m2,l1,l2,q) m2*l2*l2;

c11 = @(m1,m2,l1,l2,q) -2*m2*l1*l2*q(4)*sin(q(2));
c12 = @(m1,m2,l1,l2,q) -m2*l1*l2*q(4)*sin(q(2));
c21 = @(m1,m2,l1,l2,q) -m2*l1*l2*q(4)*sin(q(2))+m2*l1*l2*q(3)*sin(q(2));
c22 = @(m1,m2,l1,l2,q) m2*l1*l2*q(3)*sin(q(2));

M = @(m1,m2,l1,l2,q) [m11(m1,m2,l1,l2,q), m12(m1,m2,l1,l2,q); m21(m1,m2,l1,l2,q), m22(m1,m2,l1,l2,q)];
C = @(m1,m2,l1,l2,q) [c11(m1,m2,l1,l2,q), c12(m1,m2,l1,l2,q); c21(m1,m2,l1,l2,q), c22(m1,m2,l1,l2,q)];
gr = @(m1,m2,l1,l2,q) [m1*g*l1*sin(q(1))+m2*g*(l1*sin(q(1))+l2*sin(q(1)+q(2))); m2*g*l2*sin(q(1)+q(2))];
J = @(l1,l2,q) [l1*cos(q(1,1)),0;
    l1*sin(q(1,1)),0;
    l1*cos(q(1,1))+l2*cos(q(1,1)+q(2,1)),l2*cos(q(1,1)+q(2,1));
    l2*sin(q(1,1))+l2*sin(q(1,1)+q(2,1)),l2*sin(q(1,1)+q(2,1))];
J2 = @(l1,l2,q) [l1*cos(q(1,1)),0;
    l1*sin(q(1,1)),0;
    0,0;
    l1*cos(q(1,1))+l2*cos(q(1,1)+q(2,1)),l2*cos(q(1,1)+q(2,1));
    l2*sin(q(1,1))+l2*sin(q(1,1)+q(2,1)),l2*sin(q(1,1)+q(2,1));
    0,0];

dt = 1/1000;
T = 10;
t = 0:dt:(T-dt);
sz = size(t);
sz = sz(2);

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
%q0 = [3*pi/4;0;0;0];
q0 = [0;0;0;0];

S1 = [eye(3),zeros(3)];
S2 = [zeros(3),eye(3)];
B = [0,0,1];

q = zeros(2,sz);
qdot = zeros(2,sz);
q(:,1) = q0(1:2);
qdot(:,1) = q0(3:4);
L0 = get_angular_momentum(J,m1,m2,l1,l2,q0(1:2,1),q0(3:4,1));
Ls_0 = zeros(3,sz);
Ls_0(:,1) = L0;
Ls_1 = Ls_0;
%figure; hold all;
for i=1:sz-1
    %dynamics update
    mq = [q(1,i);q(2,i);qdot(1,i);qdot(2,i)];
    q_new = zeros(2,1);
    qdot_new = zeros(2,1);
    q_acc = M(m1,m2,l1,l2,mq)\(-C(m1,m2,l1,l2,mq)*mq(3:4,1)-gr(m1,m2,l1,l2,mq));
    qdot(:,i+1) = qdot(:,i)+dt*q_acc;
    q(:,i+1) = q(:,i)+dt*qdot(:,i+1);
    q_new = [q(1,i+1);q(2,i+1)];
    qdot_new = [qdot(1,i+1);qdot(2,i+1)];

    %angular momentum correction
    r = get_momentum_arm(l1,l2,q_new);
    r1 = [r(1,1);r(2,1);0];
    r2 = [r(3,1);r(4,1);0];
    A = B*(m1*get_skew_symmetric(r1)*S1+m2*get_skew_symmetric(r2)*S2)*J2(l1,l2,q_new);
    if abs(det(A*A.')) > 0.0000001
        lambda = inv(A*A.')*(-A*qdot_new+L0(3));
        qdot(:,i+1) = qdot(:,i+1)+A.'*lambda;
    end
    
    %angular momentum correction
    % Jc = [3+2*cos(q_new(2)),1+cos(q_new(2))];
    % A = Jc*Jc.';
    % if abs(det(A)) > 0.0000001
    %     lambda = inv(A)*(-Jc*qdot_new+L0(3));
    %     qdot(:,i+1) = qdot(:,i+1)+Jc.'*lambda;
    % end
    
    %compute angular momentum
    Ls_0(:,i+1) = get_angular_momentum(J,m1,m2,l1,l2,q_new,qdot_new);
    
    %visualization
    % clf; hold all;
    % plotPendulum(q(1,i+1),q(2,i+1),l1,l2);
    % pause(0.0000001);
    % axis square;
end

dynamics = @(m1,m2,l1,l2,t,q) [q(3);q(4);M(m1,m2,l1,l2,q)\(-C(m1,m2,l1,l2,q)*q(3:4,1)-gr(m1,m2,l1,l2,q))];
[t_out, q_out] = ode45(@(t,q) dynamics(m1,m2,l1,l2,t,q), t, q0);
subplot(1,3,1);
plot(t,q(1,:),t_out,q_out(:,1),LineWidth=2);
legend('explicit','ode45');
subplot(1,3,2);
plot(t,q(2,:),t_out,q_out(:,2),LineWidth=2);
legend('explicit','ode45');
subplot(1,3,3);
plot(t,Ls_0(3,:),LineWidth=2);

function plotPendulum(q1,q2,l1,l2)
  x1 = l1*sin(q1);
  y1 = -l1*cos(q1);
  x2 = l1*sin(q1)+l2*sin(q1+q2);
  y2 = -l1*cos(q1)-l2*cos(q1+q2);
  
  plot([0,x1],[0,y1],'k','linewidth',8);
  scatter(x1,y1,500,'k','filled','o');
  plot([x1,x2],[y1,y2],'k','linewidth',8);
  scatter(x2,y2,500,'k','filled','o');
  xlim([-2,2]);
  ylim([-2,2]);
  axis off;
end
function r = get_momentum_arm(l1,l2,q)
    r = zeros(4,1);
    
    r(1,1) = l1*sin(q(1,1));
    r(2,1) = -l1*cos(q(1,1));

    r(3,1) = l1*sin(q(1,1))+l2*sin(q(1,1)+q(2,1));
    r(4,1) = -l1*cos(q(1,1))-l2*cos(q(1,1)+q(2,1));
end
function L = get_angular_momentum(J,m1,m2,l1,l2,q,qdot)
    v = J(l1,l2,q)*qdot;
    v1 = [v(1,1);v(2,1);0];
    v2 = [v(3,1);v(4,1);0];

    r = get_momentum_arm(l1,l2,q);
    r1 = [r(1,1),r(2,1),0];
    r2 = [r(3,1),r(4,1),0];

    L = m1*cross(r1,v1)+m2*cross(r2,v2);
    L = L.';
end
function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end