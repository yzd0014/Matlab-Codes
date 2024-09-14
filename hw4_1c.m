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

dynamics = @(m1,m2,l1,l2,t,q) [q(3);q(4);M(m1,m2,l1,l2,q)\(-C(m1,m2,l1,l2,q)*q(3:4,1)-gr(m1,m2,l1,l2,q))];
dt = 1/1000;
T = 10;
t = 0:dt:(T-dt);

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
%q0 = [3*pi/4;0;0;0];
q0 = [0;0.5;2.1;1.2];
[t, q] = ode45(@(t,q) dynamics(m1,m2,l1,l2,t,q), t, q0);

% plot angular momentum
sz = size(t);
sz = sz(1);
H = zeros(3, sz);
AuglarMomentum = zeros(sz,1);
for i=1:sz
    v1 = zeros(3,1);
    v1(1,1) = l1*cos(q(i,1))*q(i,3);
    v1(2,1) = l2*sin(q(i,1))*q(i,3);
    r1 = zeros(3,1);
    r1(1,1) = l1*sin(q(i,1));
    r1(2,1) = -l1*cos(q(i,1));
    
    v2 = zeros(3,1);
    v2(1,1) = l1*cos(q(i,1))*q(i,3)+l2*cos(q(i,1)+q(i,2))*(q(i,3)+q(i,4));
    v2(2,1) = l1*sin(q(i,1))*q(i,3)+l2*sin(q(i,1)+q(i,2))*(q(i,3)+q(i,4));
    r2 = zeros(3,1);
    r2(1,1) = l1*sin(q(i,1))+l2*sin(q(i,1)+q(i,2));
    r2(2,1) = -l1*cos(q(i,1))-l2*cos(q(i,1)+q(i,2));
    
    H(:,i) = cross(r1,v1)+cross(r2,v2);
    AuglarMomentum(i) = 3*q(i,3)+q(i,4)+(2*q(i,3)+q(i,4))*cos(q(i,2));
end
plot(t,AuglarMomentum);
%animation
itr = 1;
x2_ic1 = zeros(2,1);
y2_ic1 = zeros(2,1);
figure; hold all;
for i=1:10:numel(t)
  clf; hold all;
  q1 = q(i,1);
  q2 = q(i,2);
  plotPendulum(q1,q2,l1,l2);
  x2_ic1(itr) = l1*sin(q1)+l2*sin(q1+q2);
  y2_ic1(itr) = -l1*cos(q1)-l2*cos(q1+q2);
  itr = itr+1;
  pause(0.0001);
  axis square;
end

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