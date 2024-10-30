clear;
import sym.*
syms alpha beta
A=sym([1 0 -sin(beta);0 cos(alpha) sin(alpha)*cos(beta);0 -sin(alpha) cos(alpha)*cos(beta)]);
A
D=inv(A);
D

%% 
clear;
syms alpha beta gamma
X=sym([1 0 0;0 cos(alpha) -sin(alpha);0 sin(alpha) cos(alpha)]);
Y=sym([cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)]);
Z=sym([cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0;0 0 1]);
X
Y
Z
T=Z*Y*X
%% 
clear;
% Euler derivative for my convention
import sym.*
syms alpha beta
A=sym([0 sin(alpha) cos(alpha)*cos(beta);1 0 sin(beta);0 cos(alpha) -sin(alpha)*cos(beta)]);
A
D=inv(A);
D
%% 
clear;
syms alpha beta gamma
X=sym([1 0 0;0 cos(alpha) sin(alpha);0 -sin(alpha) cos(alpha)]);
Y=sym([cos(beta) 0 -sin(beta);0 1 0;sin(beta) 0 cos(beta)]);
Z=sym([cos(gamma) sin(gamma) 0;-sin(gamma) cos(gamma) 0;0 0 1]);
T=X*Y*Z
%% 
clear;
%convert angular velocity to Eulder dot for my convention
alpha=0;
beta=0;
gamma=0;
M=[-cos(alpha)*sin(beta)/cos(beta) 1 sin(alpha)*sin(beta)/cos(beta);
    sin(alpha) 0 cos(alpha);
    cos(alpha)/cos(beta) 0 -sin(alpha)/cos(beta)];

w=[-1;1;0];
w=w/norm(w);
w=-w*2;
qdot=M*w;
qdot

