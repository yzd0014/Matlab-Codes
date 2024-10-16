syms A B [3 3] matrix
C = (A*B).'

%% 
clear;

% Import the Symbolic Math Toolbox
import sym.*

% Define the symbolic vectors
V1 = sym('v1', [3, 1]);  % First vector
V2 = sym('v2', [3, 1]);  % Second vector

% Compute the symbolic dot product
dotProduct = dot(V1, V2);
eye(3)
%% 
clear;
syms P [6 1] matrix
p1 = P(1:3, 1);
symmatrix2sym(p1)
%% 
clear;
syms A [1 3] matrix
syms B [3 2] matrix
syms X [2 1] matrix
t = A*sin(B*X);
Dt = diff(t,X);
symmatrix2sym(Dt)


%% 
clear;
a = [1;2;3];
b = [2,2,2];
c = a.*b;
d = a + c;
%% 
clear;
tspan = [0 5];
y0 = 0;
[t,y] = ode45(@(t,y) 2*t, tspan, y0);
%% 
clear;
syms w1 w2 w3
W = sym([0 -w3 w2;w3 0 -w1;-w2 w1 0]);
syms R [3 3]
O=W*R
