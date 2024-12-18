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
%% 
A = [ 10^-20 1; 1 1 ];
[L1,U1]=lu_cs210(A);
[L2,U2]=lu(A);
A = [ -1 2; 4 4 ];
max(abs(A(:,1)))
%% 
clear;
%% 
clear;
A=[-1 0 1 0 0 1;1 1 0 -1 0 0;0 -1 -1 0 1 0;0 0 0 1 -1 -1];
x=[1 0 1 1 1 0];
x=x';
A*x
%% 
clear;
import sym.*
syms a b c
E=[1 0 0;-a 1 0;a*c-b -c 1];
A=[1 0 0;a 1 0;b c 1];
O=E*A;
O
%% 
clear;
A=[2 1 0;6 4 2;0 3 5];
E2=[1 0 0;-3 1 0;0 0 1];
E2*A
%% 
n = 100;
A = diag(rand(n,1));
A(1,:) = rand(1,n);
A(:,1) = rand(1,n);
spy(A);
[L,U]=my_lu(A);
[L,U]=lu(A)
spy(L);
%spy(U);
%% 
A=[ 10^-20 1; 1 1 ];
[L, U]=lu(A);
L
U
function [L,U] = lu_cs210(A)
    n = size(A,1);
    L = zeros(size(A));
    U = zeros(size(A));
    A2 = A;
    for k = 1:n
        if A2(k,k) == 0
            'Encountered 0 pivot. Stopping'
            return
        end
        for i = 1:n
            L(i,k) = A2(i,k)/A2(k,k);
            U(k,i) = A2(k,i);
        end
        for i = 1:n
            for j = 1:n
                A2(i,j) = A2(i,j) - L(i,k)*U(k,j);
            end
        end
    end
end

function [L,U] = my_lu(A)
n = size(A,1);
L = zeros(size(A));
A2 = A;
for k = 1:n
if A2(k,k) == 0
'Encountered 0 pivot. Stopping'
return
end
L(k,k) = 1;
for i = k+1:n
L(i,k) = A2(i,k)/A2(k,k);
end
for i = k+1:n
for j = k+1:n
A2(i,j) = A2(i,j) - L(i,k)*A2(k,j);
end
end
end
U = triu(A2);
end



