clear;
import sym.*

syms R Ry Rz [3 3]
syms theta_y theta_z c_angle
Rz=eye(3);
Rz(1,1)=cos(theta_z);
Rz(1,2)=-sin(theta_z);
Rz(2,1)=sin(theta_z);
Rz(2,2)=cos(theta_z);
Ry=eye(3);
Ry(1,1)=cos(theta_y);
Ry(1,3)=sin(theta_y);
Ry(3,1)=-sin(theta_y);
Ry(3,3)=cos(theta_y);
R=Rz*Ry;

syms x y z [3 1]
x=zeros(3,1);
y=zeros(3,1);
z=zeros(3,1);
x(2)=-1;
y(3)=1;
z(1)=-1;
R*x

syms s [3 1]
s=-get_skew_symmetric(y)*R*x;

syms J J0 J1 J2 [1 3]
J0=x.'*R.'*get_skew_symmetric(get_skew_symmetric(y)*R*z);
J1=-x.'*R.'*get_skew_symmetric(y)*get_skew_symmetric(R*z);
J2=-cos(c_angle)/norm(s)*s.'*get_skew_symmetric(y)*get_skew_symmetric(R*x);
J=J0+J1+J2;
J
%% 
clear;
% import sym.*
% 
% syms R Rx [3 3]
% syms theta_x c_angle
Rx=eye(3);
Rx(2,2)=cos(-0.25*pi);
Rx(2,3)=-sin(-0.25*pi);
Rx(3,2)=sin(-0.25*pi);
Rx(3,3)=cos(-0.25*pi);
R=Rx;

syms x y z [3 1]
x=zeros(3,1);
y=zeros(3,1);
z=zeros(3,1);
x(2)=-1;
y(3)=1;
z(1)=-1;


% syms s [3 1]
s=-get_skew_symmetric(y)*R*x;

% syms J J0 J1 J2 [1 3]
J0=x.'*R.'*get_skew_symmetric(get_skew_symmetric(y)*R*z);
J1=-x.'*R.'*get_skew_symmetric(y)*get_skew_symmetric(R*z);
J2=-cos(0.25*pi)/norm(s)*s.'*get_skew_symmetric(y)*get_skew_symmetric(R*x);
J=J0+J1+J2;
J
function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end
