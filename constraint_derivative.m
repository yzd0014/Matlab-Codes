clear;
import sym.*

syms r [3 1]
syms p [3 1]
theta = norm(r);
alpha = sin(theta)/theta;
beta = (1-cos(theta))/(theta*theta);
Rot = cos(theta)*eye(3)+alpha*get_skew_symmetric(r)+beta*r*r.';
s = cross(Rot*p,p);
Dp = diff(s, r(1,1));
Dp

function M = get_skew_symmetric(v)
   M = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end