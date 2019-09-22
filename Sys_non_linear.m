% Author - Ambarish Prashant Chandurkar
% Non-linear systems using Newton Raphson Method (for three variables)
clc;
close all;
clear all;
syms x y z;
f=x^3 + y^3 - 53; % Enter the first function here
fx=diff(f,x);fy=diff(f,y);fz=diff(f,z); %Calculation of partial derivatives

g=2*y^3 + z^4 - 69;% Enter the second function here
gx=diff(g,x);gy=diff(g,y);gz=diff(g,z); %Calculation of partial derivatives

h=3*x^5 + 10*z^2 - 770;% Enter the third function here
hx=diff(h,x);hy=diff(h,y);hz=diff(h,z); %Calculation of partial derivatives

n=input('Enter the number of decimal places:');
epsilon = 5*10^-(n+1)

x0 = input('Enter the x intial approximation:');
y0 = input('Enter the y intial approximation:');
z0 = input('Enter the z intial approximation:');

xk = x0;yk = y0;zk = z0; %The intial approximations

for i=1:5 %Upto 5 iterations are considered
fk=vpa(subs(f,[x y z],[xk yk zk])); % Calculation of function values at each iteration
gk=vpa(subs(g,[x y z],[xk yk zk]));
hk=vpa(subs(h,[x y z],[xk yk zk]));

fxk=vpa(subs(fx,[x y z],[xk yk zk])); %Evalution of partial derivatives at each iteration
fyk=vpa(subs(fy,[x y z],[xk yk zk]));
fzk=vpa(subs(fz,[x y z],[xk yk zk]));

gxk=vpa(subs(gx,[x y z],[xk yk zk]));
gyk=vpa(subs(gy,[x y z],[xk yk zk]));
gzk=vpa(subs(gz,[x y z],[xk yk zk]));

hxk=vpa(subs(hx,[x y z],[xk yk zk]));
hyk=vpa(subs(hy,[x y z],[xk yk zk]));
hzk=vpa(subs(hz,[x y z],[xk yk zk]));

xyzk = [xk yk zk]'; %The old values of x,y,z
J = [fxk fyk fzk;gxk gyk gzk;hxk hyk hzk]; % The jacobian matrix
Fk = [fk gk hk]'; % Matrix of function values

if abs(max(eig(inv(J)))) >= 1 %Sufficient Condition of Convergence
    display('Cannot find solution !');
    exit();
end

xyzkplus1 = (xyzk - inv(J)*Fk); %The Iteration formula

if abs(xyzkplus1 - xyzk) < epsilon %Stopping Criterion for iterations
    break;
end
end

xyzans = xyzkplus1-rem(xyzkplus1,10^-n); %Displaying results upto required accuracy
fprintf('The Root Matrix is :\n');
display(xyzans);

