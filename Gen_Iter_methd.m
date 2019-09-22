% Author - Ambarish Prashant Chandurkar
%The General Iteration Method (The Fixed Point Iteration Method)
clc;
close all;
clear all;
syms x;
f=3*x^3 + 4*x^2 + 4*x +1; %Enter the Function here
g=diff(f); %The Derivative of the Function
n=input('Enter the number of decimal places:');
epsilon = 5*10^-(n+1)
x0 = input('Enter the intial approximation:');
i=1;
%Now we determine 'Interval' needed for finding iteration function I(x)
%We know that condition is:
% |1+alpha*f'(x)| < 1 , from this we find 'alpha'
%This means 1+alpha*f'(x) < 1 and 1+alpha*f'(x) > -1
%So this means alpha -> [0,(-2/f'(x))]
rng(0,'twister'); %Initialsing the random number generator
alpha = -2/vpa(subs(g,x,x0)); %This is the starting value for alpha
x_current = x0;
while i~=200
    phi_of_x_at_x_current = x_current + (alpha*vpa(subs(f,x,x_current))); %Our iteration formula
 
err=abs(phi_of_x_at_x_current-x_current); %finding the error at each iteration
if  abs(1+alpha*vpa(subs(g,x,x_current)))>=1 %If the condition that |1+alpha*f'(x)| < 1 is NOT fulfilled at
    alpha = -1*(2/vpa(subs(g,x,x0))*rand); %any x_current, then we need to use differnt alpha !
    i=1; %So we randonmly select any other alpha ONLY from our interval and set i to 0 to start the iteration again.
elseif err<epsilon %checking the amount of error at each iteration
break %if enough accuracy has been achieved the stop the iterations
end
x_current =  phi_of_x_at_x_current;
i=i+1;
end

phi_of_x_at_x_current = phi_of_x_at_x_current - rem(phi_of_x_at_x_current,10^-n); %Displaying upto required decimal places
fprintf('The Root is : %f \n',phi_of_x_at_x_current);
