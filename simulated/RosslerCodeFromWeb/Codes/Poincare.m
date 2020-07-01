% To plot Poincare map with initial conditions [1,1,0].
% a, b and c are parameter in the Rossler system
% x0 = [1 1 0]; initial condition
% 'x' is an array containing values of x,y, and z variables.
% 't' is the time variable
% RES = Resolution per cycle
% c = 1, 2.5, 3.5, 4, 4.15, 4.2, 5, 5.3, 5.7, 5.9, 5.99, 6
%--------------------------------------------------------------------------
clear; clc;
c = input('Enter the value of parameter c: ');
a = 0.2; b = 0.2;
f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)]; % Rossler system
x0 = [1; 1; 0]; 
tspan = 0:0.1:5000;
[t,x] = ode45(f,tspan,x0');
n = length(t); RES = 500; p = 0.5;
poincare_x = x(round(n*p):RES:n,1);
poincare_y = x(round(n*p):RES:n,2);
poincare_z = x(round(n*p):RES:n,3);
%--------------------------------------------------------------------------
figure(1);
plot(poincare_x, poincare_z,'r*'); xlabel('x(t)'); ylabel('z(t)');
title(sprintf('Poincare map for c = %.2f',c));
%--------------------------------------------------------------------------
grid on