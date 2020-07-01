% Projection of attractor on three plane for different parameter c
% a, b and c are parameter in Rossler system
% p = begining of steady state
% x0 = initial condition to solve Rossler system
% c = 1, 2.5, 3.5, 4, 4.15, 4.2, 5, 5.3, 5.7, 5.9, 5.99, 6
%--------------------------------------------------------------------------
clear; clc;
c = input('Enter the value of parameter c: ');
a = 0.2; b = 0.2; p = 0.5;
x0 = [1; 1; 0];
[t,x] = Rossler(a,b,c,x0);   % call the function to get the solution matrix 
n = length(t);
figure(1);
plot(x(1,round(n*p):n),x(2,round(n*p):n)); xlabel('x(t)'); ylabel('y(t)');
title(sprintf('Projection for c = %.2f',c));
figure(2);
plot(x(1,round(n*p):n),x(3,round(n*p):n)); xlabel('x(t)'); ylabel('z(t)');
title(sprintf('Projection on xz-plane for c = %.2f',c));
figure(3);
plot(x(2,round(n*p):n),x(3,round(n*p):n)); xlabel('y(t)'); ylabel('z(t)');
title(sprintf('Projection on yz-plane for c = %.2f',c));