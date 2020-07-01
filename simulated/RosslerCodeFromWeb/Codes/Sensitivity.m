% Sensitivity to initial condition
% a, b and c are fixed parameter in Rossler system
% x0 and x1 are two different initial condition
% t0 and tf are initial and final time and h is step size 
%--------------------------------------------------------------------------
clc; clear;
a = 0.2; b = 0.2; c = 5.7; t0 = 0; tf = 100; h = 0.05;
f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];
tspan = t0:h:tf;      % time interval for ode45()
%--------------------------------------------------------------------------
x0 = [1; 1; 0];    % first initial condition
[t,x] = ode45(f,tspan,x0);
plot(t,x(:,1),'-b*','MarkerSize',2); hold on
%--------------------------------------------------------------------------
x1 = [1; 1; 0.05];   % second initial condition
[t,x] = ode45(f,tspan,x1); 
plot(t,x(:,1),'-r','LineWidth',1);
hold all;
xlabel('t'); ylabel('x(t)');
legend(sprintf('x(t) with I.C [%.1f, %.1f, %.1f]',x0),...
    sprintf('x(t) with I.C [%.1f, %.1f, %.2f]',x1));
title('Sensitivity to Initial Condition');