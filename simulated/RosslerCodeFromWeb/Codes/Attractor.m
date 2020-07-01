%% To plot the Rossler Attractor with the fixed point of system
% x0 = [1;1;0] initial point (should be a column vector)
% a = 0.2; b = 0.2; c = 5.7; standard values 
% FP is the fixed point of system
% x is the solution matrix return by the function Rossler()
% t is the time variable
%--------------------------------------------------------------------------
function Attractor(a,b,c,x0)
[t,x] = Rossler(a,b,c,x0);     % the function Rossler() will solve the system
plot3(x(1,:),x(2,:),x(3,:),'-r'); hold all
FP = FixedPoint(a,b,c);     % the function will return the fixed point of the system
plot3(FP(1,1),FP(1,2),FP(1,3),'b*'); hold all   % plot of 1st fixed point
plot3(FP(2,1),FP(2,2),FP(2,3),'b*'); hold all   % plot of 2nd fixed point
xlabel('x(t)'); ylabel('y(t)'); zlabel('z(t)');
title(sprintf('Rossler Attractor with a = %.1f, b = %.1f & c = %.1f',a,b,c));
grid on;
end