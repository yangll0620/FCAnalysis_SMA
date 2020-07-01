function [t,x] = Rossler(a,b,c,x0)
% x is denoted as x(1), y as x(2) and z as x(3)
% x0 = [1; 1; 0];    Initial point
% a = 0.2; b = 0.2; c = 5.7;
% dt = The step size of t
% tmax = maximum time in second allowed to run the program in each step h
% k = number of solution excluding initial point
% x = Solution matrix 
% FP = Fixed point of the system
% Order = 1 for Euler, 4 for RK4 and 5 for RK5
Order = 5;                 % Change here for the method
tic;
f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];  % Rossler system
k = 5000;       % k+1 is the number of columns in the solution matrix M
dt = 0.05; time = 0; 
tmax = 20; p = 0;
x(:, 1) = x0 (:, 1); t(1) = 0;
while (time<tmax)
    n = 2^p;
    h = dt/n;
    t0 = 0;
    x0 = [1; 1; 0];
    C = clock;
    for j = 1:n*k
        if Order == 1
            x0 = Euler(f,h,t0,x0);         % calling the function EULER()
        else if Order == 4
                x0 = RK4(f,h,t0,x0);       % calling the function RK4()
            else x0 = RK5(f,h,t0,x0);      % calling the function RK5()
            end
        end
        t0 = j*h;
            x(:,j+1) = x0;           % saving the values 
            t(j+1) = t0;             
    end
    time = etime(clock,C);
    p = p + 1;
end
end