%% This program will generate the bifurcation map for parameter a by ploting
% b and c are fixed parameter in Rossler system while a is varying
% arange = range for the bifurcation parameter
% tspan = time interval for solving Rossler system by ode45()
% x0 = [1 1 0] initial condition for solving Rossler system
clear; clc;
a = 0.2; c = 5.7; global b; 
brange = 0.01:0.005:2;       % Range for parameter b
k = 0; tspan = 0:0.1:500;
xmax = [];                   % A matrix for storing the sorted value of x1
for b = brange 
    j = 1; k = k + 1;
    f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)]; % Rossler system
    x0 = [1 1 0];                    % initial condition for Rossler system
    k = k + 1; 
    [t,x] = ode45(f,tspan,x0); %calling ode45() to solve the Rossler system
    count = find(t>100); 
    x = x(count,:); 
    j = 1; 
    n = length(x(:,1));    % find the length of vector x1(x in our problem)
    for i=2 : n-1 
        % check for the min value in 1st column of sol matrix
        if (x(i-1,1)+eps) < x(i,1) && x(i,1) > (x(i+1,1)+eps)
            xmax(k,j)=x(i,1); % Sorting the values of x1 in increasing order
            j=j+1; 
        end 
    end 
    % generating bifurcation map by plotting j-1 element of kth row each time
    if j>1 
        plot(b,xmax(k,1:j-1),'k.','MarkerSize',1.5); 
    end 
    hold on; 
    index(k)=j-1; 
end 
xlabel('Bifurcation parameter b'); 
ylabel('x max'); 
title('Bifurcation diagram for b'); 