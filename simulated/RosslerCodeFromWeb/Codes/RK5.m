%% this functio is called by Rossler() to solve the system
function x = RK5(f,h,t,x)
K1 = h*f(t, x);
K2 = h*f(t + h/4, x + K1/4);
K3 = h*f(t + h/4, x + (K1 + K2)/8);
K4 = h*f(t + h/2, x + (-K2 + 2*K3)/2);
K5 = h*f(t + 3*h/4, x + (-3*K1 + 9*K4)/16);
K6 = h*f(t + h, x + (-3*K1 + 2*K2 + 12*K3 - 12*K4 + 8*K5)/7);
x = x + (7*K1 + 32*K3 + 12*K4 + 32*K5 + 7*K6)/90;
end