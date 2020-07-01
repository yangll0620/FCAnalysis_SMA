function x = RK4(f,h,t,x)
K1 = h*f(t,x);
K2 = h*f(t + h/2, x + K1/2);
K3 = h*f(t + h/2, x + K2/2);
K4 = h*f(t + h, x + K3);
x = x + (K1 + 2 * K2 + 2 * K3 + K4)/4;
end