function x = EULER(f,h,t,x)
x = x + h*f(t,x);
end