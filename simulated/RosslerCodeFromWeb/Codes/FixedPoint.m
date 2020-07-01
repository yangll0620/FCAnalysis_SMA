function FP = FixedPoint(a,b,c)
% This will calculate the fixed point of the system.
% There are two fixed point for Rossler system
A1 = c + sqrt(c^2 - 4 * a * b);
A2 = c - sqrt(c^2 - 4 * a * b); 
A3 = -c + sqrt(c^2 - 4 * a * b);
A4 = -c - sqrt(c^2 - 4 * a * b);
FP1 = [A1/2 A4/(2*a) A1/(2*a)];   % first fixed point
FP2 = [A2/2 A3/(2*a) A2/(2*a)];   % second fixed point
FP = [FP1; FP2];
end