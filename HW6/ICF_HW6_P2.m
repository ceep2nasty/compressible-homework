p1 = 101.3 * 10^3 ;
p2 = 1.4 * 10 ^ 6 ;
gamma = 1.4 ;
R = 287 ;
T = 15+273; 
p2_p1 = p2/p1 ;

a1 = (gamma*R*T) ^0.5 ;

% Find M1


% Find W

W_sqrt = (gamma + 1) / (2 * gamma) * (p2_p1 - 1) + 1 ;

W = a1 * W_sqrt^0.5 ;

% find u_p

u_p = (a1/gamma) * (p2_p1 - 1) * (((2*gamma)/(gamma+1))/(p2_p1 + (gamma-1)/(gamma+1)))^0.5 ;