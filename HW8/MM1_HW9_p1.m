syms lam real
A = [1 0 1 4;
     1 0 2 -1;
     2 1 3 -2];
b = [2; 1; -3];

xlam = simplify((A.'*A + lam*eye(4)) \ (A.'*b));  

x_min_norm = simplify(limit(xlam, lam, 0, 'right')) ;
