% This M-File solves an n x n tridiagonal system Ax = b
% by the tridiagonal solver given parameters
%     b: n x 1 matrix 
%     l: n-1 x 1 matrix, subdiagonal of A
%     d: n x 1 matrix, diagonal of A
%     u: n-1 x 1 matrix, superdiagonal of A
% and returns the values
%     b: the solution x of Ax=b
%     l: subdiagonal of the L of the LU factorization of A
%     d: the superdiagonal of the U of the LU factorization of A
%     u: the diagonal of the U of the LU factorization of A
% without using any additional memory

function [b, l, d, u] = TriSolve(b,l,d,u)

u(1) = u(1)/d(1);
for i = 2:size(u)
    u(i) = u(i)/(d(i) - (l(i-1).*u(i-1)));
end

b(1) = b(1)/d(1);
for i = 2:size(b)
    b(i) = (b(i)-(l(i-1).*b(i-1)))/(d(i) - (l(i-1).*u(i-1)));
end

%Back Substituion:
%b(n) = b(n);
for i = (size(b)-1):-1:1
    b(i) = b(i) - (u(i).*b(i+1));
end
end
