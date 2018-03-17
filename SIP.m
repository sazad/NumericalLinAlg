%This M-File: 
%  [e, v, r, k] = SIP(A, s, tol, M)
%     A: n x n matrix
%     s: scalar shift
%     tol: tolerance for relative error
%     M: maximum # of iterations performed
%solves using the shifted inverse power method for:
%     e: approx. eigenvalue
%     v: eigenvector corresponding to e
%     r: relative error of eigenvalue
%     k: # of iterations performed
    
function [e, v, r, k] = SIP(A, s, tol, M)
[row, col] = size(A);
Ashift = A - s*eye(row);
%Intial guess x0 choosen randomly
x = rand(row,1);
%Solve (A-sI)x = z using LUP factorization of A-sI
[L, U, p] = lu(Ashift, 'vector');
v = x/norm(x); %z0
v = v(p);
%Forward sub: y = L \ vp(p);
for i = 1:row
    for j = i+1:row
        v(j) = v(j) - L(j,i)*v(i);
    end
end

%Backwards Sub:x = U \ y;
for i = row:-1:1
    v(i) = v(i)/U(i,i);
    for j = i-1:-1:1
        v(j) = v(j) - U(j,i)*v(i);
    end
end

gamma = transpose(v) * x; %gamma
e = s + (1/gamma); %approx eigenvalue
k = 0; %# of iterations
r = Inf; %error

%while((#iterations < max iterations) && (rel. error >= tolerance))
while ((k <= M) && (r >= tol))
    v = x/norm(x);
    x = v(p);
    %y = L \ v(p);
    for i = 1:row
        for j = i+1:row
            x(j) = x(j) - L(j,i)*x(i);
        end
    end
    
    %x = U \ y;
    for i = row:-1:1
        x(i) = x(i)/U(i,i);
        for j = i-1:-1:1
            x(j) = x(j) - U(j,i)*x(i);
        end
    end
    
    eold = e; %hold the old e for calculating error
    gamma = transpose(v) * x; %put new e_k into e
    e = s + (1/gamma);
    r = abs(eold - e)/abs(e);
    k = k + 1;
end
end