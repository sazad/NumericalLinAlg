function[x, RRerr, Iter] = GS(A, b, x0, tol, M)
% This M-file solves the system Ax = b by the iterative Gauss-Seidel Method
% parameters:
%     A:      n x n matrix
%     b:      n x 1 vector
%     x0:     initial guess of solution x
%     tol:    tolerance for relative residual error
%     M:      max number of iterations allowed
% returns:
%     x:      approximate solution to Ax=b
%     RRerr:  relative residual error in the inf-norm
%     Iter:   # of iterations performed
%

%BEGIN FUNCTION:
%Intialize return values
x = x0;
RRerr = inf;
Iter = 0;
maxIter = M;
n = length(A);

% Iteration Formula: Mx^new = b + Nx^old
% Gauss-Seidel Method: A = M-N where M = L+D and N = -U
M = tril(A,0);
N = M-A;

while ((Iter < maxIter) && (RRerr > tol))
    x0 = x; %set x_old = previous approx. x found
    Iter = Iter + 1;
    R = b+N*x0; %RHS of iteration formula
    
    %Solve Mx^new = x using forward sub
    for i = 2:n
        x(i) = R(i) - (M(i, 1:i-1) * x(1:i-1));
        x(i) = x(i)/M(i,i);
    end
    % Calculate relative residul error using
    % infinity-norm of x = max(abs(x))
    RRerr = max(abs(x - x0))./max(abs(x));
end
end