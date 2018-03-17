function b = BackSub(A,b)
%% Description:
% This M-file solves upper triangular systems using backward substition.
% For A:nxn upper triangular matrix and b:nx1 column vector,
% the function call solves the system Ax = b
%     b = BackSub(A,b)
% It first checks that the folloeing initial conditions for solving an
% Upper Triangular System are satisfied:
% 1.) A is a square matrix 
% 2.) b is a column vector
% 3.) A and b have the same number of columns 
% 4.) A is an upper triangular matrix
% Then, it uses backward substitution to overwrite b with solution x

[rowa, cola] = size(A);
[rowb, colb] = size(b);

%Check that initial conditions for solving Upper Tri. System are satisfied
if (rowa ~= cola) 
    msg = sprintf('Invalid Matrix A:%dx%d must be a square matrix',...
        rowa, cola);
elseif (colb ~= 1)
    msg = sprintf('Invalid Matrix b:%dx%d, must be a column vector',...
        rowb, colb);
elseif (rowa ~= rowb)
    msg = sprintf('Matrices A:%dx%d and b:%dx%d must have same number of rows,',...
        rowa, cola, rowb, colb);
elseif ~(checkUpperTriangular(A))
    disp('Matrix A is not an Upper Triangular Matrix')
else
    %Solve Ax=b using backwards substitution
    for j = rowa:-1:1
        b(j) = b(j)/A(j,j);
        for i = (j-1):-1:1
            b(i) = b(i) - (A(i,j) * b(j));
        end
    end
end


%% 