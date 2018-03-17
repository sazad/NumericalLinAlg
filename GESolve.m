% Write a MATLAB M-file that solves Ax=b, an n x n system by
% (i) computing the LU factorization of A
% (ii) solving Ly = b by forward substitution, and 
% (iii) solving Ux=y by back substitution.


function b = GESolve(A, b)
% This M-file solves the n x n system Ax=b using Gaussian elimination
% and with no additional memory.

% The matrix A is first modified into its LU factorization where
% U is the elements of the diagonal and above of A, and 
% L is the the elements under the diagonal of A without the diagonal of 1's

% The function finds the LU in terms of block matrix multiplication
% using the following facts (Chp. 4):
% U(1,1) = A(1,1), so do not modify A
% U(1,2) = A(1,2), so do not modify A
% L(1,2) = 0 and U(2,1) = 0, but is omitted from A, so do not modify A
% L(1,1) = 1 but is ommited from A, so do not modify A

for i = 1:size(A)
    % L(2,1) = A(2,1) * (A(1,1)*(A(1,1))^(-1))
    % So, change the lower elements of A into L using
    A(i+1:size(A),i) =  A(i+1:size(A), i) * ((A(i,i))^(-1)); 
    
    % L(2,2) * U(2,2) = A(2,2) - (L(2,1)U(1,2))
    % So, change submatrix A(2,2) to be L(2,2)U(2,2)
    A(i+1:size(A),i+1:size(A))=A(i+1:size(A), i+1:size(A))-(A(i+1:size(A), i)* A(i, i+1:size(A)));
    
    % Solve Ly=b using forward substitution and the lower triangular matrix
    % b(i) = b(i)/L(i,i); but in this case, L(i,i) = 1, so not necessary
    for j = i+1:size(A)
        b(j) = b(j) - (A(j,i) * b(i));
    end
end

% After this, we have solved Ly = b with y stored in b
% Next we solve Ux = y using backwards substitution
for i = size(A):-1:1
    b(i) = b(i)/A(i,i);
    for j = i-1:-1:1
        b(j) = b(j) - (A(j,i)*b(i));
    end 
end
end