% Helper function: valid = checkUpperTriangular(A)
% Given a matrix A, this functions returns true if A is an 
% upper triangular matrix and false if not
function valid = checkUpperTriangular(A)
[rowa,cola] = size(A);
if (rowa == cola)
    valid = 1;
    i = 2;
    while (valid) && (i <= rowa)
        j = 1;
        while (valid) && (j < i)
            if (A(i,j) ~= 0)
                valid = 0;
            end
            j = j + 1;
        end
        i = i + 1;
    end
end
