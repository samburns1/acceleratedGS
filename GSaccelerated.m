clc
clear
% A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
% b = [1;1;1;1];
% A\b;
% w = 1.5;
% [solution,iter] = gs(A, b, w, 1000)


n = 500;
main_diag = 4 * ones(n, 1);          % Main diagonal with 2s
sub_diag = -1 * ones(n - 1, 1);      % Sub and super diagonals with -1s

% Creating the tridiagonal matrix
K_500 = diag(main_diag) + diag(sub_diag, 1) + diag(sub_diag, -1);
b = [1:1:n]';
w = 1.4;
disp(K_500\b)
tol = 0.01;
[solution,iter] = gs(K_500, b, w, 5000, tol)


function [solution, iterations] = gs(A, b, w, iter, tolerance)
    tol = tolerance*(ones(length(A)));
    xnow = zeros(length(A),1);
    xnext = xnow;
    for k=1:iter
        for i=1:length(A)
            
            if i == 1
                xnext(i) = (b(i) - (A(i,i+1) * xnow(i+1))) / A(i,i);
            elseif i == length(A)
                xnext(i) = (b(i) - A(i,i-1) * ((1-w) * xnow(i-1) + w * xnext(i-1))) / A(i,i);
            else
                sum = (A(i,i-1) * (1-w) * xnow(i-1)) - A(i,i+1) * (xnow(i+1) + w * xnext(i-1));
                xnext(i) = (b(i) - sum)/ A(i,i);
            end
        end
        if abs(xnext-xnow) < tol
            solution = xnext;
            iterations = k;
            fprintf('solved in %.0f iterations', k)
            break
        end
        xnow = xnext;
    end
end