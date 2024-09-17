import numpy as np

def gs(A, b, w, max_iter, tolerance):
    tol = tolerance * np.ones(len(A))
    xnow = np.zeros(len(A))
    xnext = xnow.copy()
    
    for k in range(1, max_iter + 1):
        for i in range(len(A)):
            if i == 0:
                xnext[i] = (b[i] - (A[i, i+1] * xnow[i+1])) / A[i, i]
            elif i == len(A) - 1:
                xnext[i] = (b[i] - A[i, i-1] * ((1-w) * xnow[i-1] + w * xnext[i-1])) / A[i, i]
            else:
                xnext[i] = (b[i] - (A[i, i-1] * (1-w) * xnow[i-1]) - A[i, i+1] * (xnow[i+1] + w * xnext[i-1])) / A[i, i]
        
        if np.all(np.abs(xnext - xnow) < tol):
            print(f'solved in {k} iterations')
            return xnext, k
        
        xnow = xnext.copy()
    
    return xnext, max_iter


# Creating K500
n = 500
main_diag = 2 * np.ones(n)
sub_diag = -1 * np.ones(n - 1)
K_500 = np.diag(main_diag) + np.diag(sub_diag, 1) + np.diag(sub_diag, -1)
b = np.ones((n+1))
w = 1.989
tol = 0.00001


solution, iter = gs(K_500, b, w, 50000, tol)
# print("Solution:", solution)
print("Iterations:", iter)