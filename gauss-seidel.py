import numpy as np

def gs(diag, diagsub,n, w, b, max_iter, tolerance):
    tol = tolerance * np.ones(n)
    xnow = np.zeros(n)
    xnext = xnow.copy()
    
    for k in range(1, max_iter + 1):
        for i in range(n):
            if i == 0:
                xnext[i] = (b - (diagsub * xnow[i+1])) / diag
            elif i == n - 1:
                xnext[i] = (b - diagsub * ((1-w) * xnow[i-1] + w * xnext[i-1])) / diag
            else:
                xnext[i] = (b - (diagsub * (1-w) * xnow[i-1]) - diagsub * (xnow[i+1] + w * xnext[i-1])) / diag
        
        if np.all(np.abs(xnext - xnow) < tol):
            print(f'solved in {k} iterations')
            return xnext, k
        
        xnow = xnext.copy()
    
    return xnext, max_iter


# Creating K500
n = 500
main_diag = 2 
sub_diag = -1 
b = 1
w = 1.989
tol = 0.00001


solution, iter = gs(main_diag,sub_diag, n, w,b, 5000, tol)
# print("Solution:", solution)
print("Iterations:", iter)