import numpy as np
import matplotlib.pyplot as plt

def gs(wmin,wmax,wincrement,diag, diagsub,n, b,max_iter, tolerance):
    totalpts = (wmax-wmin)/wincrement
    res = np.zeros(int(totalpts))
    err = np.zeros(int(totalpts))
    tol = tolerance * np.ones(n)
    wvals = np.arange(wmin, wmax, wincrement)
    errNminus1 = 0
    winc = 0
    for w in wvals:
        xnow = np.zeros(n)
        xnext = xnow.copy()
        for k in range(0, max_iter+1):
            for i in range(n):
                if i == 0:
                    xnext[i] = (b - (diagsub * xnow[i+1])) / diag
                elif i == n - 1:
                    xnext[i] = (b - diagsub * ((1-w) * xnow[i-1] + w * xnext[i-1])) / diag
                else:
                    xnext[i] = (b - (diagsub * (1-w) * xnow[i-1]) - diagsub * (xnow[i+1] + w * xnext[i-1])) / diag
            # if np.all(np.abs(xnext - xnow) < tol):
               
            # if np.all(np.abs(xnext - xnow) < tol):
            #     print(f'solved in {k} iterations')
            #     return xnext, k, err
            xnow = xnext.copy()
            currentres = np.zeros(n)
            for j in range(n):
                if j == 0:
                    currentres[j] = b - 2*xnext[j] + xnext[j+1]
                elif j == n - 1:
                    currentres[j] = b + xnext[j-1] - 2*xnext[j]
                else:
                    currentres[j] = b + xnext[j-1] - 2*xnext[j] + xnext[j+1]

            if k == max_iter-1:
                errNminus1 = xnext.copy()
            if k == max_iter:
                res[winc] = np.linalg.norm(currentres)
                err[winc] = np.linalg.norm(xnext-errNminus1)
                wvals[winc] = w.copy()
                print("\n----------------------")
                print('       w:', w)
                print('\n||b-Ax5000||  = ',res[winc])
                print('||x5000-x4999|| = ',err[winc])
                print("----------------------")
                winc += 1
                
    return wvals, res, err






n = 500
main_diag = 2
sub_diag = -1 
b = 1
wmin= 1.5
wmax=2
wincrement = .02
tol = 0.00001

wvals, residuals, error = gs(wmin,wmax,wincrement,main_diag,sub_diag, n,b, 5000, tol)
print("residuals:", residuals)
print("error:", error)
print('w values:', wvals)
plt.subplot(1,2,1)
plt.plot( wvals, error, 'b-o', linewidth = 3)
plt.title('error vs omega')
plt.xlabel('w')
plt.ylabel('error: ||x5000-x4999||')

plt.subplot(1,2,2)
plt.plot( wvals, residuals, 'k-o', linewidth = 3)
plt.title('residuals vs omega')
plt.xlabel('w')
plt.ylabel('resdiuals: ||b-Ax5000||')
plt.show()