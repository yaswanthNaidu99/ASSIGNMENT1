
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import shlex



N = 8   

pad = 2  
n = np.arange(N-pad)
fn = (-1/2)**n
hn1 = np.pad(fn, (0,pad), 'constant', constant_values=(0))
hn2 = np.pad(fn, (pad,0), 'constant', constant_values=(0))
h = hn1+hn2

xtemp = np.array([1.0,2.0,3.0,4.0,2.0,1.0])
x = np.pad(xtemp, (0,pad), 'constant', constant_values=(0))

def FFT_recur(x):
    N = x.shape[0] 
    
    if N ==2: 
        return np.array([[1,1],[1,-1]])@x  
    
    else:
        X_e = FFT_recur(x[0::2])
        X_o = FFT_recur(x[1::2])
        
        W = np.array([np.exp(-1j*2*np.pi*k/N) for k in range(int(N/2))])
        
        X = np.concatenate([X_e + W*X_o,X_e - W*X_o])  
        
    return X

def IFFT_recur(X):
    N= X.shape[0]
    x = (1/N)*FFT_recur(np.conj(X)) 
    return np.real(x)
    
print('x = ',np.round(x,3))
print('h = ',np.round(h,3))



X = FFT_recur(x)
H = FFT_recur(h)
Y = X*H
y = IFFT_recur(Y)

print('X[k] = ',np.round(X,3))
print('H[k] = ',np.round(H,3))
print('y = ',y)
plt.stem(range(0,N),y,use_line_collection = True)
plt.title('Filter Output using 8 point FFT')
plt.xlabel('$n$')
plt.ylabel('$y(n)$')
plt.grid()
plt.savefig('fft.eps')


