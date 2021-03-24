import time
import numpy as np
import matplotlib.pyplot as plt

def DFT(x):
    N = x.shape[0]
    W_mat = np.zeros((N,N),dtype = 'complex_')
    X = np.zeros(x.shape,dtype = 'complex_')
    
    
    for i in range(N):
        for j in range(N):
            W_mat[i][j] = np.exp(-1j*(2*np.pi*i*j/N))
            
    X = W_mat@x
    return X


def FFT_recur(x):
    N = x.shape[0] 
    X = np.zeros(x.shape,dtype = 'complex_')
    
    if N ==2:
        return np.array([[1,1],[1,-1]])@x  
    
    else:
        X_e = FFT_recur(x[0::2])
        X_o = FFT_recur(x[1::2])
        
        W = np.array([np.exp(-1j*2*np.pi*k/N) for k in range(int(N/2))])
        
        X = np.concatenate([X_e + W*X_o,X_e - W*X_o])  
        
    return X

R = 9

t_DFT = np.zeros(R)
t_FFT = np.zeros(R)

for r in range(2,R+2):
    N = 2**r
    x = np.random.randint(1,10,N)
    t0 = time.time()
    X_dft = DFT(x)
    t1 = time.time()
    X_fft = FFT_recur(x)
    t2 = time.time()
    
    t_DFT[r-2] = t1 - t0
    t_FFT[r-2] = t2 - t1
    
R_vec = np.arange(2,R+2)
plt.xlabel('\u03B3')
plt.ylabel('Run time in secs')
plt.title('Run time comparisions for N-DFT and N-FFT  wrt \u03B3 \nwhere N = 2^\u03B3')
plt.plot(R_vec,t_DFT,label = 'DFT')
plt.plot(R_vec,t_FFT,label = 'FFT')
plt.grid()
plt.legend()  
plt.savefig('dft_fft.eps')  
