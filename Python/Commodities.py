import numpy as np
from typing import Tuple, Optional
from numpy.typing import DTypeLike

def schwartz2factor(delta0:float,S0:float,r:float,sigma1:float,sigma2:float,kappa:float,alpha:float,l:float,t:np.ndarray,W:np.ndarray,
                    dtype:Optional[DTypeLike]=np.float64)\
    ->Tuple[np.ndarray,np.ndarray]:
    N=t.size
    dt=t[-1]/(N-1)

    W1=W[:,:,0]
    dW2=np.diff(W[:,:,1],1,axis=0) # shape=(time,simulations)

    "Ornstein-Uhlenbeck"
    expkt=np.exp(-kappa*t)
    expkt2=np.exp(kappa*t)
    delta = delta0*expkt+(alpha-l/kappa)*(1-expkt)+expkt*sigma2*stochInt(expkt2,dW2,dtype=dtype)

    "Geometrical Brownian Motion"
    I_delta_dt = lebInt(delta,dt,N,dtype=dtype)
    S=S0*np.exp((r-sigma1**2/2)*t-I_delta_dt+sigma1*W1)
    return S,delta


def stochInt(f:np.ndarray,dW:np.ndarray,
             dtype:Optional[DTypeLike]=np.float64)\
    -> np.ndarray:
    [N,M]=dW.shape
    I_f_dW = np.zeros((N+1,M),dtype=dtype)
    I_f_dW[1:,:]=np.cumsum(f[:-1].reshape(N,-1)*dW,axis=0)
    return I_f_dW

def lebInt(f:np.ndarray,dt:float,N:int,
             dtype:Optional[DTypeLike]=np.float64)\
    -> np.ndarray:
    [_,M]=f.shape
    I_dt = np.zeros((N,M),dtype=dtype)
    I_dt[1:,:]=np.cumsum(f[:-1,:]*dt,axis=0)
    return I_dt

if __name__=="__main__":
    from StochIntegrator import brownianMotion
    import time
    dtype=np.float32

    T=3.0
    N=3*24*10
    M=1000000
    rho=np.ones((2,2),dtype=dtype)
    rho[0,1]=.9
    rho[1,0]=.9

    tic = time.time()
    W=brownianMotion(T,N,M,d=2,rho=rho[0,1],dtype=dtype)
    toc = time.time()
    print(f'Elapsed time {toc-tic} s')
    t=np.linspace(0,T,N,endpoint=True,dtype=dtype).reshape(N,1)

    # d0,S0,r,sigma1,sigma2,kappa,alpha,lambda
    params=[.03,100,.03,1.0,1.0,1.0,1.0,0.5]

    tic = time.time()
    S,delta = schwartz2factor(*params,t,W,dtype=dtype)
    toc = time.time()
    print(f'Elapsed time {toc-tic} s')

    print(S.shape)
    print(S.dtype)