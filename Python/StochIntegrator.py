import numpy as np
from typing import Optional
from numpy.typing import DTypeLike

def brownianMotion(T:float,N:int,M:int,
                   d:Optional[int]=1,
                   rho:Optional[float]=None,
                   antithetic:Optional[bool]=False,
                   dtype:Optional[DTypeLike]=np.float64,
                   rng:Optional[np.random.Generator]=np.random.default_rng())\
    -> np.ndarray:

    dt=T/(N-1)
    dW=np.sqrt(dt)*rng.standard_normal((N-1,M,d),dtype=dtype)#4x faster than np.random.randn
    W=np.zeros((N,M,d),dtype=dtype)
    W[1:,:,:]=np.cumsum(dW,axis=0)

    if rho is not None:
        if d==2:
            if np.size(rho)==1:
                r=rho
            else:
                r=rho[0,1]
            W=np.stack([W[:,:,0],r*W[:,:,0]+np.sqrt(1-r**2)*W[:,:,1]],axis=2)
        else:
            L=np.linalg.cholesky(rho).reshape((1,1,d,d))
            W=np.squeeze(L@W.reshape((N,M,d,1)),axis=3)
    if antithetic:
        return np.concatenate([W,-W],axis=1)
    else:
        return W

if __name__=="__main__":
    import time
    dtype=np.float32

    T=3.0
    N=3*24*10
    M=100000
    rho=np.ones((2,2),dtype=dtype)
    rho[0,1]=.9
    rho[1,0]=.9
    d=rho.shape[0]

    rng = np.random.default_rng(0) 
    # tic=time.time()
    # W=brownianMotion(T,N,M,d=2,rho=rho[0,1],rng=rng,dtype=dtype)
    # toc=time.time()
    # print(f'Elapsed time {toc-tic} s')

    dt=T/(N-1)
    tic=time.time()
    dW=np.sqrt(dt)*rng.standard_normal((N-1,M,d),dtype=dtype)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')

    tic=time.time()
    W=np.zeros((N,M,d),dtype=dtype)
    W[1:,:,:]=np.cumsum(dW,axis=0)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')

    tic=time.time()
    dW=np.sqrt(dt)*rng.standard_normal((d,M,N-1),dtype=dtype)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')

    tic=time.time()
    W=np.zeros((d,M,N),dtype=dtype)
    W[:,:,1:]=np.cumsum(dW,axis=2)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')