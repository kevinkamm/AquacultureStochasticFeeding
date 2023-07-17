import numpy as np
from numpy.linalg import solve
# from itertools import combinations
from typing import Optional, Callable
from numpy.typing import DTypeLike


def lsmc(exerciseValue:np.ndarray, X:np.ndarray,discount:float,
         costs:Optional[np.ndarray]=None,
         dtype:Optional[DTypeLike]=np.float64)\
        -> np.ndarray:
    [N,M]=exerciseValue.shape
    VC=np.empty((M,1))
    V=exerciseValue[-1,:].copy()
    exerciseRegion=[]
    exercise=(N-1)*np.ones((M,1),dtype=np.int64)
    # Glassermann p. 461
    for ti in range(N-1,0,-1): #[N-1,...,1]
        if costs is not None:
            VC=discount * basis(X[ti,:,:],V,dtype=dtype) -costs[ti]
            V=discount*V -costs[ti] # Longstaff-Schwartz
            #V=VC # Tsitsiklis and Van Roy
        else:
            VC= discount * basis(X[ti,:,:],V,dtype=dtype)
            V=discount*V # Longstaff-Schwartz
            #V=VC # Tsitsiklis and Van Roy
        ind = VC <= exerciseValue[ti,:]
        exercise[ind]=ti

        V[ind]=exerciseValue[ti,ind]
        exerciseRegion.append(X[ti,ind,:])
    
    return discount*np.mean(V),exercise,exerciseRegion


def basis(X:np.ndarray,
          Y:Optional[np.ndarray]=None,
          dtype:Optional[DTypeLike]=np.float64)\
    -> np.ndarray:
    sz = X.shape
    m=sz[-1]
    nd=np.prod(sz[:-1])
    x=X.reshape((nd,m))

    # C = list(combinations(range(m),2)) #slower than numpy
    C = np.stack(np.triu_indices(m, k=1), axis=-1)   
    a = np.ones((nd,1+2*m+len(C)),dtype=dtype)
    a[:,1:m+1]=x
    a[:,m+1:2*m+1]=x**2
    a[:,2*m+1:]=x[:,C[:,0]]*x[:,C[:,1]]
    if Y is not None:
        szY=Y.shape
        coeff = solve(a.T @ a,a.T @ Y.reshape((-1,1)))
        return (a@coeff).reshape(szY)
    else:
        return a.reshape(sz[:-1]+(a.shape[-1],))

def saveDecisions(t:np.ndarray,tau:np.ndarray,S:np.ndarray,
                  seed:Optional[int] =2,
                  saveas:Optional[str]=''):
    tau = tau.flatten()
    N = t.size
    M = tau.size
    exercise=[np.array([],ndmin=1,dtype=S.dtype) for _ in range(N)]
    cont=[np.array([],ndmin=1,dtype=S.dtype) for _ in range(N)]
    totalpaths=np.arange(M,dtype=np.int32)
    paths=np.arange(M,dtype=np.int32)
    for ti in range(N):  
        wi = tau == t[ti]
        paths = np.setdiff1d(paths,totalpaths[wi])
        cont[ti]=S[ti,paths,:]
        if len(wi)>0:
            exercise[ti]=S[ti,wi,:]
            continue
    exercise=np.asarray(exercise,dtype=object)
    cont=np.asarray(cont,dtype=object)
    if saveas != '':
        np.save(saveas+'_exercise_'+str(seed)+'.npy',exercise,allow_pickle=True)
        np.save(saveas+'_cont_'+str(seed)+'.npy',cont,allow_pickle=True)
    return exercise,cont


if __name__=="__main__":
    from StochIntegrator import brownianMotion
    from Commodities import schwartz2factor
    import time

    dtype=np.float32

    T=1.0
    N=100
    M=100000
    rho=np.ones((2,2),dtype=dtype)
    rho[0,1]=.9
    rho[1,0]=.9

    W=brownianMotion(T,N,M,d=2,rho=rho,dtype=dtype)
    t=np.linspace(0,T,N,endpoint=True,dtype=dtype).reshape(N,1)
    dt=T/(N-1)

    params=[.03,100,.03,1.0,1.0,1.0,1.0,0.5]

    S,delta = schwartz2factor(*params,t,W,dtype=dtype)
    X=np.stack([S,delta],axis=2)
    exerciseValue = S-5
    discount=np.exp(-params[2]*dt).astype(dtype)

    tic=time.time()
    V=lsmc(exerciseValue,X,discount,costs=np.ones((N,1),dtype=dtype)*7,dtype=dtype)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')
    print(V)
    # X=np.stack([S[-1,:].T,delta[-1,:].T],axis=1)
    # Y=S[-1,:].T**2+delta[-1,:].T**2
    # tic=time.time()
    # B=basis(X,Y)
    # toc=time.time()
    # print(f'Elapsed time {toc-tic} s')
    # print(np.mean((B-Y)**2))