import numpy as np
import time
from pathlib import Path

from StochIntegrator import brownianMotion
from Commodities import schwartz2factor
from Regression import lsmc,saveDecisions
from DeepDecision import *

from typing import List

"Auxiliary Functions"
def coarse(X,steps):
    sz=X.shape
    return np.concatenate([X[0].reshape((1,)+sz[1:]),X[steps-1::steps]],axis=0)

class fishFarm():
    def __init__(self,
                 salmonParam:List[float],
                 soyParam:List[float],
                 r:float = 0.0303, # interest rate
                 T:float=3.0, # time horizon
                 N:int=3*24, # time steps for LSMC
                 simFactor:int=10,
                 m:float=.1, # mortality rate
                 cr:float=1.1, # conversion rate
                 n0:int=10000, # number of recruits
                 hc:float=3, # variable harvesting cost per kg (NOK/kg)
                 fc:float=10, # variable feeding cost per kg per year (NOK/kg)
                 wInf:float=6, # asymptotic weight (kg)
                 a:float=1.113,
                 b:float=1.097,
                 c:float=1.43,
                 gamma:float=0.0, # utility
                 dtype:DTypeLike=np.float32,
                 verbose:Optional[int]=1,
                 trainBoundary:Optional[bool]=False):
        self.salmonParam = salmonParam
        self.soyParam = soyParam
        self.T = T
        self.N = N
        self.simFactor = simFactor
        self.m=m
        self.cr=cr
        self.n0 = n0
        self.hc = hc
        self.fc = fc
        self.wInf = wInf
        self.a = a
        self.b = b
        self.c = c
        self.gamma = gamma
        self.r = r
        self.dtype=dtype
        self.model = '_'.join([f'{x:1.2f}' for x in salmonParam])+'-'+'_'.join([f'{x:1.2f}' for x in soyParam])
        self.verbose=verbose
        self.trainBoundary=trainBoundary

    def generateFishFarm(self,seed:int,M:int):
        M=int(M)
        np.random.seed(seed)
        rng = np.random.default_rng(seed)
        t=np.linspace(0,self.T,self.N*self.simFactor,endpoint=True,dtype=self.dtype).reshape(-1,1)
        dt=self.T/(self.N*self.simFactor-1)

        W=brownianMotion(self.T,int(self.N*self.simFactor),M,d=2,rho=self.salmonParam[-3],dtype=self.dtype,rng=rng)
        salmonP,salmonDelta = schwartz2factor(self.salmonParam[-2],self.salmonParam[-1],self.r,*self.salmonParam[1:-3],t,W,dtype=self.dtype)


        W=brownianMotion(self.T,int(self.N*self.simFactor),M,d=2,rho=self.soyParam[-3],dtype=self.dtype,rng=rng)
        soyP,soyDelta = schwartz2factor(self.soyParam[-2],1.0,self.r,*self.soyParam[1:-3],t,W,dtype=self.dtype)
        # soyP=soyP/self.soyParam[-1]

        "Joint dynamics with reduced time steps"
        S=np.stack([salmonP,salmonDelta,soyP,soyDelta],axis=2)

        "Fish Farm Functions"
        wt=self.wInf*(self.a-self.b*np.exp(-self.c*t))**3 # Bertalanffy’s growth function
        dwt=np.diff(wt,1,axis=0)/dt
        nt=self.n0*np.exp(-self.m*t) # number of fish
        Xt=nt*wt # total biomass (kg)
        CH=Xt*self.hc # harvesting cost
        CF=np.zeros((self.N*self.simFactor,1))
        CF[1:]=dwt*nt[1:]*self.cr*self.fc

        CFstoch=(CF*soyP)**(1-self.gamma)/(1-self.gamma)
        I_CF_dt_stoch = np.zeros_like(CFstoch)
        I_CF_dt_stoch[1:] = np.cumsum(np.exp(-self.r*t[1:])*CFstoch[1:]*dt,axis=0)

        CFdeterm=(CF*np.mean(soyP,axis=1,keepdims=True))**(1-self.gamma)/(1-self.gamma)
        I_CF_dt_determ = np.zeros_like(CFdeterm)
        I_CF_dt_determ[1:] = np.cumsum(np.exp(-self.r*t[1:])*CFdeterm[1:]*dt,axis=0)

        "sim points -> eval points"
        t = coarse(t,self.simFactor)
        dt = self.T/(self.N-1)
        S = coarse(S,self.simFactor)
        Xt = coarse(Xt,self.simFactor)
        CH = coarse(CH,self.simFactor)
        I_CF_dt_stoch = coarse(I_CF_dt_stoch,self.simFactor)
        I_CF_dt_determ = coarse(I_CF_dt_determ,self.simFactor)
        CFstoch = coarse(CFstoch,self.simFactor)
        CFdeterm = coarse(CFdeterm,self.simFactor)

        "LSMC global"
        VH=(S[:,:,0]*Xt-CH)**(1-self.gamma)/(1-self.gamma) # Harvesting costs
        discount=np.exp(-self.r*dt)

        costs_stoch = CFstoch*dt
        obj_stoch = np.exp(-self.r*t)*VH-I_CF_dt_stoch

        costs_determ = CFdeterm*dt
        obj_determ = np.exp(-self.r*t)*VH-I_CF_dt_determ

        
        return t,dt,S,Xt,CH,VH,discount,costs_stoch,costs_determ,obj_stoch,obj_determ

    def exerciseLSMC(self,seed:int,M:int):
        if self.verbose > 0:
            print('Generate Fish Farm')
        tic = time.time()
        t,dt,S,Xt,CH,VH,discount,costs_stoch,costs_determ,obj_stoch,obj_determ = self.generateFishFarm(seed,M)
        toc = time.time()
        if self.verbose > 0:
            print(f'\t Elapsed time {toc-tic:1.2f} s.')
            print('Compute LSMC stoch:')

        tic = time.time()
        V_stoch,exercise_stoch,exerciseRegion_stoch=lsmc(VH,S,discount,costs=costs_stoch,dtype=self.dtype)
        toc = time.time()
        tau_stoch=t[exercise_stoch]

        if self.verbose > 0:
            print(f'\t Elapsed time {toc-tic:1.2f} s.')
            print('Compute LSMC determ:')

        tic = time.time()
        V_determ,exercise_determ,exerciseRegion_determ=lsmc(VH,S[:,:,:2],discount,costs=costs_determ,dtype=self.dtype)
        toc = time.time()
        tau_determ=t[exercise_determ]

        if self.verbose > 0:
            print(f'\t Elapsed time {toc-tic:1.2f} s.')

        return t,tau_stoch,tau_determ,V_stoch,V_determ,obj_stoch,obj_determ,S 


    def generateTrainingSet(self,M:Optional[int]=int(1e6),savedir:Optional[str]=''):
        path = Path().cwd()/Path(savedir)/Path('Trainingset')/Path(self.model)
        file = path/Path('LSMC_stoch_cont_'+'0'+'.npy')
        if not path.exists() or not file.exists():
            if self.verbose>0:
                print('Generate Training Set')
            t,tau_stoch,tau_determ,V_stoch,V_determ,obj_stoch,obj_determ,S =self.exerciseLSMC(0,M)
            if savedir != '':
                path.mkdir(parents=True,exist_ok=True)
                LSMC_stoch_exercise_train,LSMC_stoch_cont_train = saveDecisions(t,tau_stoch,S,seed=0,saveas=str(path)+'/LSMC_stoch')
                LSMC_determ_exercise_train,LSMC_determ_cont_train = saveDecisions(t,tau_determ,S[:,:,:2],seed=0,saveas=str(path)+'/LSMC_determ')
        else:
            if self.verbose>0:
                print('Load Training Set')
            LSMC_stoch_cont_train = np.load(str(path)+'/LSMC_stoch_cont_'+'0'+'.npy',allow_pickle=True)
            LSMC_stoch_exercise_train = np.load(str(path)+'/LSMC_stoch_exercise_'+'0'+'.npy',allow_pickle=True)
            LSMC_determ_cont_train = np.load(str(path)+'/LSMC_determ_cont_'+'0'+'.npy',allow_pickle=True)
            LSMC_determ_exercise_train = np.load(str(path)+'/LSMC_determ_exercise_'+'0'+'.npy',allow_pickle=True)

        return LSMC_stoch_cont_train,LSMC_stoch_exercise_train,LSMC_determ_cont_train,LSMC_determ_exercise_train

    def generateValidationSet(self,M:Optional[int]=int(1e5),seed:Optional[int]=2,savedir:Optional[str]=''):
        if self.verbose>0:
            print('Generate Validation Set')
        t,tau_stoch,tau_determ,V_stoch,V_determ,obj_stoch,obj_determ,S =self.exerciseLSMC(0,M)
        LSMC_stoch_exercise_val,LSMC_stoch_cont_val = saveDecisions(t,tau_stoch,S,seed=seed,saveas='')
        LSMC_determ_exercise_val,LSMC_determ_cont_val = saveDecisions(t,tau_determ,S[:,:,:2],seed=seed,saveas='')

        return LSMC_stoch_cont_val,LSMC_stoch_exercise_val,LSMC_determ_cont_val,LSMC_determ_exercise_val,t,tau_stoch,tau_determ,V_stoch,V_determ,obj_stoch,obj_determ,S 
    
    def compareStoppingTimes(self,M:Optional[int]=1e5,seed:Optional[int]=2,savedir:Optional[str]=''):
        if self.trainBoundary:
            LSMC_stoch_cont_train,LSMC_stoch_exercise_train,LSMC_determ_cont_train,LSMC_determ_exercise_train = self.generateTrainingSet(savedir=savedir)
        LSMC_stoch_cont_val,LSMC_stoch_exercise_val,LSMC_determ_cont_val,LSMC_determ_exercise_val,t,tau_stoch,tau_determ,V_stoch,V_determ,obj_stoch,obj_determ,S  = self.generateValidationSet(savedir=savedir)
        
        print(f'Stoch LSMC:Fish pond value {V_stoch} at mean stopping time {np.mean(tau_stoch)}')
        print(f'Determ LSMC:Fish pond value {V_determ} at mean stopping time {np.mean(tau_determ)}')
        I = np.argmax(obj_stoch,axis=0)
        V_stoch_ant=np.max(obj_stoch,axis=0)
        tau_stoch_ant=t[I]
        print(f'Stoch Anticipative: Fish pond value {np.mean(V_stoch_ant)} at mean stopping time {np.mean(tau_stoch_ant)}')
        I = np.argmax(obj_determ,axis=0)
        V_determ_ant=np.max(obj_determ,axis=0)
        tau_determ_ant=t[I]
        print(f'Determ Anticipative: Fish pond value {np.mean(V_determ_ant)} at mean stopping time {np.mean(tau_determ_ant)}')

        if savedir == '':
            traindir=''
        else:
            traindir=savedir+'/Models/'+self.model
        if self.trainBoundary:
            print('Neural Comparison')
            DDCstoch=trainDecision(LSMC_stoch_cont_train,LSMC_stoch_exercise_train,4,Cont_val=LSMC_stoch_cont_val,Exercise_val=LSMC_stoch_exercise_val,batch_size = 128,dtype = self.dtype,savedir=traindir,model='stoch')
            DDCdeterm=trainDecision(LSMC_determ_cont_train,LSMC_determ_exercise_train,2,Cont_val=LSMC_determ_cont_val,Exercise_val=LSMC_determ_exercise_val,batch_size = 128,dtype = self.dtype,savedir=traindir,model='determ')

            tau_stoch_ddc_stoch,V_stoch_ddc_stoch=deepDecision(t,obj_stoch,DDCstoch,S)
            tau_stoch_ddc_determ,V_stoch_ddc_determ=deepDecision(t,obj_stoch,DDCdeterm,S[:,:,:2])
            tau_determ_ddc_determ,V_determ_ddc_determ=deepDecision(t,obj_determ,DDCdeterm,S[:,:,:2])

            print(f'\tStoch-Stoch Deep Decision: Fish pond value {np.mean(V_stoch_ddc_stoch)} at mean stopping time {np.mean(tau_stoch_ddc_stoch)}')
            print(f'\tStoch-Determ Deep Decision: Fish pond value {np.mean(V_stoch_ddc_determ)} at mean stopping time {np.mean(tau_stoch_ddc_determ)}')
            print(f'\tDeterm-Determ Deep Decision: Fish pond value {np.mean(V_determ_ddc_determ)} at mean stopping time {np.mean(tau_determ_ddc_determ)}')
            print(f'\tStoch Revenue = {np.mean(V_stoch_ddc_stoch)/np.mean(V_stoch_ddc_determ)}* Determ Revenue')

        print('Pathwise comparsion')
        V_stoch_stoch=0
        for wi in range(tau_stoch.size):
            tau = tau_stoch[wi]
            ti=np.argwhere(t==tau)[0,0]
            V_stoch_stoch+=obj_stoch[ti,wi]
        V_stoch_stoch/=obj_stoch.shape[1]
        print(f'\tStoch-Stoch pathwise: Fish pond value {V_stoch_stoch} at mean stopping time {np.mean(tau_stoch)}')

        V_stoch_determ=0
        for wi in range(tau_determ.size):
            tau = tau_determ[wi]
            ti=np.argwhere(t==tau)[0,0]
            V_stoch_determ+=obj_stoch[ti,wi]
        V_stoch_determ/=obj_stoch.shape[1]
        print(f'\tStoch-Determ pathwise: Fish pond value {V_stoch_determ} at mean stopping time {np.mean(tau_determ)}')
        print(f'\tStoch Revenue = {V_stoch_stoch/V_stoch_determ}* Determ Revenue')


if __name__=="__main__":
    import time
    'Salmon'
    # mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
    salmonParam=[0.818, 0.219, 0.163, 0.495, 1.286, 0.630, 0.921, 0.0303, 27.81]

    'Soy'
    # mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
    soyParam=[0.1,2.0,1.0,1.5,0.05,0.01,0.6,0.0303,1]

    farm1=fishFarm(salmonParam,soyParam)
    tic=time.time()
    farm1.generateFishFarm(0,1000000)
    toc=time.time()
    print(f'Elapsed time {toc-tic} s')