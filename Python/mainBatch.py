from FishFarm import fishFarm
import os
import tensorflow as tf
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt 
# from scipy.stats import norm 
import scipy.stats as st
import time
from joblib import Parallel, delayed

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
os.environ["TF_ENABLE_ONEDNN_OPTS"]="1"
tf.config.set_visible_devices([], 'GPU')

M=10000 # LSMC simulations (batch size)
m=10000 # trials 
threads=1
'Salmon'
salmonParams=[[0.12, 0.23, 0.75, 2.6, 0.02, 0.01, 0.9, 0.57, 95],[0.12, 0.23, 0.75, 2.6, 0.02, 0.2, 0.9, 0.57, 95],[0.12, 0.23, 0.75, 2.6, 0.02, 0.6, 0.9, 0.57, 95]]
soyParams=[[0.15, 0.5, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500],[0.15, 1, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500],[0.15, 2, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500]]

for salmonParam in salmonParams:
    for soyParam in soyParams:

        'Correlation Matrix'
        rho = None

        "Fish feeding 25% of production cost, disease 30%, harvest 10%. Total production cost = 50% of price = labor, smolt, ..."
        salmonPrice=salmonParam[-1] #NOK/KG
        harvestingCosts=salmonPrice*0.5*0.1 # roughly 10%
        feedingCosts=salmonPrice*0.5*0.25
        initialSalmon=0.5*salmonPrice+feedingCosts+harvestingCosts #we add the costs to salmon price since they are respected in the model, other costs are fixed and thus removed
        salmonParam2=salmonParam.copy()
        salmonParam2[-1]=initialSalmon
        print(f'Feeding costs {feedingCosts} and Harvesting costs {harvestingCosts}')
        soyParam[-1]=feedingCosts # to save the right dataset, since initial price is not relevant for soy model

        farm1=fishFarm(salmonParam2,soyParam,fc=soyParam[-1],hc=harvestingCosts,trainBoundary=False,rho=rho,verbose=1)

        V_ss=[]
        tau_s_m=[]
        V_sd=[]
        tau_d_m=[]
        ri=[]
        V_path=[]

        #takes ca. 20min
        tic = time.time()
        res=Parallel(n_jobs=threads)(delayed(farm1.compareStoppingTimes)(M=M,savedir='Python',seed=2+i) for i in range(m))
        V_ss,tau_s_m,V_sd,tau_d_m,ri,V_path=zip(*res)
        print(f'Elapsed time {time.time()-tic}')

        d={'Stoch-Stoch Value':V_ss,'Stoch-Stoch Mean Stopping Time':tau_s_m,'Stoch-Determ Value':V_sd,'Stoch-Determ Mean Stopping Time':tau_d_m,'Rel Improvement':ri}
        name='salmon_'+'_'.join(['{0:1.2f}'.format(x) for x in salmonParam])+'_soy_'+'_'.join(['{0:1.3f}'.format(x) for x in soyParam])+'_rho_'+'_'.join( ['{0:1.3f}'.format(rho[0,2]) if rho is not None else ''])
        saveDir='Python/Statistics/'+name

        df = pd.DataFrame(data=d)
        dfstat = df.describe()

        print('Save DataFrames')

        ci95Lower = []
        ci95Upper = []
        for i in range(df.shape[1]):
            data=df.iloc[:,i]
            ci=st.norm.interval(alpha=0.95, loc=np.mean(data), scale=st.sem(data))
            ci95Lower.append(ci[0])
            ci95Upper.append(ci[1])
        dfstat.loc['CI 95 Lower']=ci95Lower
        dfstat.loc['CI 95 Upper']=ci95Upper

        df.to_csv(saveDir+'.csv',sep=';',decimal='.')
        dfstat.to_csv(saveDir+'_stats'+'.csv',sep=';',decimal=',')

        #del res,df,dfstat,V_ss,tau_d_m,tau_s_m,V_ss,V_path,farm1
        'Correlation Matrix'

        rho=np.ones((4,4),dtype=np.float32)*(0.2)
        np.fill_diagonal(rho,1)
        rho[0,1]=salmonParam[6]
        rho[1,0]=salmonParam[6]
        rho[2,3]=soyParam[6]
        rho[3,2]=soyParam[6]



        "Fish feeding 25% of production cost, disease 30%, harvest 10%. Total production cost = 50% of price = labor, smolt, ..."
        salmonPrice=salmonParam[-1] #NOK/KG
        harvestingCosts=salmonPrice*0.5*0.1 # roughly 10%
        feedingCosts=salmonPrice*0.5*0.25
        initialSalmon=0.5*salmonPrice+feedingCosts+harvestingCosts #we add the costs to salmon price since they are respected in the model, other costs are fixed and thus removed
        salmonParam2=salmonParam.copy()
        salmonParam2[-1]=initialSalmon
        print(f'Feeding costs {feedingCosts} and Harvesting costs {harvestingCosts}')
        soyParam[-1]=feedingCosts # to save the right dataset, since initial price is not relevant for soy model

        farm1=fishFarm(salmonParam2,soyParam,fc=soyParam[-1],hc=harvestingCosts,trainBoundary=False,rho=rho,verbose=-1)

        V_ss=[]
        tau_s_m=[]
        V_sd=[]
        tau_d_m=[]
        ri=[]
        V_path=[]

        #takes ca. 20min
        tic = time.time()
        res=Parallel(n_jobs=threads)(delayed(farm1.compareStoppingTimes)(M=M,savedir='Python',seed=2+i) for i in range(m))
        V_ss,tau_s_m,V_sd,tau_d_m,ri,V_path=zip(*res)
        print(f'Elapsed time {time.time()-tic}')

        d={'Stoch-Stoch Value':V_ss,'Stoch-Stoch Mean Stopping Time':tau_s_m,'Stoch-Determ Value':V_sd,'Stoch-Determ Mean Stopping Time':tau_d_m,'Rel Improvement':ri}
        name='salmon_'+'_'.join(['{0:1.2f}'.format(x) for x in salmonParam])+'_soy_'+'_'.join(['{0:1.3f}'.format(x) for x in soyParam])+'_rho_'+'_'.join( ['{0:1.3f}'.format(rho[0,2]) if rho is not None else ''])
        saveDir='Python/Statistics/'+name

        df = pd.DataFrame(data=d)
        dfstat = df.describe()

        print('Save DataFrames')

        ci95Lower = []
        ci95Upper = []
        for i in range(df.shape[1]):
            data=df.iloc[:,i]
            ci=st.norm.interval(alpha=0.95, loc=np.mean(data), scale=st.sem(data))
            ci95Lower.append(ci[0])
            ci95Upper.append(ci[1])
        dfstat.loc['CI 95 Lower']=ci95Lower
        dfstat.loc['CI 95 Upper']=ci95Upper

        df.to_csv(saveDir+'.csv',sep=';',decimal='.')
        dfstat.to_csv(saveDir+'_stats'+'.csv',sep=';',decimal=',')

        #del res,df,dfstat,V_ss,tau_d_m,tau_s_m,V_ss,V_path,farm1
        'Correlation Matrix'

        rho=np.ones((4,4),dtype=np.float32)*(-0.2)
        np.fill_diagonal(rho,1)
        rho[0,1]=salmonParam[6]
        rho[1,0]=salmonParam[6]
        rho[2,3]=soyParam[6]
        rho[3,2]=soyParam[6]


        "Fish feeding 25% of production cost, disease 30%, harvest 10%. Total production cost = 50% of price = labor, smolt, ..."
        salmonPrice=salmonParam[-1] #NOK/KG
        harvestingCosts=salmonPrice*0.5*0.1 # roughly 10%
        feedingCosts=salmonPrice*0.5*0.25
        initialSalmon=0.5*salmonPrice+feedingCosts+harvestingCosts #we add the costs to salmon price since they are respected in the model, other costs are fixed and thus removed
        salmonParam2=salmonParam.copy()
        salmonParam2[-1]=initialSalmon
        print(f'Feeding costs {feedingCosts} and Harvesting costs {harvestingCosts}')
        soyParam[-1]=feedingCosts # to save the right dataset, since initial price is not relevant for soy model

        farm1=fishFarm(salmonParam2,soyParam,fc=soyParam[-1],hc=harvestingCosts,trainBoundary=False,rho=rho,verbose=-1)

        V_ss=[]
        tau_s_m=[]
        V_sd=[]
        tau_d_m=[]
        ri=[]
        V_path=[]

        #takes ca. 20min
        tic = time.time()
        res=Parallel(n_jobs=threads)(delayed(farm1.compareStoppingTimes)(M=M,savedir='Python',seed=2+i) for i in range(m))
        V_ss,tau_s_m,V_sd,tau_d_m,ri,V_path=zip(*res)
        print(f'Elapsed time {time.time()-tic}')

        d={'Stoch-Stoch Value':V_ss,'Stoch-Stoch Mean Stopping Time':tau_s_m,'Stoch-Determ Value':V_sd,'Stoch-Determ Mean Stopping Time':tau_d_m,'Rel Improvement':ri}
        name='salmon_'+'_'.join(['{0:1.2f}'.format(x) for x in salmonParam])+'_soy_'+'_'.join(['{0:1.3f}'.format(x) for x in soyParam])+'_rho_'+'_'.join( ['{0:1.3f}'.format(rho[0,2]) if rho is not None else ''])
        saveDir='Python/Statistics/'+name

        df = pd.DataFrame(data=d)
        dfstat = df.describe()

        print('Save DataFrames')

        ci95Lower = []
        ci95Upper = []
        for i in range(df.shape[1]):
            data=df.iloc[:,i]
            ci=st.norm.interval(alpha=0.95, loc=np.mean(data), scale=st.sem(data))
            ci95Lower.append(ci[0])
            ci95Upper.append(ci[1])
        dfstat.loc['CI 95 Lower']=ci95Lower
        dfstat.loc['CI 95 Upper']=ci95Upper

        df.to_csv(saveDir+'.csv',sep=';',decimal='.')
        dfstat.to_csv(saveDir+'_stats'+'.csv',sep=';',decimal=',')


