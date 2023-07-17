import numpy as np
import tensorflow as tf
from tensorflow.keras import Model
from tensorflow.keras.layers import BatchNormalization, Dense
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.optimizers import Adam

from joblib import Parallel,delayed
from tqdm import tqdm
from pathlib import Path
import time
import matplotlib.pyplot as plt
from typing import Dict, Optional
from numpy.typing import DTypeLike


class DeepDecisionClassifier(Model):
    def __init__(self,d):
        super().__init__()
        # use batch_size = None for variable batch size
        # d stochastic variables at time i
        # Input(shape=[d], batch_size=None, name='input')
        self.d=d
        self.B1 = BatchNormalization()
        # self.B2 = BatchNormalization()
        # self.B3 = BatchNormalization()
        # self.B4 = BatchNormalization()
        self.B5 = BatchNormalization()
        self.D1 = Dense(d*16, activation = 'relu')
        self.D2 = Dense(d*32, activation = 'relu')
        # self.D3 = Dense(d*128, activation = 'tanh')
        self.D4 = Dense(d*16, activation = 'sigmoid')
        self.D5 = Dense(2, activation = 'softmax')

    def call(self,inputs,training=False):
        x = self.B1(inputs)
        x = self.D1(x)
        # x = self.B2(x)
        x = self.D2(x)
        # x = self.B3(x)
        # x = self.D3(x)
        # x = self.B4(x)
        x = self.D4(x)
        x = self.B5(x)
        x = self.D5(x)
        return x
    
    def train_step(self,data):
        input = data[:,:self.d]
        label = data[:,self.d:]
        with tf.GradientTape() as tape:
            p = self(input,training=True)
            loss = self.compiled_loss(p,label,regularization_losses=self.losses)

        var_list = self.trainable_variables
        gradients = tape.gradient(loss, var_list)
        self.optimizer.apply_gradients(zip(gradients, var_list))

        self.compiled_metrics.update_state(p,label)

        return {m.name: m.result() for m in self.metrics}

def loadDecisionModel(path:str,model:str):
    files = Path(path).glob(model+'_*')
    DDCs = {}

    for f in files:
        tmp = f.name
        key = int(f.name.replace(model+'_',''))
        value = tf.keras.models.load_model(str(f.parent)+'/'+tmp)
        DDCs[key]=value

    return DDCs

def deepDecision(t:np.ndarray,obj:np.ndarray,DDCs:Dict[int,tf.keras.Model],S:np.ndarray):
    [N,M,d]=S.shape
    paths = np.arange(M,dtype=np.int64)
    totalPaths = np.arange(M,dtype=np.int64)
    V=obj[-1].copy()
    tau = t[-1]*np.ones_like(V)
    ex = np.zeros((M,1),dtype=np.bool8)
    for i in range(1,N-1):
        if i in DDCs:
            DDC = DDCs[i]
            p = DDC(S[i])
            ex = (1 - np.argmax(p,axis=1)).astype(np.bool8)
        else:
            continue

        wiEx = np.intersect1d(paths, totalPaths[ex])
        paths = np.setdiff1d(paths,totalPaths[ex])
        V[wiEx]=obj[i,wiEx]
        tau[wiEx]=t[i]
    return tau,V


# def trainSingleDecision(Cont_train:np.ndarray,
#                   Exercise_train:np.ndarray,
#                   d:int,
#                   i:int,
#                   batch_size:int = 128,
#                   dtype:DTypeLike = np.float32,
#                   path:Optional[str]='',
#                   savedir:Optional[str]='',
#                   model:Optional[str]=''):
    
#     exData=Exercise_train[i]
#     if len(exData)==0 or exData.shape[0]<batch_size:
#         return None
#     contData=Cont_train[i]
#     epochs=int(contData.shape[0]/exData.shape[0])
#     steps_per_epoch=int(exData.shape[0]/batch_size)

#     def gen():
#         j=0
#         stride=exData.shape[0]
#         while True:
#             if (j+1)*stride >= contData.shape[0]:
#                 j=0
#             cont=contData[j*stride:(j+1)*stride,:]
#             ex = exData[:,:]
#             for jj in range(int(stride/batch_size)):
#                 c = cont[jj*batch_size:(jj+1)*batch_size,:]
#                 e = ex[jj*batch_size:(jj+1)*batch_size,:]
#                 dataInput = np.concatenate([e,c],axis=0)
#                 dataClass = np.concatenate([np.stack([np.ones((e.shape[0]),dtype=dtype),np.zeros((e.shape[0]),dtype=dtype)],axis=1),
#                                     np.stack([np.zeros((c.shape[0]),dtype=dtype),np.ones((c.shape[0]),dtype=dtype)],axis=1)],axis=0)
#                 yield np.concatenate([dataInput,dataClass],axis=1)
#             j=j+1

#     dataset = tf.data.Dataset.from_generator(gen,output_signature=(tf.TensorSpec(shape=(2*batch_size,d+2))))

#     DDC = DeepDecisionClassifier(d)
#     DDC.compile(optimizer='adam',loss=CategoricalCrossentropy(from_logits=False),metrics='categorical_accuracy', jit_compile=False, run_eagerly=False)

#     tic = time.time()
#     DDC.fit(dataset,epochs=epochs,steps_per_epoch=steps_per_epoch,verbose=0)
#     ctime = time.time() - tic
    
#     if savedir != '':
#         DDC.save(str(path)+'/'+model+'_'+f'{i:02d}')

#     return DDC
    
# def trainDecision(Cont_train:np.ndarray,
#                   Exercise_train:np.ndarray,
#                   d:int,
#                   Cont_val:Optional[np.ndarray]=None,
#                   Exercise_val:Optional[np.ndarray]=None,
#                   batch_size:int = 128,
#                   dtype:DTypeLike = np.float32,
#                   seed:int = 1234,
#                   savedir:Optional[str]='',
#                   model:Optional[str]=''):
#     path = Path().cwd()/Path(savedir)
#     models = path.glob(model+'_*')
#     decision={}
#     if not Path.exists(path) or (savedir == '') or len(list(models))==0:
#         print('Train Classifier')
#         N=Cont_train.shape[0]
#         DDClist = Parallel(n_jobs=12)(delayed(trainSingleDecision(Cont_train,Exercise_train,d,i,batch_size=batch_size,dtype=dtype,savedir=savedir,model=model,path=path)) for i in range(1,N-1)) 
#         for j in range(0,len(DDClist)):
#             DDC = DDClist[j]
#             if len(DDC)>0:
#                 decision[j+1]=DDC
#     else:
#         print('Load Classifier')
#         decision=loadDecisionModel(path,model)
#     return decision
    

def trainDecision(Cont_train:np.ndarray,
                  Exercise_train:np.ndarray,
                  d:int,
                  Cont_val:Optional[np.ndarray]=None,
                  Exercise_val:Optional[np.ndarray]=None,
                  batch_size:int = 128,
                  dtype:DTypeLike = np.float32,
                  seed:int = 1234,
                  savedir:Optional[str]='',
                  model:Optional[str]=''):
    
    path = Path().cwd()/Path(savedir)
    models = path.glob(model+'_*')
    decision={}
    if not Path.exists(path) or (savedir == '') or len(list(models))==0:
        print('Train Classifier')
        if savedir != '':
            path.mkdir(exist_ok=True,parents=True)
        np.random.seed(seed)
        tf.random.set_seed(seed)


        N=Cont_train.shape[0]

        for i in (pbar:=tqdm(range(1,N-1))):

            exData=Exercise_train[i]
            if len(exData)==0 or exData.shape[0]<batch_size:
                continue
            contData=Cont_train[i]
            epochs=int(contData.shape[0]/exData.shape[0])
            steps_per_epoch=int(exData.shape[0]/batch_size)

            def gen():
                j=0
                stride=exData.shape[0]
                while True:
                    if (j+1)*stride >= contData.shape[0]:
                        j=0
                    cont=contData[j*stride:(j+1)*stride,:]
                    ex = exData[:,:]
                    for jj in range(int(stride/batch_size)):
                        c = cont[jj*batch_size:(jj+1)*batch_size,:]
                        e = ex[jj*batch_size:(jj+1)*batch_size,:]
                        dataInput = np.concatenate([e,c],axis=0)
                        dataClass = np.concatenate([np.stack([np.ones((e.shape[0]),dtype=dtype),np.zeros((e.shape[0]),dtype=dtype)],axis=1),
                                            np.stack([np.zeros((c.shape[0]),dtype=dtype),np.ones((c.shape[0]),dtype=dtype)],axis=1)],axis=0)
                        yield np.concatenate([dataInput,dataClass],axis=1)
                    j=j+1

            dataset = tf.data.Dataset.from_generator(gen,output_signature=(tf.TensorSpec(shape=(2*batch_size,d+2))))

            DDC = DeepDecisionClassifier(d)
            DDC.compile(optimizer='adam',loss=CategoricalCrossentropy(from_logits=False),metrics='categorical_accuracy', jit_compile=False, run_eagerly=False)

            tic = time.time()
            DDC.fit(dataset,epochs=epochs,steps_per_epoch=steps_per_epoch,verbose=0)
            ctime = time.time() - tic

            if Exercise_val is not None:
                exDataVal=Exercise_val[i]
                contDataVal=Cont_val[i]

                ex_pred=DDC(exDataVal).numpy()
                cont_pred=DDC(contDataVal).numpy()
                cont_err=np.mean(np.argmax(cont_pred,axis=1))
                ex_err=np.mean(np.argmax(ex_pred,axis=1))
                ex_inc=np.sum(np.argmax(cont_pred,axis=1)==0)/ex_pred.size
                pbar.set_postfix({'Training':ctime,'Cont Err':cont_err,'Ex Err':ex_err,'Ex increase':ex_inc})
            else:
                pbar.set_postfix({'Training':ctime})

            decision[i]=DDC
            if savedir != '':
                DDC.save(str(path)+'/'+model+'_'+f'{i:02d}')
    else:
        print('Load Classifier')
        decision=loadDecisionModel(path,model)
    return decision