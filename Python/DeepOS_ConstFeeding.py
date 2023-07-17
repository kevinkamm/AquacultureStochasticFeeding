"""
This code is adapted from from arXiv:1908.01602v3

"""

import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras import Model, Sequential, Input
from tensorflow.keras.layers import BatchNormalization, Dense, Reshape, Activation
# from tensorflow.keras.losses import BinaryCrossentropy, MeanSquaredError
from tensorflow.keras.optimizers import Adam

import time as timer
from tqdm import tqdm, trange

from typing import List,Tuple
from numpy.typing import DTypeLike

class DeepOSNet(Model):
    def __init__(self,
                 hiddenDims: List[int],
                 outputDims: int,
                 d : int,
                 N : int) \
            -> None:
        self.hiddenDims = hiddenDims
        self.outputDims = outputDims
        self.d = d
        self.N = N

    def build(self) \
            -> Model:
        model = Sequential(name='DeepOS')
        model.add(Reshape((self.N*self.d,),name='Reshape_0_1',input_shape=(self.N , self.d)))
        model.add(BatchNormalization(epsilon=1e-6,axis=-1,momentum=0.9))
        model.add(Reshape((self.N , self.d), name='Reshape_0_2'))
        for i, hD in enumerate(self.hiddenDims):
            model.add(Dense(units=hD,
                            activation=None,
                            name=f'deep_os_hidden_{i + 1}'))
            model.add(Reshape((self.N * hD,), name=f'Reshape_{i + 1}_1', input_shape=(self.N, hD)))
            model.add(BatchNormalization(epsilon=1e-6, axis=-1, momentum=0.9))
            model.add(Reshape((self.N, hD), name=f'Reshape_{i + 1}_2'))
            model.add(Activation('relu'))
        i += 1
        model.add(Dense(units=self.outputDims,
                        activation=None,
                        name='deep_os_out'))
        model.add(Reshape((self.N * self.outputDims,), name=f'Reshape_{i + 1}_1', input_shape=(self.N, self.outputDims)))
        model.add(BatchNormalization(epsilon=1e-6, axis=-1, momentum=0.9))
        model.add(Reshape((self.N, self.outputDims), name=f'Reshape_{i + 1}_2'))
        model.add(Activation('sigmoid'))

        return model


class DeepOS:
    def __init__(self,
                 batch_size : int,
                 hiddenDims : List[int],
                 d : int,
                 N : int,
                 T : float,
                 m,cr,n0,hc,fc,wInf,a,b,c,r,kappa,alpha,sigma1,sigma2,rho,l,delta0,P0):
        """

        :param batch_size:
        :param hiddenDims:
        :param d:
        :param N:
        :param T:
        """

        "Parameters for fish farming"
        self.m = m
        self.cr = cr
        self.n0 = n0
        self.hc = hc
        self.fc = fc
        self.wInf = wInf

        "Parameters for Bertalanffyos growth function"
        self.a = a
        self.b = b
        self.c = c

        "Model parameters"
        self.r = r
        self.kappa = kappa
        self.alpha = alpha
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.rho = rho
        self.l = l
        self.delta0 = delta0
        self.P0 = P0

        "Parameters"
        # batch size
        self.batch_size = batch_size
        # dimension of stochastic process
        self.d = 2  # legacy parameter
        # time steps
        self.N = N
        self.Nsim = 10*N  # for the simulation of Schwartz-2-factor
        # # time
        self.T = T

        "DNN Configuration"
        # Placeholder for real data
        X = Input(shape=[N-1, d+1], batch_size=batch_size, name='StochProcess_Objective')
        self.net = DeepOSNet(hiddenDims,1,d+1,N-1).build()
        X_tilde = self.net(X)
        self.net_model=Model(inputs=X,
                             outputs=X_tilde,
                             name='DeepOSModel')
        # print(self.net_model.summary())
        # print(self.net.summary())



    @tf.function
    def train_deep_os(self,p,x,opt):
        # self.initialize()
        # print(opt.learning_rate(opt.iterations))
        with tf.GradientTape() as tape:
            nets = self.net_model(tf.transpose(tf.concat([x[:, :, :-1], p[:, :, :-1]], axis=1),(0,2,1)),training=True)
            # nets = self.net_model(tf.transpose(tf.concat([x[:, :, :-1], p[:, :, :-1]], axis=1), (0, 2, 1)),
            #                       training=False)
            nets = tf.transpose(nets,(0,2,1))
            u_list = [nets[:, :, 0]]
            u_sum = u_list[-1]
            for k in range(1, self.N - 1):
                u_list.append(nets[:, :, k] * (1. - u_sum))
                u_sum += u_list[-1]

            u_list.append(1. - u_sum)
            u_stack = tf.concat(u_list, axis=1)
            p = tf.squeeze(p, axis=1)
            loss = tf.reduce_mean(tf.reduce_sum(-u_stack * p, axis=1))

        var_list = self.net.trainable_variables
        gradients = tape.gradient(loss,var_list)
        opt.apply_gradients(zip(gradients,var_list))

        idx = tf.argmax(tf.cast(tf.cumsum(u_stack, axis=1) + u_stack >= 1,
                                dtype=tf.uint8),
                        axis=1,
                        output_type=tf.int32)
        stopped_payoffs = tf.reduce_mean(tf.gather_nd(p, tf.stack([tf.range(0, self.batch_size, dtype=tf.int32), idx],
                                                                  axis=1)))
        return loss, stopped_payoffs

    # @tf.function
    def objective(self)\
            -> Tuple[tf.Tensor,tf.Tensor]:
        """

        :return:
        """
        t = tf.linspace(0.,T,self.Nsim)
        dt=self.T / (self.Nsim-1)
        "correlated Brownian motions"
        W = tf.concat([tf.zeros((self.batch_size,1,2)),
                       tf.cumsum(tf.random.normal(
                              shape=(self.batch_size, self.Nsim-1, 2),
                              stddev=tf.math.sqrt(dt)),
                            axis=1)],
                       1)
        W1 = W[:,:,0]
        W2 = self.rho * W1 + tf.math.sqrt(1-self.rho**2) * W[:,:,1]
        dW2 = W2[:,1:]-W2[:,0:-1]

        "convenience yield"
        expkt = tf.reshape(tf.exp(-self.kappa*t),(1,-1))
        expkt2 = tf.reshape(tf.exp(self.kappa*t),(1,-1))
        I_expkt_dW2 = tf.concat([tf.zeros((self.batch_size,1)),
                                 tf.cumsum(expkt2[:,:-1]*dW2,1)],
                                1)
        delta = self.delta0*expkt+ \
                (self.alpha-self.l/self.kappa)*(1-expkt)+\
                expkt*self.sigma2*I_expkt_dW2

        "spot price"
        I_delta_dt = tf.concat([tf.zeros((self.batch_size,1)),
                                tf.cumsum(delta[:,:-1]*dt,1)],
                               1)


        # <--todo improve, slice definition of P
        P = self.P0 * tf.exp((self.r-self.sigma1**2/2)*t-I_delta_dt+self.sigma1*W1)

        delta_t = tf.reshape(delta[:,int(self.Nsim/self.N)-1::int(self.Nsim/self.N)],(self.batch_size,1,self.N))
        P_t = tf.reshape(P[:, int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)],
                             (self.batch_size, 1, self.N))
        # todo improve, slice definition of P-->

        "Bertalanffy’s growth function"
        wt = wInf * (self.a - self.b * tf.exp(-self.c*t))**3
        dwt = tf.concat([tf.constant([0.]),wt[1:]-wt[:-1]],0)/dt

        "Number of fish"
        nt = self.n0 * tf.exp(-self.m*t)

        "Total biomass (kg)"
        Xt = nt * wt

        "Harvesting cost"
        CH = Xt * self.hc

        "Feeding cost"
        CF = nt * dwt * self.cr * self.fc
        I_CF_dt = tf.cumsum(tf.concat([tf.constant([0.]), tf.exp(-self.r * t[1:]) * CF[1:] * dt], 0), 0)

        "Only use coarse grid for optimization"
        wt = tf.reshape(wt[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))
        nt = tf.reshape(nt[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))
        Xt = tf.reshape(Xt[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))
        CH = tf.reshape(CH[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))
        CF = tf.reshape(CF[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))
        I_CF_dt = tf.reshape(I_CF_dt[int(self.Nsim / self.N) - 1::int(self.Nsim / self.N)], (1, 1, -1))


        "Joint process"
        x = tf.concat([delta_t,P_t],1)

        "Value of objective"
        tcoarse=tf.reshape(t[int(self.Nsim/self.N)-1::int(self.Nsim/self.N)],(1,1,-1))
        p = tf.exp(-self.r*tcoarse)*(P_t*Xt-CH)-I_CF_dt

        return x,p



    def train_model(self,
                    epochs : int,
                    lr_boundaries,
                    lr_values,
                    beta1=0.9,
                    beta2=0.999,
                    epsilon=1e-8):
        learning_rate_fn = tf.keras.optimizers.schedules.PiecewiseConstantDecay(lr_boundaries, lr_values)
        opt=Adam(learning_rate_fn,
                   beta_1=beta1,
                   beta_2=beta2,
                   epsilon=epsilon)
        # opt = Adam(learning_rate=1e-4,
        #            beta_1=beta1,
        #            beta_2=beta2,
        #            epsilon=epsilon)
        for _ in (pbar := tqdm(range(epochs), desc='Train DeepOS')):
            x,p = self.objective()
            loss, stopped_payoffs = self.train_deep_os(p,x,opt)
            pbar.set_postfix({'loss': loss, 'stopped payoffs': stopped_payoffs})

    def simulate_price(self,
                       mc_runs : int,
                       timeBoundary : int = -1):
        px_mean=0
        t = tf.linspace(0.,T,self.Nsim)
        salmonP={i:[] for i in range(0,self.N-1)} 
        salmonDelta={i:[] for i in range(0,self.N-1)} 
        tau_mean=0
        for _ in (pbar := tqdm(range(mc_runs), desc='Simulate Opt Stopping')):
            x,p = self.objective()
            nets = self.net(tf.transpose(tf.concat([x[:,:, :-1], p[:, :,:-1]], axis=1), (0,2, 1)),training=False)
            # nets = self.net(tf.transpose(tf.concat([x[:, :, :-1], p[:, :, :-1]], axis=1), (0, 2, 1)), training=False)
            nets = tf.transpose(nets, (0,2, 1))
            u_list = [nets[:, :, 0]]
            u_sum = u_list[-1]
            for k in range(1, self.N - 1):
                u_list.append(nets[:, :, k] * (1. - u_sum))
                u_sum += u_list[-1]

            u_list.append(1. - u_sum)
            u_stack = tf.concat(u_list, axis=1)
            p = tf.squeeze(p, axis=1)

            idx = tf.argmax(tf.cast(tf.cumsum(u_stack, axis=1) + u_stack >= 1,
                                    dtype=tf.uint8),
                            axis=1,
                            output_type=tf.int32)
            
            if timeBoundary>-1:
                for i in range(0,self.N-1):
                    exercise = tf.cast(idx==i,tf.int32)
                    tmpP=tf.boolean_mask(x[:,1,i],exercise)
                    salmonP[i].append(tmpP.numpy())
                    tmpDelta=tf.boolean_mask(x[:,0,i],exercise)
                    salmonDelta[i].append(tmpDelta.numpy())

            stopped_payoffs = tf.reduce_mean(tf.gather_nd(p, tf.stack([tf.range(0, self.batch_size, dtype=tf.int32), idx],
                                                                      axis=1)))

            tau = tf.reshape(t[int(self.Nsim/N)-1::int(self.Nsim/N)],(1,-1))*tf.ones((self.batch_size,1))
            stopped_time = tf.reduce_mean(
                tf.gather_nd(tau, tf.stack([tf.range(0, self.batch_size, dtype=tf.int32), idx],
                                         axis=1)))
            px_mean += stopped_payoffs
            tau_mean += stopped_time
            pbar.set_postfix({'stopped payoffs': stopped_payoffs, 'mean stopping time': stopped_time})

        if timeBoundary>-1:
            for i in range(0,self.N-1):
                salmonP[i]=np.concatenate(salmonP[i])
                salmonDelta[i]=np.concatenate(salmonDelta[i])

            return px_mean / mc_runs, tau_mean.numpy() / mc_runs, [salmonP,salmonDelta]
        else:
            return px_mean / mc_runs, tau_mean.numpy() / mc_runs


if __name__ == "__main__":
    # tf.config.set_visible_devices([], 'GPU')
    # print(tf.config.experimental.get_synchronous_execution())
    # print(tf.config.experimental.list_physical_devices())
    # print(tf.config.threading.get_inter_op_parallelism_threads())
    # print(tf.config.threading.get_intra_op_parallelism_threads())
    # from tensorflow.python.framework.ops import disable_eager_execution
    # disable_eager_execution()
    import matplotlib.pyplot as plt
    import os
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
    os.environ["TF_ENABLE_ONEDNN_OPTS"]="1"
    tf.random.set_seed(1234)

    # Simulation parameters
    N = 3*24  # internal simulation is 10*N time steps
    batch_size = int(8192/4)  # 2^i for max performance on CUDA device
    d = 2  # legacy parameter

    # Parameters for fish farming, see Table 4 Ewald2017
    T = 3
    m = .1
    cr = 1.1
    n0 = 10000
    hc = 3
    # fc = 7
    fc = 10
    wInf = 6

    # Bertalanffyo's growth function, see Footnote 20 Ewald2017
    a = 1.113
    b = 1.097
    c = 1.43

    # Model parameters for Panel B, see Table 2 Ewald2017
    r = 0.0303  # see page 5 Norwegian average interest rate

    # salmonParam=[0.121539, 0.319185, 1.29164, 3.98723, 0.0351972, 0.000327077, 0.915016, 0, 27.81]
    # salmonParam=[0.174557, 0.303309, 0.185905, 5.50823e-06, 0.0349834, 0.0648319, 0.897994, 0, 27.81]
    salmonParam=[0.818, 0.219, 0.163, 0.495, 1.286, 0.630, 0.921, 0.0303, 27.81]
    mu,sigma1,sigma2,kappa,alpha,l,rho,delta0,P0 = salmonParam

    lr_values = [0.05, 0.005, 0.0005]
    mc_runs = 500

    ticDeepOS = timer.time()
    neurons = [d + 50, d + 50]
    train_steps = 3000 + d
    lr_boundaries = [int(500 + d / 5), int(1500 + 3 * d / 5)]

    deepOS = DeepOS(batch_size,neurons,d,N,T,m,cr,n0,hc,fc,wInf,a,b,c,r,kappa,alpha,sigma1,sigma2,rho,l,delta0,P0)
    deepOS.train_model(train_steps,lr_boundaries,lr_values)
    px_mean, tau_mean, exerciseRegion = deepOS.simulate_price(mc_runs,timeBoundary=1)

    Nsim=N*10
    t = np.linspace(0,T,Nsim,endpoint=True)
    tCoarse = t[int(Nsim/N)-1::int(Nsim/N)]
    tn = (tCoarse>=tau_mean).nonzero()[0][0]

    for i in range(-3,3):
        if tn+i<N-1 and tn+i>=0:
            plt.figure()
            plt.scatter(exerciseRegion[0][tn+i],exerciseRegion[1][tn+i])
            plt.title(f'at time t={tCoarse[tn+i]}, i={tn+i}')
            plt.savefig(f'Figures/const {i}.png')
            plt.close()
