from FishFarm import fishFarm
import os
import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
os.environ["TF_ENABLE_ONEDNN_OPTS"]="1"
tf.config.set_visible_devices([], 'GPU')

'Production Costs in Table 4. Norce report'

'Salmon'
# mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
# salmonParam=[0.12, 0.23, 0.75, 2.6, 0.02, 0.01, 0.9, 0.57, 95] # down,down
salmonParam=[0.12, 0.23, 0.75, 2.6, 0.02, 0.2, 0.9, 0.57, 95] # down,up
# salmonParam=[0.12, 0.23, 0.75, 2.6, 0.02, 0.6, 0.9, 0.57, 95] # up,up

"Fish feeding 25% of production cost, disease 30%, harvest 10%. Total production cost = 50% of price = labor, smolt, ..."
salmonPrice=salmonParam[-1] #NOK/KG
harvestingCosts=salmonPrice*0.5*0.1 # roughly 10%
feedingCosts=salmonPrice*0.5*0.5
initialSalmon=0.5*salmonPrice+feedingCosts+harvestingCosts #we add the costs to salmon price since they are respected in the model, other costs are fixed and thus removed
salmonParam[-1]=initialSalmon
print(f'Feeding costs {feedingCosts} and Harvesting costs {harvestingCosts}')

'Soy'
# mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
# soyParam=[0.15, 0.5, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500] # low vol
soyParam=[0.15, 1, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500] # medium vol
# soyParam=[0.15, 2, 0.4, 1.2, 0.06, 0.14, 0.44, 0.0, 1500] # high vol
soyParam[-1]=feedingCosts

farm1=fishFarm(salmonParam,soyParam,fc=soyParam[-1],hc=harvestingCosts,trainBoundary=True)
farm1.compareStoppingTimes(savedir='Paper')
