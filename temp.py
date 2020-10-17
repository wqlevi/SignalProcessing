# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import scipy.io
import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from statsmodels.tsa.arima_model import ARIMA
from statsmodels.tsa.stattools import adfuller
from pandas import DataFrame
from sklearn.metrics import mean_squared_error

TR = 0.1 # sampling time fMRI
fmri_dummy = np.ones((1,6400))

os.chdir('/big_data/qi/05232019/05232019_analyzed/15/Results_Slice3/')
fmri_3_tmp = scipy.io.loadmat('line_scanning_data_s3.mat')
fmri_3 = np.absolute(fmri_3_tmp['total_cortical_depth_map'])

t_fmri = np.linspace(0,np.size(fmri_dummy,1)*TR,np.size(fmri_dummy,1))

# 2D heat-map fMRI
sns.heatmap(fmri_3,cmap = 'Greys_r')

# 1D timecourse fMRI
plt.figure(figsize = (8,2))
plt.plot(t_fmri,np.mean(fmri_3,0))
plt.xlim(0,640)

plt.ylabel('Raw')
plt.show()

ts = np.mean(fmri_3,0)
# auto-correlation
plt.figure(figsize = (8,2))
plt.acorr(ts,maxlags = 20,usevlines = False)
# Applying ARIMA
model = ARIMA(ts,order=(0,1,0))
model_fit = model.fit(disp=0)
print(model_fit.summary())
residuals = DataFrame(model_fit.resid)
residuals.plot(kind = 'kde')
plt.show()
# ARIMA prediction
train, test = [ts[:int(len(ts)/4)],ts[int(len(ts)/4):]]
history = [x for x in train]
prediction = list()
for t in range (len(test)):
    model = ARIMA(history, order  = (0,1,0))
    model_fit = model.fit(disp = 0)
    output = model_fit.forecast()
    yhat = output[0]
    prediction.append(yhat)
    obs = test[t]
    history.append(obs)
    print('predicted = %f, expected = %f' % (yhat,obs))
error = mean_squared_error(test,prediction)
print('Test MSE: %.df' % error)

plt.plot(test+1,label = 'raw')
plt.plot(prediction,color = 'red',label = 'prediction')
plt.legend(loc = 'best')
plt.show()
# Rolling mean and std
ts_1 = DataFrame(ts)
rl_mean = ts_1.rolling(window = 100).mean()
rl_std = ts_1.rolling(window = 100).std()
plt.figure(figsize = (16,8))
plt.plot(ts_1,color = 'blue',label = 'raw')
plt.plot(rl_mean,color = 'red',label = 'Rolling Mean')
plt.plot(rl_std,color = 'green',label = 'Rolling std')
plt.legend(loc = 'best')
plt.show()
# adfuller test
results = adfuller(ts)
print('ADF Stats: {}'.format(results[0]))
print('p-Value: {}'.format(results[1]))
print('Critical Values:')
for key, value in results[4].items():
    print ('\t{}: {}'.format(key, value))
if results[1] < 0.05:
    print ('Null rejected,series are stationary')
else:
    print ('Null accepted,series are not stationary')
