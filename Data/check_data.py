# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:24:24 2014

@author: istvan
"""
import csv
import numpy as np
from scipy import special
from pandas import *
import matplotlib 
from scipy.stats import norm
from matplotlib import dates
import matplotlib.pyplot as plt
import datetime
import re 


plt.rc('font', family='serif')
  

def CheckData(DataFiles):
    #loading data
    for f in DataFiles: 
        print('XXXXXXXXX '+f+' XXXXXXXXXXX')
        #Reading in the file
        #df=DataFrame.from_csv('RawData/'+f).reset_index()
        df=DataFrame.from_csv(f).reset_index()
    
        #Dropping unnecessary columns 
        #df=df[['#RIC' , 'Date[G]', 'Time[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Exch Time']]
        df=df[['#RIC' , 'Date[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Time[G]','Exch Time']]
        #df=df[['#RIC' , 'Date[G]', 'Time[G]','GMT Offset', 'Type','Price', 'Volume']]
        # Forward filling missing values fro Bid Price and Ask Price 

        print('Number of rows: '+str(len(df)))
       
        #Filtering on trades
        trades=df[df['Type']=='Trade']
        print('Number of trades: '+str(len(trades)))
        
        trades=trades[(np.isnan(trades['Price'])==0) & ( np.isnan(trades['Volume'])==0) ]
        print('Number of trades without price: '+str(len(trades)))

        
        #trades=trades[:5000]


        trades_timeg=trades.sort_values(['Date[G]','Time[G]'])
        trades_timee=trades.sort_values(['Date[G]','Exch Time'])       
        
        price_timeg=np.array(trades_timeg['Price'])
        price_timee=np.array(trades_timee['Price'])
        
        price_diff=price_timeg-price_timee
        sum_diff=np.sum(price_diff)
        print(sum_diff)
        
 
DataFiles=['NYSE_201004_data_11.csv'] 
CheckData(DataFiles)