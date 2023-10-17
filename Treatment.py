# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:53:43 2023

@author: Marta Via
"""

import pandas as pd
import numpy as np
import glob
import os as os
import datetime as dt
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy import stats
import pymannkendall as mk

#Si una funció no té self, no es pot cridar des de fora.

class Basics:
    
    def __init__(self,x):
        self.x = 0
        pass
    
    def Hello(self):
        print("Hello. These are treatment basics.")
    
    def R2(self, a, b):
        c = pd.DataFrame({"a": a, "b": b})
        cm = c.corr(method='pearson')
        r = cm.iloc[0, 1]
        return (r**2).round(2)
    def slope(self, b, a):
        c = pd.DataFrame({"a": a, "b": b})
        mask = ~np.isnan(a) & ~np.isnan(b)
        a1 = a[mask]
        b1 = b[mask]
        if (a1.empty) or (b1.empty):
            s = np.nan
        else:
            s, intercept, r_value, p_value, std_err = linregress(a1, b1)
        return s.round(2), intercept.round(2)
    # Select dates
    def SelectDates(self, df, datein, dateout): #datein and dateout must be stringsº
        date_in=pd.to_datetime(datein, dayfirst=True)
        date_out=pd.to_datetime('01/05/2023', dayfirst=True)
        mask_dates= (df['datetime']> date_in) & (df['datetime'] <= date_out)
        df_out=df[mask_dates]
        return df_out
    #
    #
    def Data_Matcher(self, df_small, df_big, dt_name_small, dt_name_big):
        li=[]
        for i in range(len(df_small)):
            for j in range(0,len(df_big)):
                if df_big[dt_name_big].iloc[j] == df_small[dt_name_small].iloc[i]:
                    li.append(df_big.iloc[j])
        df_out=pd.DataFrame(li)
        return df_out

                    
    def averages(self, range_in, range_out, vals_in): # not only for dates
        vals_out=[]
        for i in range(1,len(range_out)):
            acum=[]
            for j in range(0,len(range_in)):
                if range_in[j]>range_out[i-1] and range_in[j]<= range_out[i]:
                    acum.append(vals_in[j])
            vals_out.append(pd.Series(acum).sum()/float(len(acum)))
        return vals_out           
    
    def averaging(self, df_big, dt_small, dt_big): #Needs to be revised
        dfi=[]
        for i in range(1, len(dt_small)):
           x = df_big[(dt_big>dt_small.iloc[i-1]) & (dt_big<=dt_small.iloc[i])].mean()
           dfi.append(x)
        df_small=pd.DataFrame(dfi)
        return df_small            
                
    def MannKendall(self, a, alpha):
        b=mk.original_test(a, alpha)
        # print(b)
        return b
    
    def MannKendall_seasonal(self, a, alpha):
        b=mk.seasonal_test(a, alpha)
        c=mk.seasonal_sens_slope(a, alpha)
        print(b)
        print('Slope: ', c)
        return b,c     
            
            
            
            
            
            
            
            