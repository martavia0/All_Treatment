# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 10:11:34 2023

@author: Marta Via
"""

#%%import pandas as pd
import numpy as np
import glob
import os as os
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy import stats
path_py="C:/Users/maria/Documents/Marta Via/1. PhD/F. Scripts/Python Scripts"
os.chdir(path_py)
from Treatment import *
#%%
trt = Basics(5)
trt.Hello()
print(trt.x)
#%%
pm1 = ['Org', 'SO4', 'NO3', 'NH4', 'Chl', 'BClf', 'BCsf'] #labels_pm1
nr =  ['Org', 'SO4', 'NO3', 'NH4', 'Chl'] #labels_nrpm1
factors = ['COA', 'HOA', 'BBOA', 'Amine-OA', 'LO-OOA', 'MO-OOA'] #labels_factors
c_df=['green', 'red', 'blue', 'orange', 'fuchsia', 'grey', 'saddlebrown'] #colors_nrpm1
c_oa=['mediumpurple', 'grey', 'saddlebrown', 'steelblue', 'lightgreen', 'darkgreen'] #colors_factors
#%%
path_treatment="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/All Series/Treatment"
os.chdir(path_treatment)
df=pd.read_csv('Chemistry_PM1.txt', sep='\t')
#%% Fixing time
dt=pd.to_datetime(df['time_utc'], dayfirst=True)
df.index=dt
df['datetime']=df.index
#%% Overall time series
fig, axs=plt.subplots(figsize=(20,4))
df.plot(y=nr, ax=axs, color=c_df)
axs.grid('x')
axs.legend(loc='upper left')
#%% Replacing the negatives below their DL (0.148, 0.0224, 0.012, 0.284, 0.011)
df.Org[df.Org < -0.0] = np.nan
df.SO4[df.SO4 < -0.0] = np.nan
df.NO3[df.NO3 < -0.0] = np.nan
df.NH4[df.NH4 < -0.0] = np.nan
df.Chl[df.Chl < -0.0] = np.nan
# Replotting
fig, axs=plt.subplots(figsize=(15,4))
df.plot(y=['Org', 'SO4', 'NO3', 'NH4', 'Chl'], ax=axs, color=['green', 'red', 'blue', 'orange', 'fuchsia','grey','saddlebrown'],subplots=True)
axs.grid('x')
axs.legend(loc='upper left')
#%% Diel yearly plots PM1 compounds
df['Hour']=df['datetime'].dt.hour
df['Year']=df['datetime'].dt.year
year_diel=df[['Org', 'SO4', 'NO3', 'NH4', 'Chl']].groupby([df['Year'], df['Hour']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
for i in range(0,10):
    axs.plot([i*24,i*24],[0,8], c='grey')
year_diel.plot(color=['green', 'red', 'blue', 'orange', 'fuchsia'], ax=axs)
axs.grid('x')
axs.set_xticks(range(0,216,12))
# axs.set_xticklabels([0,'2017\n00:00','2018\n00:00','2019\n00:00','2020\n00:00','2021\n00:00','2022\n00:00','2023\n00:00','',''])
axs.set_xticklabels(['2014\n00:00','','2015\n00:00','','2017\n00:00','', '2018\n00:00','','2019\n00:00','','2020\n00:00','','2021\n00:00','',
                     '2022\n00:00','','2023\n00:00','',])
axs.set_xlabel('Year - Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% Yearly monthly plots OA compounds
df['Month']=df['datetime'].dt.month
year_monthly=df[nr].groupby([df['Year'], df['Month']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
year_monthly.plot.bar(y=nr, color=c_df, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('(Year, Month)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% BP by months
bp = dict(linestyle='-', linewidth=0.6)
mp = dict(marker='o', linewidth=0.6,markeredgecolor='black', markerfacecolor='black')
mdp = dict(color='k',  linewidth=0.6)
wp = dict(linestyle='-', linewidth=0.6)

col = 'SO4'
fig, axs =plt.subplots()
df.boxplot(by=['Month'], column=col, showfliers=False, showmeans=True, ax=axs, boxprops=bp, whiskerprops=wp,
           medianprops = mdp, meanprops=mp, fontsize=14)
axs.set_xlabel('Month', fontsize=15)
plt.suptitle(col, fontsize=16)
plt.title('')

#%% BP by years
fig, axs =plt.subplots()
col='Org'
df.boxplot(by=['Year'], column=col, showfliers=False, showmeans=True, ax=axs, boxprops=bp, whiskerprops=wp,
           medianprops = mdp, meanprops=mp, fontsize=12)
plt.suptitle(col, fontsize=14)
plt.title('')
axs.set_xlabel('Year', fontsize=13)

#%% Time series""
oa=pd.read_csv('OA_factors_TS.txt', sep='\t')
del oa['LO-OOA + MO-OOA']
oa['datetime']=pd.to_datetime(oa['PMF_Time_Stop'], dayfirst=True)
oa.index=oa['datetime']
oa=oa[oa.index >pd.to_datetime('01/01/2017 00:00', dayfirst=True)]
#%%
fig, axs=plt.subplots(figsize=(20,4))
oa.plot(y=factors, ax=axs, color=c_oa)
axs.grid('x')
axs.legend(loc='upper left')
#%% Replotting
fig,axs=plt.subplots(figsize=(15,4))
oa.plot(y=factors, ax=axs, color=c_oa,subplots=True)
axs.grid('x')
axs.legend(loc='upper left')
#%% Season plots NR compounds
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
df['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in df.index]
season=df[nr].groupby([df['Season']]).mean()

fig, axs=plt.subplots(figsize=(4,4))
season.plot.bar(y=nr, color=c_df, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('Season', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% Pie NRPM1
df_mean = df.mean()
fig, axs=plt.subplots(figsize=(20,4))
year_pie=df[nr].groupby([df['Year']]).mean()
year_pie.T.plot.pie(subplots=True, legend=False,labeldistance=.5,  pctdistance=1.25, 
                    fontsize=9, ax=axs, colors=c_df, autopct='%1.0f%%', startangle=90,
                    counterclock=False,)
#%% Yearly monthly plots NRPM1 compounds
df['Month']=df['datetime'].dt.month
monthly_year=df[nr].groupby([df['Month'],df['Year']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
monthly_year.plot.bar(y=nr, color=c_df, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('(Month, Year)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% Monthly plots OA compounds
monthly=df[nr].groupby([df['Month']]).mean()
fig, axs=plt.subplots(figsize=(8,4))
monthly.plot.bar(y=nr, color=c_df, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%%
variable='Org'
ym=pd.DataFrame(year_monthly['Org'])
monthly_years=ym.sort_values(by=[ 'Year','Month'])
my=pd.DataFrame({'2014':[],  '2015':})
#%%+++++++++++++++++++++++++++++++++++ OA SA +++++++++++++++++++++++++++++++++++++++++++++++++++++++
#%% OA reconstruction
oa_sum=oa.sum(axis=1)[:-1]
oa_rec=pd.DataFrame()
oa_rec['OA']=org_avg[0]
oa_sum.index=range(0,len(oa_sum))
oa_rec['OA_app']=oa_sum
oa_rec.index=oa.index[:-1]
oa_rec['dt']=oa_rec.index
oa_rec['Year']=oa_rec['dt'].dt.year
fig, axs=plt.subplots(figsize=(5,5))
a=axs.scatter(x=org_avg, y=oa_sum,  c=oa_rec['Year'], vmin=2017)
axs.set_ylabel('OA ($\mu g·m^{-3}$)', fontsize=14)
axs.set_xlabel('OA apportionment ($\mu g·m^{-3}$)', fontsize=14)
m, n = slope(oa_rec['OA'], oa_rec['OA_app'])
cb = plt.colorbar(a)
r=str(R2(oa_rec['OA'], oa_rec['OA_app']))
axs.text(x=3,y=45, s='R$^2$ = '+r+'\n'+'y='+str(m)+' · x + '+str(n))
axs.set_xlim([0,50])
axs.set_ylim([0,50])

#%% Diel yearly plots OA compounds
oa['Hour']=oa['datetime'].dt.hour
oa['Year']=oa['datetime'].dt.year
year_diel=oa[factors].groupby([oa['Year'], oa['Hour']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
year_diel.plot(y=factors, color=c_oa, ax=axs)
axs.grid('x')
axs.set_xticks(range(0,168,12))
# axs.set_xticklabels([0,'2017\n00:00','2018\n00:00','2019\n00:00','2020\n00:00','2021\n00:00','2022\n00:00','2023\n00:00','',''])
axs.set_xticklabels(['2017\n00:00','', '2018\n00:00','','2019\n00:00','','2020\n00:00','','2021\n00:00','',
                     '2022\n00:00','','2023\n00:00','',])
axs.set_xlabel('Year - Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% Pie
oa_mean = oa.mean()
fig, axs=plt.subplots(figsize=(20,4))
year_pie=oa[factors].groupby([oa['Year']]).mean()
year_pie.T.plot.pie(subplots=True, legend=False,labeldistance=.5,  pctdistance=1.25, 
                    fontsize=9, ax=axs, colors=c_oa, autopct='%1.0f%%', startangle=90,
                    counterclock=False,)
axs.set_title('2017')
#%% Yearly monthly plots OA compounds
oa['Month']=oa['datetime'].dt.month
year_monthly=oa[factors].groupby([oa['Year'], oa['Month']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
year_monthly.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('(Year, Month)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%%
#%% Monthly plots OA compounds
monthly=oa[factors].groupby([oa['Month']]).mean()
fig, axs=plt.subplots(figsize=(8,4))
monthly.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%% Season plots OA compounds
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
oa['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in oa.index]
season=oa[factors].groupby([oa['Season']]).mean()

fig, axs=plt.subplots(figsize=(4,4))
season.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('Season', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%%
year_season=oa[factors].groupby([ oa['Season'], oa['Year'],]).mean()

fig, axs=plt.subplots(figsize=(8,4))
year_season.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('(Season, Year)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)

#%%%Inspection of the PMF files!
s=pd.read_csv('Specs.txt', sep='\t', header=None)
e=pd.read_csv('Errors.txt', sep='\t', header=None)
t=pd.read_csv('acsm_utc_time.txt', sep='\t')
mz=pd.read_csv('amus.txt', sep='\t')
s.columns=mz['mz']
e.columns=mz['mz']
s.index=pd.to_datetime(t['acsm_utc_end_time'], dayfirst=True)
e.index=pd.to_datetime(t['acsm_utc_end_time'], dayfirst=True)
#%%
ratios=pd.DataFrame()
ratios['mz57mz55']=s[57]/s[55]
ratios['mz43mz44']=s[43]/s[44]
ratios['mz57mz44']=s[57]/s[44]
ratios['mz60mz44']=s[60]/s[44]
ratios['mz58mz44']=s[58]/s[44]
#%%Filtering
for i in s.columns:
    s[i][s[i]> 1] = np.nan
    e[i][e[i]> 1] = np.nan
    s[i][s[i]< -0.05] = np.nan

#%%
#********************************* BOUNDARY LAYER ******************************************
#%%
path_bl="D:\Data co-located instruments\Boundary Layer"
os.chdir(path_bl)
bl_df=pd.read_csv('BL_Barcelona_2015-2020.txt', sep='\t',dtype = {'Datetime': str, 'BL': float})
bl_df['dt']=pd.to_datetime(bl_df['Datetime'], dayfirst=True)
bl_df.index=bl_df['dt']
# bl_df['BL'].hist(bins=1000)
bl_df.boxplot(column='BL')
#%% BL TS!
fig, axs=plt.subplots(figsize=(15,4))
bl_df['BL'].plot()
axs.grid('x')
#%% Monthly, Yearly, Hourly plots OA compounds
bl_df['Month']=bl_df['dt'].dt.month
bl_df['Year']=bl_df['dt'].dt.year
bl_df['Hour']=bl_df['dt'].dt.hour

fig, axs=plt.subplots(figsize=(10,4))
bl_df.boxplot(column='BL', by='Month', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_title('')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('BL height (m)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
bl_df.boxplot(column='BL', by='Year', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_title('')
axs.set_xlabel('Year', fontsize=12)
axs.set_ylabel('BL height (m)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
bl_df.boxplot(column='BL', by='Hour', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_title('')
axs.set_xlabel('Hour', fontsize=12)
axs.set_ylabel('BL height (m)', fontsize=12)
#%% Averaging to original timestamps
bl_avg=trt.averaging(bl_df['BL'], df['datetime'], bl_df['dt'])
bl_avg.index=dt.iloc[1:]
bl_avg.columns=['BL']
bl_avg.plot()
#%% Normalisation by BL
bl_avg['BL_norm']=bl_avg['BL']/bl_avg['BL'].median()
df_norm=pd.DataFrame()
df_norm['Org']=df['Org']*bl_avg['BL_norm']
df_norm['SO4']=df['SO4']*bl_avg['BL_norm']
df_norm['NO3']=df['NO3']*bl_avg['BL_norm']
df_norm['NH4']=df['NH4']*bl_avg['BL_norm']
df_norm['Chl']=df['Chl']*bl_avg['BL_norm']
df_norm['BC']=df['BC']*bl_avg['BL_norm']
df_norm['BClf']=df['BClf']*bl_avg['BL_norm']
df_norm['BCsf']=df['BCsf']*bl_avg['BL_norm']

#%%
fig, axs=plt.subplots(figsize=(20,4))
df_norm.plot(y=nr, ax=axs, color=c_df)
axs.grid('x')
axs.legend(loc='upper left')
#%%
df_norm['dt']=df['datetime']
df_norm['Month']=df_norm['dt'].dt.month
df_norm['Year']=df_norm['dt'].dt.year
df_norm['Hour']=df_norm['dt'].dt.hour
col='BC'
fig, axs=plt.subplots(figsize=(10,4))
df_norm.boxplot(column=col, by='Month', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
df_norm.boxplot(column=col, by='Year', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Year', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
df_norm.boxplot(column=col, by='Hour', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#%%
df['dt']=df['datetime']
df['Month']=df['dt'].dt.month
df['Year']=df['dt'].dt.year
df['Hour']=df['dt'].dt.hour
col='BC'
fig, axs=plt.subplots(figsize=(10,4))
df.boxplot(column=col, by='Month', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
df.boxplot(column=col, by='Year', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Year', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
df.boxplot(column=col, by='Hour', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#
#%%+++++++++++++++++++++++++++++++++++ OA SA BL +++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  to OA time series timestamps
bl_avg_oa=trt.averaging(bl_df['BL'], oa['datetime'], bl_df['dt'])
bl_avg_oa.index=oa['datetime'].iloc[1:]
bl_avg_oa.columns=['BL']
bl_avg_oa.plot()
#%% Normalisation by BL
bl_avg_oa['BL_norm']=bl_avg_oa['BL']/bl_avg_oa['BL'].median()
oa_norm=pd.DataFrame()
oa_norm['COA']=oa['COA']*bl_avg_oa['BL_norm']
oa_norm['HOA']=oa['HOA']*bl_avg_oa['BL_norm']
oa_norm['Amine-OA']=oa['Amine-OA']*bl_avg_oa['BL_norm']
oa_norm['MO-OOA']=oa['MO-OOA']*bl_avg_oa['BL_norm']
oa_norm['LO-OOA']=oa['LO-OOA']*bl_avg_oa['BL_norm']
oa_norm['BBOA']=oa['BBOA']*bl_avg_oa['BL_norm']
#%%
oa_norm['dt']=oa['datetime']
oa_norm['Month']=oa_norm['dt'].dt.month
oa_norm['Year']=oa_norm['dt'].dt.year
oa_norm['Hour']=oa_norm['dt'].dt.hour
col='MO-OOA'
fig, axs=plt.subplots(figsize=(10,4))
oa_norm.boxplot(column=col, by='Month', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
oa_norm.boxplot(column=col, by='Year', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Year', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
fig, axs=plt.subplots(figsize=(10,4))
oa_norm.boxplot(column=col, by='Hour', showfliers=False, showmeans=True, ax=axs)
axs.grid('x')
axs.set_xlabel('Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
#



#%%#%% Mann-Kendall (Accepted: there is no trend. REjected: tehre is a significant trend.)
#%%
mk=pd.DataFrame()
mk['HOA']=oa['HOA']
mk['dt']=oa['datetime'].dt.date
mk_daily=mk['HOA'].groupby(mk['dt']).mean()

print(trt.MannKendall(mk_daily, 0.05))


#%%
#Check this!!!
# Stacked barplot in Period A and B for NR+BC
fw_tt=fw_t.transpose()
fw_tt=fw_tt.round(1)
OA=fw_tt['Org']+fw_tt['SO4']+fw_tt['NO3']+fw_tt['NH4']+fw_tt['Chl']+fw_tt['BC']
OA=OA.round(1)
bars=fw_tt.plot.bar(stacked=True,figsize=(8.5, 5), grid=True, color = ['lawngreen','red', 'blue', 'gold', 'fuchsia','black'])
for i in range(0,len(OA)):
    plt.text(i-0.15, OA.iloc[i]+0.05, OA.iloc[i], weight='bold',fontsize=13)
ac=0
for j in range(0,len(fw_tt)):
    ac=0
    ac=ac+fw_tt.Org.iloc[j]
    plt.text(j-0.15, ac-1.2, fw_tt.Org.iloc[j],color="white",weight='bold')
    ac=ac+fw_tt.SO4.iloc[j]
    plt.text(j-0.15, ac-0.8, fw_tt.SO4.iloc[j],color="white",weight='bold')
    ac=ac+fw_tt.NO3.iloc[j]
    plt.text(j-0.15, ac-0.7, fw_tt.NO3.iloc[j],color="white",weight='bold')
    ac=ac+fw_tt.NH4.iloc[j]
    plt.text(j-0.15, ac-0.7, fw_tt.NH4.iloc[j],color="white",weight='bold')
    ac=ac+fw_tt.Chl.iloc[j]
    plt.text(j-0.15, ac-0.2,fw_tt.Chl.iloc[j],color="white", weight='bold')
    ac=ac+fw_tt.BC.iloc[j]
    plt.text(j-0.15, ac-0.8,fw_tt.BC.iloc[j],color="white", weight='bold')    
plt.legend(fontsize=9.5)
plt.show()

#%%
bp = dict(linestyle='-', linewidth=1, color='k')
mdp = dict(linestyle='-', linewidth=1.5, color='darkgrey')
mp = dict(marker='o',linewidth=1, markeredgecolor='black', markerfacecolor='k')
wp = dict(linestyle='-', linewidth=1, color='k')
#%%

"""INTERCOMP"""

#%% Importing data
path="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/All Series/Intercomp/"
os.chdir(path)
nr=pd.read_csv('Chemistry_PM1_clean.txt', sep='\t',dayfirst=True)
bc=pd.read_csv('BC_2019_2023.txt', sep='\t', dayfirst=True)
pm1=pd.read_csv('PM1_2019_2023.txt', sep='\t', dayfirst=True)
#%%Averaging
bc['datetime']=pd.to_datetime(bc['date'], dayfirst=True) 
bc.index=bc['datetime']
bc_h=bc.resample('H').mean()

nr['datetime'] = pd.to_datetime(nr['acsm_utc_time'], dayfirst=True)
nr.index=nr['datetime']
nr_h=nr.resample('H').mean()

pm1['datetime'] = pd.to_datetime(pm1['Horarios'], dayfirst=True)
pm1.index=pm1['datetime']

time_range=pd.date_range(start="01/01/2019",end="31/12/2022", freq='1H' )       
#%% Homogenising times
nr_ac, bc_ac, pm1_ac = [],[],[]

for i in range(len(time_range)):
    timestamp = time_range[i]
    if timestamp in nr_h.index:
        row = nr_h.loc[timestamp]
        nr_ac.append(row)
    else:
        nr_ac.append(pd.Series([np.nan]*len(nr_h.columns)))
    if timestamp in bc_h.index:
        row = bc_h.loc[timestamp]
        bc_ac.append(row)
    else:
        bc_ac.append(pd.Series([np.nan]*len(bc_h.columns)))
    if timestamp in pm1.index:
        row = pm1.loc[timestamp]
        pm1_ac.append(row)
    else:
        pm1_ac.append(pd.Series([np.nan]*len(pm1.columns)))

nr_ac2=pd.DataFrame(nr_ac)
nr_ac2.index=time_range
nr_ac2.drop(['Hour',  'BC'], inplace=True, axis=1)
bc_ac2=pd.DataFrame(bc_ac)
bc_ac2.index=time_range
bc_ac2.drop([0], axis=1, inplace=True)
bc_ac2.columns=['BC']
pm1_ac2=pd.DataFrame(pm1_ac)
pm1_ac2.index=time_range
#%% Sum of ACSM + BC
rpm1 = pd.DataFrame()
rpm1.index = time_range
rpm1=pd.concat([rpm1, nr_ac2, bc_ac2, pm1_ac2], axis=1)
rpm1['ACSM_BC']=rpm1['Chl']+rpm1['NH4']+rpm1['NO3']+rpm1['SO4']+rpm1['Org']+rpm1['BC']
#%% ScatterPlot
fig, axs=plt.subplots(figsize=(6,6))
axs.scatter(x=rpm1['PM1'], y=rpm1['ACSM_BC'], c=time_range)
axs.set_xlabel('PM$_1$ GRIMM ($\mu g·m^{-3}$)', fontsize=16)
axs.set_ylabel('PM$_1$ ACSM + BC ($\mu g·m^{-3}$)', fontsize=16)
slope= trt.slope(rpm1['PM1'], rpm1['ACSM_BC'])[0]
interc=trt.slope(rpm1['ACSM_BC'],rpm1['PM1'])[1]
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
axs.text(x=2,y=53, s="y = "+str(slope)+'x + '+str(interc) + '\nR$^2$ = '+str(trt.R2(rpm1['PM1'], rpm1['ACSM_BC'])), fontsize=16)
axs.grid()
axs.set_ylim(-1,65)
axs.set_xlim(-1,65)
#%% Time series plot
fig2, axs2=plt.subplots(figsize=(12,3))
rpm1['datetime_range']=time_range
rpm1.plot(x='datetime_range', y=['PM1', 'ACSM_BC'], ax=axs2, alpha=0.6)
axs2.set_ylabel('PM$_1$ ($\mu g·m^{-3}$)', fontsize=16)
axs2.set_xlabel('Date', fontsize=16)
slope= trt.slope(rpm1['PM1'], rpm1['ACSM_BC'])[0]
interc=trt.slope(rpm1['PM1'], rpm1['ACSM_BC'])[1]
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
axs2.grid()
axs2.set_ylim(-1,70)
axs2.legend(['PM$_1$ GRIMM', 'PM$_1$ ACSM + BC'], loc='upper right')
rpm1.to_csv("Hourly_data.txt", sep='\t')
#%% Years
years=[2019, 2020, 2021, 2022]
for i in years:
    print('Year: ', i)
    pm1_year=pd.DataFrame()
    mask = (rpm1['Year']==i)
    pm1_year = rpm1[mask]
    pm1_year.to_csv('PM1_intercomp_'+str(i)+'.txt', sep='\t')
#%% Yearly treatment
pm1_2019=rpm1[(rpm1['Year']==2019)]
nrpm1=pm1_2019
#%%
#*************** IONS ANALYSIS *********************
#
#%% Importing files
path_ions="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/All Series/PMF_Matrices"
os.chdir(path_treatment)
specs=pd.read_csv('Specs.txt', sep='\t', header=None)
specs_time=pd.read_csv('acsm_utc_time.txt', sep='\t')
specs_amus=pd.read_csv('amus.txt', sep='\t')
specs.columns=specs_amus['mz']
dt_specs=pd.to_datetime(specs_time['acsm_utc_end_time'], dayfirst=True)
specs['dt_specs']=dt_specs
specs.index=dt_specs
#Ions definitions
oa=specs.sum(axis=1)
ions=pd.DataFrame()
ions['f44']=specs[44]/oa
ions['f43']=specs[43]/oa
ions['f60']=specs[60]/oa
ions['f55']=specs[55]/oa
ions['f57']=specs[57]/oa
ions['f73']=specs[73]/oa
ions['dt']=specs.index
ions_labels=['f44', 'f43', 'f60', 'f55', 'f57', 'f73']
#%%filtering ions 
ions.plot()
mask = (ions['f44']>=0.0) & (ions['f44']<=1.0) & (ions['f43']>=0.0) & (ions['f43']<=1.0) & (ions['f73']>=0.0) & (ions['f60']>=0.0) & (ions['f55']>=0.0) & (ions['f57']>=0.0)
ions_f = ions.loc[mask]
ions_f.plot()
ions_f['Year']=ions['dt'].dt.year
ions_f['Month']=ions['dt'].dt.month
ions_f.to_csv('Ions.txt', sep='\t')
ions_ym=pd.DataFrame(ions_f[ions_labels].groupby([ions_f['Year'], ions_f['Month'],]).mean())
#%%Rearrange ions yearly month
li, li_std=[],[]
yearmonth=pd.DataFrame(pd.date_range(start = '01/01/2014', end = "07/01/2023").to_period('M').unique(), columns=['ym'])
for i in range(0,len(yearmonth)):
    year = yearmonth['ym'].iloc[i].year
    month = yearmonth['ym'].iloc[i].month
    print(year, month)
    mask = (ions_f['Year']==year) & ((ions_f['Month']==month))
    prova = ions_f.loc[mask]
    li.append(prova.mean())
    li_std.append(prova.std())
ym_mean=pd.DataFrame(li)
ym_std=pd.DataFrame(li_std)
#%%  Plotting year/months ions
ym_mean.index, ym_std.index = yearmonth['ym'], yearmonth['ym']
os.chdir(path_treatment)
ion=['f44', 'f43', 'f55', 'f57', 'f60', 'f73']  
fig, axs=plt.subplots(figsize=(10,10))
pl=ym_mean[ion].plot(ax=axs, x_compat=True, color='grey', subplots=True, sharex=True, lw=2, marker='o', grid=True)   
axs.set_xlabel('Monthly means', fontsize=11)
plt.savefig('Ions_monthyear.png', bbox_inches='tight')
#%% Ions Year intercomp
cols=ym_mean.columns[0:6]
ym_mean['6073']=ym_mean['f73'] + ym_mean['f60']
ym_mean['5557']=ym_mean['f55'] + ym_mean['f57']
ym_mean['4344']=ym_mean['f43'] + ym_mean['f44']
cols=['f44','f43', 'f60','f55', 'f57', 'f73','6073', '5557', '4344']
j=0
toplot=cols[j]
titles_cols=['f44', 'f43', 'f60', 'f55', 'f57', 'f73', 'f60 + f73', 'f55 + f57', 'f43 + f44']
yearly_ions=pd.DataFrame()
liy=[]
ym_mean=ym_mean.reset_index(drop=True)
years=['2014', '2015','2016', '2017', '2018', '2019', '2020','2021', '2022', '2023']
for i in range(0,len(years)):
    print(i,i*12+1,(i+1)*12+1)
    liy.append(ym_mean[toplot].iloc[i*12+1:(i+1)*12])
    df_temp=pd.DataFrame(ym_mean[toplot].iloc[i*12+1:(i+1)*12+1])
    yearly_ions=pd.concat([yearly_ions.reset_index(drop=True), df_temp.reset_index(drop=True)], axis=1, ignore_index=True)
yearly_ions.columns=years
months=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
yearly_ions.index=months
# Yearly plot
greys=['gainsboro','lightgrey','lightgrey','silver', 'darkgrey',  '#8F8E8E', 'grey' , 'dimgrey', '#3C3C3C', 'k']
rainbow=['mediumpurple', 'hotpink', 'r', 'orange', 'gold', 'yellowgreen', 'green', 'skyblue', 'dodgerblue', 'slateblue']
markers=['o', 'v','X','s','^', 'p', '>', '*', '<', 'P']
fig, axs=plt.subplots(figsize=(8,4))
# pl=yearly_ions[years[9]].plot.line()
for i in range(0, len(yearly_ions.columns)):
    pl=yearly_ions[years[i]].plot.line( fontsize=13,ax=axs, color=greys[i], grid=True, markersize=8, marker=markers[i], legend=True, lw=2)
axs.set_xticks(range(len(months)))
axs.set_xticklabels(months, fontsize=14)
axs.legend(bbox_to_anchor=(1.005, 1.005))
axs.set_ylabel(titles_cols[j]+'\n (adim.)', fontsize=15)
plt.savefig(toplot+'_monthyear.png', bbox_inches='tight')
#%% Ratios definition
ratios=pd.DataFrame()
ratios['dt']=ions_f['dt']
ratios['SOA_freshness']=ions['f43']/ions['f44']
ratios['POA_SOA']=(ions['f55']+ions['f57']+ions['f60']+ions['f73'])/(ions['f43']+ions['f44'])
ratios['Traffic']=(ions['f55']+ions['f57'])/(ions['f43']+ions['f44']+ions['f60']+ions['f73']+ions['f55']+ions['f57'])
ratios['BB']=(ions['f60']+ions['f73'])/(ions['f43']+ions['f44']+ions['f60']+ions['f73']+ions['f55']+ions['f57'])
ratios['OC']=0.079+4.31*ions['f44']
ratios['OAOC']=1.29*ratios['OC']+1.17
ratios['OAOC_nop']=oa/ratios['OC']
ratios['OA']=oa
mask_ratios = (ratios['SOA_freshness']>=0.0) & (ratios['SOA_freshness']<=2.0) & (ratios['POA_SOA']>=0.0) #& (ratios['OAOC_nop']<=20.0)
ratios_f=ratios[mask_ratios]
ratios_f[ratios_f.columns[1:8]].plot(figsize=(12,12), legend=True, subplots=True, lw=2, color='grey')
ratios_f.to_csv('ratios.txt', sep='\t')
#%% Ratios averages year/months
li, li_std=[],[]
ratios_f['Month'], ratios_f['Year']=ratios_f['dt'].dt.month, ratios_f['dt'].dt.year
yearmonth=pd.DataFrame(pd.date_range(start = '01/01/2014', end = "07/01/2023").to_period('M').unique(), columns=['ym'])
for i in range(0,len(yearmonth)):
    year = yearmonth['ym'].iloc[i].year
    month = yearmonth['ym'].iloc[i].month
    print(year, month)
    mask = (ratios_f['Year']==year) & ((ratios_f['Month']==month))
    prova = ratios_f.loc[mask]
    li.append(prova.mean())
    li_std.append(prova.std())
ym_ratios_mean=pd.DataFrame(li)
ym_ratios_std=pd.DataFrame(li_std)  
#%%  Plotting year/months ions
ym_ratios_mean.index, ym_ratios_std.index = yearmonth['ym'], yearmonth['ym']
os.chdir(path_treatment)
toplot=['SOA_freshness', 'POA_SOA', 'Traffic', 'BB', 'OC', 'OAOC']  
fig, axs=plt.subplots(figsize=(10,10))
pl=ym_ratios_mean[toplot].plot(ax=axs, x_compat=True, color='grey', subplots=True, sharex=True, lw=2, marker='o', grid=True)   
axs.set_xlabel('Monthly means', fontsize=11)
plt.savefig('IonRatios_monthyear.png', bbox_inches='tight')
#%%Studying further OA:OC
fig, axs=plt.subplots(figsize=(12,2))
ym_ratios_mean['OAOC'].plot(ax=axs, color='darkgreen', marker='o', grid=True)
axs.set_ylabel('Parametrised OA:OC\n(adim.)')
axs.set_ylim(1.5,2.5)
#%%Studying further Others
fig, axs=plt.subplots(figsize=(12,2))
ym_ratios_mean['POA_SOA'].plot(ax=axs, color='violet', marker='o', grid=True)
axs.set_ylabel('POA/SOA ($\mu g·m^{-3}$.)', fontsize=12)
# axs.set_ylim(1.5,2.5)
axs.set_xlabel('Monthly means')
#%% Year intercomp
cols=ym_ratios_mean.columns[0:6]
j=0
toplot=cols[j]
titles_cols=['SOA freshness', 'POA / SOA', 'Traffic', 'Biomass burning', 'OC', 'OA / OC']
yearly=pd.DataFrame()
liy=[]
ym_ratios_mean=ym_ratios_mean.reset_index(drop=True)
years=['2014', '2015','2016', '2017', '2018', '2019', '2020','2021', '2022', '2023']
for i in range(0,len(years)):
    print(i,i*12+1,(i+1)*12+1)
    liy.append(ym_ratios_mean[toplot].iloc[i*12+1:(i+1)*12])
    df_temp=pd.DataFrame(ym_ratios_mean[toplot].iloc[i*12+1:(i+1)*12+1])
    yearly=pd.concat([yearly.reset_index(drop=True), df_temp.reset_index(drop=True)], axis=1, ignore_index=True)
yearly.columns=years
months=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
yearly.index=months
# Yearly plot
greys=['gainsboro','lightgrey','lightgrey','silver', 'darkgrey',  '#8F8E8E', 'grey' , 'dimgrey', '#3C3C3C', 'k']
rainbow=['mediumpurple', 'hotpink', 'r', 'orange', 'gold', 'yellowgreen', 'green', 'skyblue', 'dodgerblue', 'slateblue']
fig, axs=plt.subplots(figsize=(8,4))
markers=['o', 'v','X','s','^', 'p', '>', '*', '<', 'P']
for i in range(0,len(years)):
    pl=yearly.plot.line(y=[years[i]], fontsize=13,ax=axs, color=greys[i], grid=True, markersize=8, marker=markers[i], legend=True, lw=2)
axs.set_xticks(range(len(yearly)))
axs.set_xticklabels(months, fontsize=14)
axs.legend(bbox_to_anchor=(1.005, 1.005))
axs.set_ylabel(titles_cols[j]+'\n (adim.)', fontsize=15)
plt.savefig(toplot+'_monthyear.png', bbox_inches='tight')
yearly_mean=yearly.mean(axis=1)
fig2, axs2=plt.subplots(figsize=(8,4))
yearly_mean.plot(color='grey', title=titles_cols[j], ax=axs2)
plt.savefig(toplot+'_average.png', bbox_inches='tight')













    
    
    
    
    


