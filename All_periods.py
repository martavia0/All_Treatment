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
factors = ['HOA', 'COA', 'Amine-OA', 'BBOA', 'LO-OOA', 'MO-OOA' ]
ions = ['mz57mz55', 'mz43mz44', 'mz57mz44','mz60mz44', 'mz58mz44' ]
c_df=['green', 'red', 'blue', 'orange', 'fuchsia', 'grey', 'saddlebrown'] #colors_nrpm1
#%%
path_treatment="C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/All Series/Treatment"
os.chdir(path_treatment)
df=pd.read_csv('Chemistry_PM1.txt', sep='\t')
#%% Fixing time
dt=pd.to_datetime(df['acsm_utc_end_time'], dayfirst=True)
df.index=dt
df['datetime']=df.index
#%% Overall time series
fig, axs=plt.subplots(figsize=(20,4))
df.plot(y=nr, ax=axs, color=c_df)
axs.grid('x')
axs.legend(loc='upper left')
#%% Replacing the negatives below their DL (0.148, 0.0224, 0.012, 0.284, 0.011)
df.Org[df.Org < -0.148] = np.nan
df.SO4[df.SO4 < -0.0224] = np.nan
df.NO3[df.NO3 < -0.012] = np.nan
df.NH4[df.NH4 < -0.284] = np.nan
df.Chl[df.Chl < -0.011] = np.nan
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

col = 'NO3'
fig, axs =plt.subplots()
df.boxplot(by=['Month'], column=col, showfliers=False, showmeans=True, ax=axs, boxprops=bp, whiskerprops=wp,
           medianprops = mdp, meanprops=mp, fontsize=14)
axs.set_xlabel('Month', fontsize=15)
plt.suptitle(col, fontsize=16)
plt.title('')

#%% BP by years
fig, axs =plt.subplots()
col='Chl'
df.boxplot(by=['Year'], column=col, showfliers=False, showmeans=True, ax=axs, boxprops=bp, whiskerprops=wp,
           medianprops = mdp, meanprops=mp, fontsize=12)
plt.suptitle(col, fontsize=14)
plt.title('')
axs.set_xlabel('Year', fontsize=13)

#%% Season plots NR compounds
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
df['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in df.index]
season=df[nr].groupby([df['Season']]).mean()

fig, axs=plt.subplots(figsize=(4,4))
season.plot.bar(y=nr, color=c_df, ax=axs, stacked=True, zorder=3)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Season', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.63))
#%% Season yearly NR compounds
season_year=df.groupby(by=['Year', 'Season']).mean()

fig, axs=plt.subplots(figsize=(6,4))
season_year.plot.bar(y=nr, color=c_df, ax=axs, stacked=True, legend=False, zorder=3)
axs.grid(axis='y', zorder=0)
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
#%% SUPERATIONS!!

limit_who_25_daily = 15
dfi=df.copy(deep=True)
dfi['date']=dfi['datetime'].dt.date
dfi=dfi.groupby('date').mean()
dfi['nr']=dfi['Org']+dfi['NO3']+dfi['NH4']+dfi['SO4']+dfi['Chl']
dfi['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in dfi.index]
dfi2 = dfi.copy(deep=True) 
mask=dfi['nr'].iloc[:]>=limit_who_25_daily
dfi=dfi.loc[mask]
#
dfi['Year'] =dfi['Year'].astype(int)
dfi['SY'] = dfi['Year'].astype(str)+dfi['Season']
dfi2['Year'] =dfi2['Year'].astype(int)
dfi2['SY'] = dfi2['Year'].astype(str)+dfi2['Season']

sup_season_year = dfi.groupby(['Year', 'Season']).count()#/nb_days_season
index_ndays=dfi2.groupby(['Year', 'Season']).count()

b= pd.DataFrame()
b['date'] = pd.date_range(start='1/1/2014', end='31/12/2023')
b['Year'] =b['date'].dt.year
b['Month']=b['date'].dt.month
b['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in b['date']]
c = b.groupby(['Year', 'Season']).count()
d = c.index

sup1 = sup_season_year['Org']
sup2 = index_ndays['Org']
# sup3 = pd.merge(left=sup2, right=sup1, left_index=True, right_index=True, how='outer')
# sup3['Perc'] = 100.0*sup3['Org_y']/sup3['Org_x']
sup4 = pd.merge(left=c, right=sup1, left_index=True, right_index=True, how='outer')
sup4=sup4['Org']
sup4=sup4.reset_index()

for i in range(0, len(d)):
    if d[i] not in sup2.index:
        print(i, d[i])
        sup4.loc[i]=999.9
sup4.index=d
sup4=sup4.replace(np.nan, 0.0)
sup4=sup4.replace(999.9, np.nan)

sup5 = pd.merge(left=sup4, right=sup2, left_index=True, right_index=True, how='outer')
sup5['Perc'] = 100.0*sup5['Org_x']/sup5['Org_y']


fig, axs =plt.subplots(figsize=(8,4))
axs.grid(axis ='y', color='lightgrey')
sup5['Perc'].plot(kind='line', marker='D', color='k', zorder=3.5, ax=axs, lw=0, markersize=4)
sup5['Perc'].plot(kind='bar', stacked=True, color='dimgrey', zorder=3, ax=axs, align='center')
# sup5['Perc'].plot(kind='line', marker='D', color='dimgrey', zorder=3, ax=axs, lw=0)

axs.set_ylabel('% of days superating WHO \n daily PM$_{2.5}$ threshold', fontsize=12)
x1, y1 = [3.5, 3.5], [0, 50]
x2, y2 = [7.5, 7.5], [0, 50]
x3, y3 = [11.5, 11.5], [0, 50]
x4, y4 = [15.5, 15.5], [0, 50]
x5, y5 = [19.5, 19.5], [0, 50]
x6, y6 = [23.5, 23.5], [0, 50]
x7, y7 = [27.5, 27.5], [0, 50]
x8, y8 = [31.5, 31.5], [0, 50]
x9, y9 = [35.5, 35.5], [0, 50]


axs.plot(x1, y1, x2, y2, x3,y3, x4,y4, x5,y5, x6,y6,x7,y7, x8,y8, x9,y9, marker = '', color='grey', ls=':')
axs.set_ylim(0,35)
axs.grid(axis ='y', color='lightgrey', zorder=0)
axs.set_title('NR-PM$_1$', fontsize=14)

#%% Superations stacked Season, Year
limit_who_25_daily = 15
dfi=df.copy(deep=True)
dfi['date']=dfi['datetime'].dt.date
dfi=dfi.groupby('date').mean()
dfi['nr']=dfi['Org']+dfi['NO3']+dfi['NH4']+dfi['SO4']+dfi['Chl']
dfi['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in dfi.index]
dfi2 = dfi.copy(deep=True) 
mask=dfi['nr'].iloc[:]>=limit_who_25_daily
dfi=dfi.loc[mask]

fig, axs =plt.subplots(figsize=(8,4))

x1, y1 = [9.5,9.5], [0, 100]
x2, y2 = [19.5,19.5], [0, 105]
x3, y3 = [29.5,29.5], [0, 105]
x4, y4 = [39.5,39.5], [0, 105]
x5, y5 = [49.5,49.5], [0, 105]

axs.set_ylim(0,100)
axs.plot(x1, y1, x2, y2, x3,y3, x4,y4, x5,y5, marker = '', color='k', ls=':')

from matplotlib.patches import Rectangle
axs.add_patch(Rectangle((-0.5, 0), 10,100, facecolor='#DCF3FF', fill=True))
axs.add_patch(Rectangle((9.5,0), 10,100,  facecolor='#FFFFDC', fill=True))
axs.add_patch(Rectangle((19.5,0), 10,100,  facecolor='#E4FFDC', fill=True))
axs.add_patch(Rectangle((29.5,0), 10,100,  facecolor='#F1E0CD', fill=True))

dfi2=dfi.groupby([dfi['Season'], dfi['Year']]).mean()
c = b.groupby(['Season', 'Year']).count()
dfi3 = pd.merge(left=c, right=dfi2, left_index=True, right_index=True, how='outer')
dfi_idx = dfi3.index
dfi.reset_index(inplace=True)
dfi4=pd.DataFrame({'OA':dfi3['Org']*100/dfi3[nr].sum(axis=1),'SO4':dfi3['SO4']*100/dfi3[nr].sum(axis=1),
                    'NO3':100*dfi3['NO3']/dfi3[nr].sum(axis=1), 'NH4':100*dfi3['NH4']/dfi3[nr].sum(axis=1), 
                    'Cl':100*dfi3['Chl']/dfi3[nr].sum(axis=1)})
dfi4.plot(kind='bar', stacked=True, ax=axs, color=c_df)
axs.set_xticks(range(0,40))
axs.set_xticklabels([int(i) for i in range(2014, 2024)]*4)
axs.set_xlabel('\nSeasons, Years', fontsize=12)
axs.set_ylabel('Mean composition ($\mu g·m^{-3}$)', fontsize=12)
axs.set_title('Days with daily WHO limits superation')
axs.legend(loc=[1.01, 0.64])

axs.text(x = 4, y=-22, s='DJF')
axs.text(x = 14, y=-22, s='JJA')
axs.text(x = 23, y=-22, s='MAM')
axs.text(x = 34, y=-22, s='SON')

#%%
'''
+++++++++++++++++++++++++++++++++++ OA SA +++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''
#%% OA SA:  Import and time series""
path_oa = "C:/Users/maria/Documents/Marta Via/1. PhD/A. Data/All Series/Treatment"
oa=pd.read_csv('SA_TS.txt', sep='\t')
oa['datetime']=pd.to_datetime(oa['PMF Time (UTC)'], dayfirst=True)
oa.index=oa['datetime']

c_oa=[ 'grey','mediumpurple','steelblue', 'saddlebrown',  'lightgreen', 'darkgreen'] #colors_factors
factors = ['HOA','COA', 'Amine-OA','BBOA', 'LO-OOA', 'MO-OOA'] #labels_factors

fig, axs=plt.subplots(figsize=(20,4))
oa.plot(y=factors, ax=axs, color=c_oa)
axs.grid('x')
axs.legend(loc='upper left')
#%% OA SA: Simple pie
oa[factors].mean().plot.pie(colors=c_oa, autopct='%1.0f%%')
#%% OA SA:  Replotting
oa['Year']=oa['datetime'].dt.year
oa['Month']=oa['datetime'].dt.month
oa['Hour']=oa['datetime'].dt.hour
y= oa.groupby('Year').mean()
yy=oa.groupby('Year').std()
m=oa.groupby('Month').mean()
mm=oa.groupby('Month').std()
h=oa.groupby('Hour').mean()
hh=oa.groupby('Hour').std()

fig,axs=plt.subplots(figsize=(13,10), nrows=len(factors), sharex='col', ncols=4, width_ratios=[3, 1,1, 1])
for i in range(0,len(factors)):
    axs[i,0].plot(oa[factors[i]], c=c_oa[i])
    axs[i,1].plot(y[factors[i]], c=c_oa[i], marker='o')
    axs[i,1].errorbar(x=y.index,y=y[factors[i]], yerr=yy[factors[i]], color= c_oa[i])
    axs[i,2].plot(m[factors[i]], c=c_oa[i], marker='o')
    axs[i,2].errorbar(x=m.index,y=m[factors[i]], yerr=mm[factors[i]], color= c_oa[i])
    axs[i,3].plot(h[factors[i]], c=c_oa[i] , marker='o')
    axs[i,3].errorbar(x=h.index,y=h[factors[i]], yerr=hh[factors[i]], color= c_oa[i])
    axs[i,0].set_ylabel(factors[i])
axs[0,1].set_xlim(2014,2024)
axs[5,0].set_xlabel('Time series')
axs[5,1].set_xlabel('Years')
axs[5,2].set_xlabel('Months')
axs[5,3].set_xlabel('Hours')

# axs.grid(axis='x')
#%% OA SA:  OA reconstruction
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

#%% OA SA:  Diel yearly plots OA compounds
oa['Hour']=oa['datetime'].dt.hour
oa['Year']=oa['datetime'].dt.year
year_diel=oa[factors].groupby([oa['Year'], oa['Hour']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
year_diel.plot(y=factors, color=c_oa, ax=axs)
axs.grid('x')
axs.set_xticks(range(0,210,12))
# axs.set_xticklabels([0,'2017\n00:00','2018\n00:00','2019\n00:00','2020\n00:00','2021\n00:00','2022\n00:00','2023\n00:00','',''])
axs.set_xticklabels(['2014\n00:00','', '2015\n00:00','','2017\n00:00','', '2018\n00:00','','2019\n00:00','','2020\n00:00','','2021\n00:00','',
                     '2022\n00:00','','2023\n00:00','',])
axs.set_xlabel('Year - Hour', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.56))
#%% OA SA:  Pie by years
oa_mean = oa.mean()
fig, axs=plt.subplots(figsize=(22,4))
year_pie=oa[factors].groupby([oa['Year']]).mean()
year_pie.T.plot.pie(subplots=True, legend=False,  pctdistance=.75, 
                    fontsize=9, ax=axs, colors=c_oa, autopct='%1.0f%%', startangle=90,
                    labels=None,counterclock=False)
#%% OA SA:  Stacked yearly
oa_mean = oa.mean()
fig, axs=plt.subplots(figsize=(22,4))
fig, axs= plt.subplots(figsize=(6,6))
y[factors].plot(kind='bar', stacked=True, color=c_oa, ax=axs)
axs.legend(loc=(1.03, 0.70))
axs.set_ylabel('Absolute concentrations ($\mu g·m^{-3}$)')
#%% OA SA:  Yearly monthly plots OA compounds
oa['Month']=oa['datetime'].dt.month
year_monthly=oa[factors].groupby([oa['Year'], oa['Month']]).mean()
fig, axs=plt.subplots(figsize=(10,4))
year_monthly.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True)
axs.grid('x')
axs.set_xlabel('(Year, Month)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)

#%% OA SA:  Monthly plots OA compounds
monthly=oa[factors].groupby([oa['Month']]).mean()
fig, axs=plt.subplots(figsize=(5,4))
monthly.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True, zorder=2)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Month', fontsize=12)
axs.set_ylabel('Absolute concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

#%% OA SA:  Season plots OA compounds
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
oa['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in oa.index]
season=oa[factors].groupby([oa['Season']]).mean()

fig, axs=plt.subplots(figsize=(4,4))
season.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True, zorder=2)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Season', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

#%% OA SA:  Yearly season
rdm_dates = pd.DataFrame()
rdm_dates['datetime']=pd.date_range(start='01/01/2014', end='31/12/2023')
rdm_dates['Month']=rdm_dates['datetime'].dt.month
rdm_dates['Year']=rdm_dates['datetime'].dt.year
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
rdm_dates['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in rdm_dates['datetime']]
ys_all = rdm_dates.groupby( [rdm_dates['Year'], rdm_dates['Season']]).count()

year_season=oa[factors].groupby([ oa['Year'],oa['Season']]).mean()
ys=pd.DataFrame()
ys = pd.merge(left=ys_all, right=year_season, left_index=True, right_index=True, how='outer')

fig, axs=plt.subplots(figsize=(12,4))
ys.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True, zorder=3)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('(Year, Season)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

x1, y1 = [3.5,3.5], [0, 50]
x2, y2 = [7.5,7.5], [0, 50]
x3, y3 = [11.5,11.5], [0, 50]
x4, y4 = [15.5, 15.5], [0, 50]
x5, y5 = [19.5, 19.5], [0, 50]
x6, y6 = [23.5, 23.5], [0, 50]
x7, y7 = [27.5, 27.5], [0, 50]
x8, y8 = [31.5, 31.5], [0, 50]
x9, y9 = [35.5, 35.5], [0, 50]

axs.set_ylim(0,8)
axs.plot(x1, y1, x2, y2, x3,y3, x4,y4, x5,y5, x6,y6,x7,y7, x8,y8, x9,y9, marker = '', color='grey', ls=':')
#%% Season yearly
rdm_dates = pd.DataFrame()
rdm_dates['datetime']=pd.date_range(start='01/01/2014', end='31/12/2023')
rdm_dates['Month']=rdm_dates['datetime'].dt.month
rdm_dates['Year']=rdm_dates['datetime'].dt.year
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
rdm_dates['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in rdm_dates['datetime']]
ys_all = rdm_dates.groupby( [rdm_dates['Season'], rdm_dates['Year']]).count()

year_season=oa[factors].groupby([ oa['Season'],oa['Year']]).mean()
ys=pd.DataFrame()
ys = pd.merge(left=ys_all, right=year_season, left_index=True, right_index=True, how='outer')

fig, axs=plt.subplots(figsize=(12,4))
ys.plot.bar(y=factors, color=c_oa, ax=axs, stacked=True, zorder=3)
axs.grid(axis='y', zorder=1)
axs.set_xlabel('(Year, Season)', fontsize=12)
axs.set_ylabel('Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

x1, y1 = [9.5,9.5], [0, 50]
x2, y2 = [19.5, 19.5], [0, 50]
x3, y3 = [29.5, 29.5], [0, 50]

from matplotlib.patches import Rectangle
axs.add_patch(Rectangle((-0.5, 0), 10, 8, facecolor='#DCF3FF', fill=True))
axs.add_patch(Rectangle((9.5,0), 10,10,  facecolor='#FFFFDC', fill=True))
axs.add_patch(Rectangle((19.5,0), 10,10,  facecolor='#E4FFDC', fill=True))
axs.add_patch(Rectangle((29.5,0), 10,10,  facecolor='#F1E0CD', fill=True))

#E4FFDC
axs.set_ylim(0,8)
axs.plot(x1, y1, x2, y2, x3,y3, marker = '', color='grey', ls=':')
#%% Monthly yearly
rdm_dates = pd.DataFrame()
rdm_dates['datetime']=pd.date_range(start='01/01/2014', end='31/12/2023')
rdm_dates['Month']=rdm_dates['datetime'].dt.month
rdm_dates['Year']=rdm_dates['datetime'].dt.year
ym_all = rdm_dates.groupby( [rdm_dates['Year'], rdm_dates['Month']]).count()

year_month=oa[factors].groupby([ oa['Year'],oa['Month']]).mean()
ym = pd.merge(left=ym_all, right=year_month, left_index=True, right_index=True, how='outer')
ym =ym.replace(0, np.nan)
# ym =ym.replace(np.na)


fig, axs=plt.subplots(figsize=(10,5))
ym.plot(y=factors, color=c_oa, ax=axs, zorder=3, marker='o', markersize=3)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Year', fontsize=12)
axs.set_ylabel('Monthly Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

axs.set_xticks(range(0,120, 12))
axs.set_xticklabels([int(i) for i in range(2014, 2024)])

x0, y0 = [0.0, 0.0], [0, 30]
x1, y1 = [12, 12], [0, 30]
x2, y2 = [24, 24], [0, 30]
x3, y3 = [36, 36], [0, 30]
x4, y4 = [48, 48], [0, 30]
x5, y5 = [60, 60], [0, 30]
x6, y6 = [72, 72], [0, 30]
x7, y7 = [84, 84], [0, 30]
x8, y8 = [96, 96], [0, 30]
x9, y9 = [108, 108], [0, 30]
x10, y10 = [120, 120], [0, 30]

x11, y11 = [6, 6], [0, 30]
x12, y12 = [18, 18], [0, 30]
x13, y13 = [30, 30], [0, 30]
x14, y14 = [42, 42], [0, 30]
x15, y15 = [54, 54], [0, 30]
x16, y16 = [66, 66], [0, 30]
x17, y17 = [78, 78], [0, 30]
x18, y18 = [90, 90], [0, 30]
x19, y19 = [102, 102], [0, 30]
x20, y20 = [114, 114], [0, 30]

axs.set_ylim(0,4)
axs.plot(x0, y0, x1, y1, x2, y2, x3,y3, x4,y4, x5,y5, x6,y6,x7,y7, x8,y8, x9,y9, x10,y10, marker = '', color='grey', ls=':')
axs.plot(x11, y11, x12, y12, x13,y13, x14,y14, x15,y15, x16,y16, x17,y17, x18,y18, x19,y19, x20,y20, marker = '', color='lightgrey', ls=':')


#%% Yearly monthly
rdm_dates = pd.DataFrame()
rdm_dates['datetime']=pd.date_range(start='01/01/2014', end='31/12/2023')
rdm_dates['Month']=rdm_dates['datetime'].dt.month
rdm_dates['Year']=rdm_dates['datetime'].dt.year
ym_all = rdm_dates.groupby( [rdm_dates['Month'], rdm_dates['Year']]).count()

year_month=oa[factors].groupby([ oa['Month'],oa['Year']]).mean()
ym = pd.merge(left=ym_all, right=year_month, left_index=True, right_index=True, how='outer')
ym =ym.replace(0, np.nan)
# ym =ym.replace(np.na)


fig, axs=plt.subplots(figsize=(12,5))
ym.plot(kind='bar', width=0.9,  y=factors, color=c_oa, ax=axs, zorder=3, stacked=True)#, marker='o', markersize=3, lw=0)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('\nMonth, year', fontsize=12)
axs.set_ylabel('Yearly Concentration ($\mu g·m^{-3}$)', fontsize=12)
axs.legend(loc=(1.03, 0.55))

axs.set_xticks(range(0,120, 2))
axs.set_xticklabels(['2014',  '',  '2018', '', '2022']*12, rotation =90)

x0, y0 = [0.0, 0.0], [0, 30]
x1, y1 = [10.5, 10.5], [0, 30]
x2, y2 = [20.5, 20.5], [0, 30]
x3, y3 = [30.5, 30.5], [0, 30]
x4, y4 = [40.5, 40.5], [0, 30]
x5, y5 = [50.5, 50.5], [0, 30]
x6, y6 = [60.5, 60.5], [0, 30]
x7, y7 = [70.5, 70.5], [0, 30]
x8, y8 = [80.5, 80.5], [0, 30]
x9, y9 = [90.5, 90.5], [0, 30]
x10, y10 = [100.5, 100.5], [0, 30]
x11, y11 = [110.5, 110.5], [0, 30]
x12, y12 = [120.5, 120.5], [0, 30]


axs.set_ylim(0,10)
axs.plot(x0, y0, x1, y1, x2, y2, x3,y3, x4,y4, x5,y5, x6,y6,x7,y7, x8,y8, x9,y9, x10,y10, marker = '', color='k', ls=':')
axs.plot(x11, y11, x12, y12, marker = '', color='k', ls=':')#, x13,y13, x14,y14, x15,y15, x16,y16, x17,y17, x18,y18, x19,y19, x20,y20, marker = '', color='grey', ls=':')
axs.add_patch(Rectangle((0, 0), 10.5, 11, facecolor='#D0F2FF', fill=True))
axs.add_patch(Rectangle((10.5, 0), 10, 11, facecolor='#B7EBFE', fill=True))
axs.add_patch(Rectangle((20.5, 0), 10, 11, facecolor='#DAFFDC', fill=True))
axs.add_patch(Rectangle((30.5, 0), 10, 11, facecolor='#C7FFCA', fill=True))
axs.add_patch(Rectangle((40.5, 0), 10, 11, facecolor='#BAFFBE', fill=True))
axs.add_patch(Rectangle((50.5, 0), 10, 11, facecolor='#FEFFD5', fill=True))
axs.add_patch(Rectangle((60.5, 0), 10, 11, facecolor='#FEFFC3', fill=True))
axs.add_patch(Rectangle((70.5, 0), 10, 11, facecolor='#FEFFB3', fill=True))
axs.add_patch(Rectangle((80.5, 0), 10, 11, facecolor='#FEF3D3', fill=True))
axs.add_patch(Rectangle((90.5, 0), 10, 11, facecolor='#FFF0C5', fill=True))
axs.add_patch(Rectangle((100.5, 0), 10, 11, facecolor='#FEECB6', fill=True))
axs.add_patch(Rectangle((110.5, 0), 10, 11, facecolor='#E5F8FF', fill=True))

axs.text(x=2.5, y=-1.65, s='JAN')
axs.text(x=12.5, y=-1.65, s='FEB')
axs.text(x=22.5, y=-1.65, s='MAR')
axs.text(x=32.5, y=-1.65, s='APR')
axs.text(x=42.5, y=-1.65, s='MAY')
axs.text(x=52.5, y=-1.65, s='JUN')
axs.text(x=62.5, y=-1.65, s='JUL')
axs.text(x=72.5, y=-1.65, s='AUG')
axs.text(x=82.5, y=-1.65, s='SEP')
axs.text(x=92.5, y=-1.65, s='OCT')
axs.text(x=102.5, y=-1.65, s='NOV')
axs.text(x=112.5, y=-1.65, s='DEC')
#%% OA SA:  SUPERATIONS COMPOSITION

superations_days= dfi.copy(deep=True)
del superations_days['nr']
superations_days['nr']=superations_days[nr].sum()

oa_d=pd.DataFrame()
oa['date']=oa['datetime'].dt.date
oa_d_mean =oa.groupby('date').mean()

oa_d['datetime']=pd.date_range("01/01/2014", '31/12/2023')
oa_d.index=oa_d['datetime']
oa_d['date']=oa_d['datetime'].dt.date
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
oa_d['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in oa_d['date']]
oa_ds=oa_d.groupby('Season').count()

oa_daily = pd.merge(left = oa_d, right= oa_d_mean, left_index=True, right_index=True, how='outer')
oa_daily['datetime']=oa_daily.index
oa_daily['date']=oa_daily['datetime'].dt.date

mask_sup = oa_daily['date'].isin(superations_days['date'])
oa_sup=oa_daily.loc[mask_sup]

month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
oa_sup['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in oa_sup['date']]

oa_sup_my = oa_sup.groupby(['Season']).mean()
oa_sup_my_count=oa_sup.groupby(['Season']).count()
oa_sup_my['oa'] = oa_sup_my[factors].sum(axis=1)
oa_sup_my2 = pd.DataFrame(data = {'HOA':oa_sup_my['HOA']/oa_sup_my['oa'], 'COA':oa_sup_my['COA']/oa_sup_my['oa'], 
                                  'Amine-OA':oa_sup_my['Amine-OA']/oa_sup_my['oa'], 'BBOA':oa_sup_my['BBOA']/oa_sup_my['oa'],
                                  'LO-OOA':oa_sup_my['LO-OOA']/oa_sup_my['oa'],'MO-OOA':oa_sup_my['MO-OOA']/oa_sup_my['oa']})
oa_sup_my2 = 100.0*oa_sup_my2
fig, axs=plt.subplots(figsize=(4,4))
oa_sup_my2.plot(kind='bar', stacked=True, ax=axs, color=c_oa, zorder=3)
axs.legend(loc=(1.01, 0.55))
axs.set_ylabel('Average composition of OA (%)', fontsize=12)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Season', fontsize=12)
axs.set_ylim(0,108)
for i in range(0,4):
    axs.text(x=i-0.18, y=102, s=str((oa_sup_my_count['HOA'].iloc[i]*100/oa_ds['date'].iloc[i]).round(1))+'%')
axs.set_title('Superations of WHO daily limit\n', fontsize=13)

#%%
'''
++++++++++++++++++++++++  IONSSSSS  ++++++++++++++++++++++++++

'''
#%%% IONS: Inspection of the PMF files!
s=pd.read_csv('Specs.txt', sep='\t', header=None)
e=pd.read_csv('Errors.txt', sep='\t', header=None)
t=pd.read_csv('acsm_utc_time_mats.txt', sep='\t')
mz=pd.read_csv('amus.txt', sep='\t')
s.columns=mz['mz']
e.columns=mz['mz']
s.index=pd.to_datetime(t['Time (UTC)'], dayfirst=True)
e.index=pd.to_datetime(t['Time (UTC)'], dayfirst=True)
#Filtering specs
for i in s.columns:
    s[i][s[i]> 1] = np.nan
    e[i][e[i]> 1] = np.nan
    s[i][s[i]< -0.05] = np.nan
#%% IONS: Filtering ratios
ratios=pd.DataFrame()
ratios['datetime']=ratios.index
ratios['mz57mz55']=s[57]/s[55]
ratios['mz43mz44']=s[43]/s[44]
ratios['mz57mz44']=s[57]/s[44]
ratios['mz60mz44']=s[60]/s[44]
ratios['mz58mz44']=s[58]/s[44]
for i in ratios.columns:
    ratios[i][ratios[i]> 50] = np.nan
    ratios[i][ratios[i]< -50] = np.nan
#%% IONS: TS ratios
c_ions = ['darkviolet', 'forestgreen', 'gray', 'chocolate','cadetblue']
ratios['datetime'] =ratios.index
fig, axs = plt.subplots(figsize=(16,8), nrows = len(ions), sharex ='col', ncols =4, width_ratios=[3, 1,1, 1])
ratios['Year']=ratios['datetime'].dt.year
ratios['Month']=ratios['datetime'].dt.month
ratios['Hour']=ratios['datetime'].dt.hour
ratios_y, ratios_m, ratios_h=ratios.groupby('Year').mean(), ratios.groupby('Month').mean(), ratios.groupby('Hour').mean()
ratios_ye, ratios_me, ratios_he=ratios.groupby('Year').std(), ratios.groupby('Month').std(), ratios.groupby('Hour').std()

for i in range(0,len(ions)):
    ratios[ions[i]].plot( ax=axs[i,0], color=c_ions[i], legend=False)
    ratios_y[ions[i]].plot( ax=axs[i,1], color=c_ions[i], legend=False, marker='o', ms=4)
    # axs[i,1].errorbar(x=ratios_y.index, y=ratios_y[ions[i]], yerr=ratios_ye[ions[i]],color=c_ions[i])
    ratios_m[ions[i]].plot( ax=axs[i,2], color=c_ions[i], legend=False, marker='o', ms=4)
    # axs[i,2].errorbar(x=ratios_m.index, y=ratios_m[ions[i]], yerr=ratios_me[ions[i]],color=c_ions[i])
    ratios_h[ions[i]].plot( ax=axs[i,3], color=c_ions[i], legend=False, marker='o', ms=4)
    # axs[i,3].errorbar(x=ratios_h.index, y=ratios_h[ions[i]], yerr=ratios_he[ions[i]],color=c_ions[i])

    axs[i,0].set_ylabel(ratios.columns[i]+'\n(adim.)')
    axs[i,1].grid(axis='x'), axs[i,2].grid(axis='x'), axs[i,3].grid(axis='x')


axs[0,0].set_ylim(0,5)
axs[1,0].set_ylim(0,10)
axs[2,0].set_ylim(0,6)
axs[3,0].set_ylim(0,3)
axs[4,0].set_ylim(0,5)
axs[0,0].grid(axis='x'), axs[1,0].grid(axis='x'), axs[2,0].grid(axis='x'), axs[3,0].grid(axis='x'),axs[4,0].grid(axis='x')
fig.suptitle('Ions ratios')
#%% IONS:  SUPERATIONS COMPOSITION

superations_days= dfi.copy(deep=True)
del superations_days['nr']
superations_days['nr']=superations_days[nr].sum()

ratiosions=pd.DataFrame()
ratiosions=ratios.copy(deep=True)
ratiosions['date']=pd.to_datetime(ratios['datetime']).dt.date
ratiosions =ratiosions.groupby('date').mean()

r_d=pd.DataFrame()
r_d['datetime']=pd.date_range("01/01/2014", '31/12/2023')
r_d.index=r_d['datetime']
r_d['date']=r_d['datetime'].dt.date
month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
r_d['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in r_d['date']]

r_daily = pd.merge(left = r_d, right= ratiosions, left_index=True, right_index=True, how='outer')
r_daily['datetime']=r_d.index
r_daily['date']=oa_daily['datetime'].dt.date

mask_sup = r_daily['date'].isin(superations_days['date'])
r_sup=r_daily.loc[mask_sup]

month_to_season_dct = {1: 'DJF', 2: 'DJF',3: 'MAM', 4: 'MAM', 5: 'MAM',6: 'JJA', 7: 'JJA', 8: 'JJA',9: 'SON', 10: 'SON', 11: 'SON',12: 'DJF'}
r_sup['Season'] = [month_to_season_dct.get(t_stamp.month) for t_stamp in r_sup['date']]

r_sup_my = r_sup.groupby(['Season']).mean()
r_sup_my_count=r_sup.groupby(['Season']).count()

fig, axs=plt.subplots(figsize=(4,4))
r_sup_my[ions].plot(marker='o', ax=axs, color=c_ions, zorder=3 )
axs.legend(loc=(1.01, 0.55))
axs.set_ylabel('OA Ion ratios (adim.)', fontsize=12)
axs.grid(axis='y', zorder=0)
axs.set_xlabel('Season', fontsize=12)

axs.set_title('Superations of WHO daily limit\n', fontsize=13)

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



#%%
'''********************************* BOUNDARY LAYER ******************************************'''
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
#%%
'''+++++++++++++++++++++++++++++++++++ OA SA BL +++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

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
specs_time=pd.read_csv('acsm_utc_time_mats.txt', sep='\t')
specs_amus=pd.read_csv('amus.txt', sep='\t')
specs.columns=specs_amus['mz']
dt_specs=pd.to_datetime(specs_time['Time (UTC)'], dayfirst=True)
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
mask = [(ions['f44']>=0.0) & (ions['f44']<=1.0) & (ions['f43']>=0.0) & (ions['f43']<=1.0) & (ions['f73']>=0.0) & (ions['f60']>=0.0) & (ions['f55']>=0.0) & (ions['f57']>=0.0)]
ions_f = ions.loc[mask]
ions_f=ions.copy(deep=True)
ions_f.plot()
ions_f['Year']=ions['dt'].dt.year
ions_f['Month']=ions['dt'].dt.month
ions_f.to_csv('Ions.txt', sep='\t')
ions_ym=pd.DataFrame(ions_f[ions_labels].groupby([ions_f['Year'], ions_f['Month'],]).mean())
#%%  Plotting year/months ions
os.chdir(path_treatment)
ion=['f44', 'f43', 'f55', 'f57', 'f60', 'f73']  
fig, axs=plt.subplots(figsize=(10,10))
pl=ions_ym[ion].plot(ax=axs, x_compat=True, color='grey', subplots=True, sharex=True, lw=2, marker='o', grid=True)   
axs.set_xticks(range(0,))
axs.set_xlabel('Monthly means', fontsize=11)
plt.savefig('Ions_monthyear.png', bbox_inches='tight')
#%% Ions Year intercomp
ym_mean=ions_ym#.columns[0:6]
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
for k in range(0, len(yearly_ions.columns)):
    print(greys[k])
    pl=yearly_ions[years[k]].plot.line(fontsize=13,ax=axs,  grid=True, markersize=8, marker=markers[i], legend=True, lw=2)
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
#%%

''' BACKTRAJECTORIES '''


#%% IMproting and daily
bt = pd.read_csv('Retrotrajectories.txt', sep='\t')
bt['datetime']=pd.to_datetime(bt['Date'], dayfirst=True)
bt.index = bt['datetime']
#%% Grouping by basic backtrajectories only
bt_basic =[]
for i in range(0,len(bt)):
    # print(bt['Retro'].iloc[i][:3])
    if bt['Retro'].iloc[i][:3] == 'REG':
        bt_basic.append('REG')
    elif bt['Retro'].iloc[i][:3] == 'MED':
        bt_basic.append('MED')
    elif bt['Retro'].iloc[i][:3] == 'ANT':
        bt_basic.append('ANT')
    elif bt['Retro'].iloc[i][:3] == 'NAF':
        bt_basic.append('NAF')
    elif bt['Retro'].iloc[i][:3] == 'ANW':
        bt_basic.append('ANW')
    elif bt['Retro'].iloc[i][:3] == 'REG':
        bt_basic.append('REG')
    elif bt['Retro'].iloc[i][:3] == 'MED':
        bt_basic.append('MED')
    elif bt['Retro'].iloc[i][:3] == 'AMT':
        bt_basic.append('ANT')
    elif (bt['Retro'].iloc[i][:3] == 'EU') or (bt['Retro'].iloc[i][:3] == 'EU ') or (bt['Retro'].iloc[i][:3] == 'EU-') or (bt['Retro'].iloc[i][:3] == 'EU?'):
        bt_basic.append('EU')
    elif (bt['Retro'].iloc[i][:3] == 'AW') or (bt['Retro'].iloc[i][:3] == 'AW ') or (bt['Retro'].iloc[i][:3] == 'AW-') or (bt['Retro'].iloc[i][:3] == 'AW/'):
        bt_basic.append('AW')
    elif (bt['Retro'].iloc[i][:3] == 'AN') or (bt['Retro'].iloc[i][:3] == 'AN-') or (bt['Retro'].iloc[i][:3] == 'AN ') or (bt['Retro'].iloc[i][:3] == 'AN/') or (bt['Retro'].iloc[i][:3] == 'AN?') or (bt['Retro'].iloc[i][:3] == 'WAN'):
        bt_basic.append('AN')
    else:
        bt_basic.append(bt['Retro'].iloc[i][:3])

bt['Basic_retro']=bt_basic
#%% Daily nr and plotting
df['date']=df['datetime'].dt.date
dfd = df.groupby(df['date']).mean()

nrd = pd.merge(left = dfd, right= bt, left_index=True, right_index=True, how='outer')
nrbt= nrd.groupby('Basic_retro').mean()
nrbt_c= nrd.groupby('Basic_retro').count()
nrbt_c_p = pd.DataFrame(data = 100*nrbt_c['Org'] / nrbt_c['Org'].sum())

fig, axs=plt.subplots(figsize=(10,4))
nrbt.plot(kind='bar', stacked=True, ax=axs, color=c_df)
axs2=axs.twinx()
nrbt_c_p.plot(marker ='o', lw=0, ax=axs2, markersize=10, color='k', legend=False)
axs.set_ylabel('Daily concentrations $(\mu g·m^{-3})$', fontsize=12)
axs2.set_ylabel('Frequency (%)', fontsize=12)
x1, y1 = [-5,20], [15,15]
axs.plot(x1, y1, ls='--', color='grey')

axs.text(x=-0.35, y=15.5, s='WHO PM$_{2.5}$ daily limit', c='dimgrey', fontsize=12)
axs.text(x=-0.35, y=17.5, s='Mean', c='k', fontsize=12)
axs.set_xlabel('Backtrajectories scenarios', fontsize=12)
axs.set_ylim(0, 17)
axs.legend(loc=(1.08, 0.6))



#%%
limit_who_25_daily=15.5
df['date']=df['datetime'].dt.date
dfd = df.groupby(df['date']).mean()

nrd = pd.merge(left = dfd, right= bt, left_index=True, right_index=True, how='outer')
nrd['nr']=nrd['Org']+nrd['SO4']+nrd['NO3']+nrd['NH4']+nrd['Chl']

mask =nrd['nr']>=limit_who_25_daily
nrd_high=nrd.loc[mask]
nrd_high_c=nrd_high.groupby('Basic_retro').count()
nrd_high_m=nrd_high.groupby('Basic_retro').mean()

fig, axs=plt.subplots(figsize=(10,4))
nrd_high_m[nr].plot(kind='bar', stacked=True, ax=axs, color=c_df)
axs2=axs.twinx()
nrd_high_c['Org'].plot(marker ='o', lw=0, ax=axs2, markersize=10, color='k', legend=False)
axs.set_ylabel('Daily concentrations $(\mu g·m^{-3})$', fontsize=12)
axs2.set_ylabel('Frequency (%)', fontsize=12)
x1, y1 = [-5,20], [15,15]
axs.plot(x1, y1, ls='--', color='grey')

axs.text(x=-0.35, y=15.5, s='WHO PM$_{2.5}$ daily limit', c='grey', fontsize=12)
axs.text(x=-0.35, y=26, s='Days obove WHO PM$_{2.5}$ daily limit', c='k', fontsize=12)

axs.set_xlabel('Backtrajectories scenarios', fontsize=12)
# axs.set_ylim(0, 17)
axs.legend(loc=(1.08, 0.6))    
    
    
    
    


