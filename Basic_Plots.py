# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 15:43:32 2023

@author: Marta Via
"""

import pandas as pd
import numpy as np
import glob
import os as os
import datetime as dt
import matplotlib.pyplot as plt

#Si una funció no té self, no es pot cridar des de fora.

class Profiles:
    
    def __init__(self,x):
        self.x = 0
        pass
    
    def Hello(self):
        print("Hello Profile plotter")
    
    def Profiles_OA(self, df, relf, name, nf):
        fig,axes=plt.subplots(nrows=nf,ncols=2, figsize=(38,20), sharex='col',gridspec_kw=dict(width_ratios=[3, 1]))
        for i in range(0,nf):
            axes[i,0].bar(mz['lab'][:73], df.iloc[:,i][:73], color='grey')
            ax2=axes[i,0].twinx()
            axes[i,0].tick_params(labelrotation=0)
            ax2.plot(mz.lab[:73], relf.iloc[:,i][:73], marker='o', linewidth=False,color='black')
            axes[i,0].grid(axis='x')
            axes[i,1].bar(mz['lab'][73:], df.iloc[:,i][73:], color='grey')
            axes[i,1].set_yscale('log')
            ax4=axes[i,1].twinx()
            ax4.plot(mz.lab[73:], relf.iloc[:,i][73:], marker='o', linewidth=False,color='black')
            axes[i,1].grid(axis='x')
            axes[i,0].tick_params('x', labelrotation=90, labelsize=16)
            axes[i,0].tick_params('y', labelrotation=0, labelsize=16)
            axes[i,1].tick_params(labelrotation=0, labelsize=16)
            axes[i,0].set_ylabel(df.columns[i], fontsize=18)
            ax4.tick_params(labelrotation=0, labelsize=16)
            ax2.tick_params(labelrotation=0, labelsize=16)
            axes[i,0].set_xticklabels(mz.lab[:73], fontsize=17)
            axes[i,1].set_xticklabels(mz.lab[73:], fontsize=17)
            ax2.set_ylim(0,max(relf.iloc[:,i][:73]))
   

    def Profiles_MT(self, df, relf, name, nf):
            fig,axes=plt.subplots(nrows=nf,ncols=2, figsize=(22,15), sharex='col',gridspec_kw={'width_ratios': [1.5, 1]})
            for i in range(0,nf):
                axes[i,0].bar(mz['lab'][:73], df.iloc[:,i][:73], color='grey')
                OA_perc=(100.0*df.iloc[:73, i].sum()/df.iloc[:,i].sum()).round()
                ax2=axes[i,0].twinx()
                axes[i,0].tick_params(labelrotation=90)
                ax2.plot(mz.lab[:73], relf.iloc[:,i][:73], marker='o', linewidth=False,color='black')
                axes[i,0].grid(axis='x')
                axes[i,1].bar(mz['lab'][73:], df.iloc[:,i][73:], color='grey')
                axes[i,1].set_yscale('log')
                ax4=axes[i,1].twinx()
                ax4.plot(mz.lab[73:], relf.iloc[:,i][73:], marker='o', linewidth=False,color='black')
                # ax4.set_yscale('log')
                # axes[i,0].set_ylabel(labels_2[i], fontsize=17)
                axes[i,1].grid(axis='x')
                axes[i,0].set_ylabel(labels_2[i], fontsize=18)
                axes[i,0].tick_params(labelrotation=90, labelsize=16)
                axes[i,1].tick_params(labelrotation=90, labelsize=16)
                axes[i,0].tick_params('x', labelrotation=90, labelsize=16)
                axes[i,0].tick_params('y', labelrotation=0, labelsize=16)
                ax4.tick_params(labelrotation=0, labelsize=16)
                ax2.tick_params(labelrotation=0, labelsize=16)
                axes[i,0].set_xticklabels(mz.lab[:73], fontsize=17)
                axes[i,1].set_xticklabels(mz.lab[73:], fontsize=17)
                axes[i,0].text(x=-3,y=max(df.iloc[:73,i])*0.8,s='OA = '+str(OA_perc)+'%', fontsize=16 )
                # axes[i,0].tick_params(axis='both', which='major', labelsize=10)
                fig.tight_layout(pad=2.0)
                # axes[i,0].set_ylabel(labels_pie[i] + '\n')
                # fig.suptitle(name)  
    
    def Profiles_MT_std(self, df,df_std, relf, name, nf):
            fig,axes=plt.subplots(nrows=nf,ncols=2, figsize=(22,15), sharex='col',gridspec_kw={'width_ratios': [1.5, 1]})
            for i in range(0,nf):
                axes[i,0].bar(mz['lab'][:73], df.iloc[:,i][:73], color='grey')
                axes[i,0].errorbar(mz['lab'][:73], df.iloc[:,i][:73], yerr=df_std.iloc[:,i][:73],  fmt='.', color='k')#, uplims=True, lolims=True)
                OA_perc=(100.0*df.iloc[:73, i].sum()/df.iloc[:,i].sum()).round(0)
                ax2=axes[i,0].twinx()
                axes[i,0].tick_params(labelrotation=90)
                ax2.plot(mz.lab[:73], relf.iloc[:,i][:73], marker='o', linewidth=False,color='black')
                axes[i,0].grid(axis='x')
                axes[i,1].bar(mz['lab'][73:], df.iloc[:,i][73:], color='grey')
                axes[i,1].errorbar(mz['lab'][73:], df.iloc[:,i][73:], yerr=df_std.iloc[:,i][73:],  fmt='.', color='k')#, uplims=True, lolims=True)            
                axes[i,1].set_yscale('log')
                ax4=axes[i,1].twinx()
                ax4.plot(mz.lab[73:], relf.iloc[:,i][73:], marker='o', linewidth=False,color='black')
                # ax4.set_yscale('log')
                # axes[i,0].set_ylabel(labels_2[i], fontsize=17)
                axes[i,1].grid(axis='x')
                axes[i,0].set_ylabel(labels_2[i], fontsize=18)
                axes[i,0].tick_params(labelrotation=90, labelsize=16)
                axes[i,1].tick_params(labelrotation=90, labelsize=16)
                axes[i,0].tick_params('x', labelrotation=90, labelsize=16)
                axes[i,0].tick_params('y', labelrotation=0, labelsize=16)
                ax4.tick_params(labelrotation=0, labelsize=16)
                ax2.tick_params(labelrotation=0, labelsize=16)
                axes[i,0].set_xticklabels(mz.lab[:73], fontsize=17)
                axes[i,1].set_xticklabels(mz.lab[73:], fontsize=17)
                axes[i,0].text(x=-3,y=max(df.iloc[:73,i])*0.9,s='OA = '+str(OA_perc)+'%', fontsize=16 )
                # axes[i,0].tick_params(axis='both', which='major', labelsize=10)
                fig.tight_layout(pad=2.0)
                # axes[i,0].set_ylabel(labels_pie[i] + '\n')
                # fig.suptitle(name)   

    def Profiles_Filters(self df, relf, name, nf):
        fig,axes=plt.subplots(nrows=nf, figsize=(8,8), sharex='col')
        for i in range(0,nf):
            axes[i].bar(mz['lab'], df.iloc[:,i], color='gray')
            ax2=axes[i].twinx()
            axes[i].tick_params(labelrotation=90)
            ax2.plot(mz['lab'], relf.iloc[:,i], marker='o', linewidth=False,color='black')
            axes[i].grid()
            axes[i].set_ylabel(df.columns[i]+'\n(μg·m$^{-3}$)')
            axes[i].set_yscale('log')
        # fig.suptitle(name, fontsize=16)
        # plt.text(x=0, y=10, s='Concentration (μg·m$^{-3}$)', fontsize=14)# \t\t\t\t
        axes[num_f-1].set_xlabel('Species', fontsize=14) 

class PolarPlots()

    def __init__(self,x):
        self.x = 0
        pass
    
    def Hello(self):
        print("Hello Pollution rose plotter")

    def Pollution_rose(self, df, factor) :
        #df is the meteo file, which must include the wd, ws and the factors at the same timestamps.
        theta = np.linspace(0,360,17)
        r=np.linspace(0,10,11)
        li_theta=[]
        for i in range(0,len(theta)-1):
            mask=(df['wd']>= theta[i]) &  (df['wd']<theta[i+1])
            df1=df[mask]
            li_r=[]
            for j in range(0,len(r)-1):
                mask2=(df1['ws']>r[j])&(df1['ws']<r[j+1])
                df2=df1[mask2]
                li_r.append(df2[factor].mean())
            li_theta.append(li_r)
        theta_rad=(theta*np.pi/180.0).round(5)
        poll=pd.DataFrame(li_theta)#, columns=r[:-1],index=theta_rad[:-1])   
        fig, ax = plt.subplots(subplot_kw={"projection":"polar"}, figsize=(4,4))
        R, Theta=np.meshgrid(r,theta*np.pi/180.0)
        plot=ax.pcolormesh(Theta, R, poll, cmap='Greys', vmin=-0.19,vmax=(meteo2[factor].quantile(0.80)))
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location("N")  # theta=0 at the top
        ax.set_rticks(r)  # Less radial ticks
        ax.grid()
        ax.set_title('AS + Heavy oil combustion' +' \n$(\mu g·m^{-3})$', fontsize=15)
        cb=fig.colorbar(plot, ax=ax, orientation='horizontal')


class Residuals:
    def __init__(self,x):
        self.x = 0
        pass
    
    def Hello(self):
        print("Hello Pollution rose plotter")

    def Residuals_MT_histogram(nbins,res_1, res_2):
        path_py="C:/Users/maria/Documents/Marta Via/1. PhD/F. Scripts/Python Scripts"
        os.chdir(path_py)
        from Histograms import *
        hist = Histogram_treatment(zero)

        hist_intersection_r1=[]
        
        fig, axs=plt.subplots(figsize=(5,5), dpi=100)

        values_a, bins_a, patches_a = plt.hist(res_1, bins=4401, range=[-50,50])
        values_b, bins_b, patches_b = plt.hist(res_2, bins=21, range=[-50,50])
        
        fig2, axs2=plt.subplots(figsize=(4,4))
        values_a3=averages(list(bins_a[:-1]), list(bins_b), values_a)
        values_a3=pd.Series(values_a3).replace(np.nan,0)
        hist_int=histogram_intersection_new(values_a3.values, values_b, bins_b).round(1)
        hist_intersection_r1.append(hist_int)
        axs2.plot(bins_b[:-1],values_a3,  color='gray', label='HR')
        axs2.plot(bins_b[:-1],values_b, color='silver', label='LR')
        axs2.set_ylabel('Frequency (adim.)', fontsize=13)
        axs2.set_xlabel('Scaled residuals', fontsize=13)
        axs2.text(x=-50,y=14, s='F$^{*}$= '+str(hist_int.round(1))+'%', fontsize=12)
        axs2.legend(['HR', 'LR'], loc='upper left')    



#÷To investigate:
    #ax14.xaxis.set_major_formatter(mdates.DateFormatter('%m/%y'))
    #ax11.xaxis.set_major_locator(mdates.MonthLocator())
