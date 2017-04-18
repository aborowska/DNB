# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:44:13 2014

@author: istvan
"""

import csv
import math
import numpy as np
import re
import matplotlib 
from matplotlib import dates
import matplotlib.pyplot as plt
import datetime
#from sklearn.neighbors import KernelDensity
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from scipy.stats import skew



plt.rc('font', family='serif')

def acovf(x, unbiased=False, demean=True):
    '''
    Autocovariance for 1D

    Parameters
    ----------
    x : array
        Time series data. Must be 1d.
    unbiased : bool
        If True, then denominators is n-k, otherwise n
    demean : bool
        If True, then subtract the mean x from each element of x
   
    Returns
    -------
    acovf : array
        autocovariance function
    '''
    x = np.squeeze(np.asarray(x))
    if x.ndim > 1:
        raise ValueError("x must be 1d. Got %d dims." % x.ndim)
    n = len(x)

    if demean:
        xo = x - x.mean()
    else:
        xo = x
    if unbiased:
        xi = np.arange(1, n + 1)
        d = np.hstack((xi, xi[:-1][::-1]))
    else:
        d = n * np.ones(2 * n - 1)
   
    return (np.correlate(xo, xo, 'full') / d)[n - 1:]


def acf(x, nlags=40, unbiased=True):
  
    avf = acovf(x, unbiased=unbiased, demean=True)
    acf = avf[:nlags + 1] / avf[0]
 
    return acf

def CalculateIneff(Data):
    iN=len(Data)    
    iLag=iN-1   
    a=acf(Data,iLag,False)
    diff=a-np.divide( 2*np.ones(iLag+1), np.sqrt(np.array([iN]*(iLag+1))-np.array(range(0,iLag+1))))
    
    tempsum=0
    for i in range(0,iN):
        if(diff[i]>0):
            tempsum=tempsum+2*abs(a[i])
        else:
            return tempsum
        
      
    
def kde_sklearn(x, x_grid, bandwidth, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    
    # score_samples() returns the log-likelihood of the samples
   
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)
 
def FormatNum(Num, NumOfDig):
    Text=str(round(Num,NumOfDig))
   
    if(Text[len(Text)-2:len(Text)]!='.0'):    
        NumArray=Text.split( '.')
        Int=NumArray[0]
    else:
        NumArray=Text
        Int=Text[:len(Text)-2]
        
        
    Out=''
    for i in range(0,len(Int.lstrip('-'))/3):
        if i==0:
            Out=Int[len(Int)-(i+1)*3:len(Int)-i*3]
        else:    
            Out=Int[len(Int)-(i+1)*3:len(Int)-i*3]+' '+Out
        
    if(len(Int.lstrip('-'))/3==0):
        Out=Int[:len(Int)-(len(Int.lstrip('-'))/3)*3]
    else:
        Out=Int[:len(Int)-(len(Int.lstrip('-'))/3)*3]+' '+Out
 
    if len(NumArray)==2:    
        Out=Out+'.'+NumArray[1]
    
    return Out

def ResultTable(Date, RowNames,ColNames,aMean, aLow, aHigh,aInef,aESS):
    iNumOfCols=2*len(ColNames)+1
    iNumOfRows=len(RowNames)+2
    print iNumOfCols
    print iNumOfRows
    iNumOfInt=3    
    sTable='\\begin{tabular}{'
    for i in range(0, iNumOfCols):
        sTable=sTable+'c'
    sTable=sTable+'} \\toprule \n ' 
    for i in range(0,iNumOfRows):
        for j in range(0,iNumOfCols):        
            if i==0:
                if( (j>0) & (j<=len(ColNames))):
                    sTable+=' & \multicolumn{2}{c}{'+ColNames[j-1]+'}' 
            elif i==1:
                if((j>0) & (j%2==1) ):
                    sTable+='&  \\multicolumn{1}{c}{Skellam}'  
                if((j>0) & (j%2==0) ):
                    sTable+='&  \\multicolumn{1}{c}{$\\Delta$NB}'
            else:
                if j==0:
                    sTable+=RowNames[i-2]
                else:
                    if(aMean[i-2,j-1]!=0):
                        sTable+=' & '+FormatNum(aMean[i-2,j-1],iNumOfInt)
                    else:
                        sTable+=' & '
         
            
                    
        if i==iNumOfRows-1:
            sTable+=' \\\\ \n'
            for j in range(0,iNumOfCols):       
                if j>0:
                    if(aMean[i-2,j-1]!=0):
                        sTable+=' & \\begin{scriptsize} ['+FormatNum(aLow[i-2,j-1],iNumOfInt)+','+FormatNum(aHigh[i-2,j-1],iNumOfInt) +'] \\end{scriptsize} '   
                    else:
                        sTable+=' & '
            sTable=sTable+' \\\\  \\midrule \n'  
            for j in range(0,iNumOfCols):        
                if j==0:
                    sTable+='Inefficiency'
                else:        
                    sTable+=' & '+FormatNum(aInef[j-1],2)
            sTable+=' \\\\ \n'
            for j in range(0,iNumOfCols):        
                if j==0:
                    sTable+='ESS'
                else:        
                    sTable+=' & '+FormatNum(aESS[j-1],2)
            sTable+=' \\\\  \\bottomrule \n'
        elif i==1:
            sTable=sTable+' \\\\ \\midrule \n '
        elif i==0:
            sTable=sTable+' \\\\ \n'
            for j in range(0,iNumOfCols):
                if (j%2==1):
                   sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+2)+'} ' 
        else:            
            sTable=sTable+' \\\\ \n'
            for j in range(0,iNumOfCols):       
                if j>0:
                    if(aMean[i-2,j-1]!=0):
                        sTable+=' & \\begin{scriptsize} ['+FormatNum(aLow[i-2,j-1],iNumOfInt)+','+FormatNum(aHigh[i-2,j-1],iNumOfInt) +'] \\end{scriptsize} '   
                    else:
                        sTable+=' & '       
            sTable=sTable+' \\\\ \n'     
    
    
    
    sTable=sTable+'\n\\end{tabular}'

    f = open("../Paper/Tables/ResultTable"+Date +".tex", "w")
    f.write(sTable)
    f.close()       
    
def ResultTable2(Date, RowNames,ColNames,aMean, aLow, aHigh,aInef,aESS,TickPerRow):
    print ColNames
    iNumOfTables=int(np.ceil(float(len(ColNames))/TickPerRow))
    iNumOfCols=4*min(len(ColNames),TickPerRow)+1
    iNumOfRows=len(RowNames)+2
    print iNumOfTables
    print iNumOfCols
    print iNumOfRows
    iNumOfInt=3    
    sTable='\\begin{tabular}{'
    for i in range(0, iNumOfCols):
        sTable=sTable+'c'
    sTable=sTable+'} \\toprule \n '
    for k in range(0,iNumOfTables):    
        for i in range(0,iNumOfRows):
            for j in range(0,iNumOfCols):        
                if i==0:
                    if( (j>0) & (j<=TickPerRow) ):
                        print j                        
                        sTable+='& \\multicolumn{4}{c}{'+ColNames[k*TickPerRow +j-1]+'}'
                elif i==1:
                    if((j>0) & ((j)%4==1) ):
                        sTable+='&  \\multicolumn{1}{c}{Skellam}'  
                    if((j>0) & ((j)%4==2) ):
                        sTable+='&  \\multicolumn{1}{c}{Ord Norm}'  
                    if((j>0) & ((j)%4==0) ):
                        sTable+='&  \\multicolumn{1}{c}{Ord t}'      
                    if((j>0) & ((j)%4==3) ):
                        sTable+='&  \\multicolumn{1}{c}{$\\Delta$NB}'
                else:
                    if j==0:
                        sTable+=RowNames[i-2]
                    else:
                        if(aMean[i-2,(k*(iNumOfCols-1)+j)-1]!=0):
                            sTable+=' & '+FormatNum(aMean[i-2,(k*(iNumOfCols-1)+j)-1],iNumOfInt)
                        else:
                            sTable+=' & '
             
                
                        
            if i==iNumOfRows-1:
                sTable+=' \\\\ \n'
                for j in range(0,iNumOfCols):       
                    if j>0:
                        if(aMean[i-2,(k*(iNumOfCols-1)+j)-1]!=0):
                            sTable+=' & \\begin{scriptsize} ['+FormatNum(aLow[i-2,(k*(iNumOfCols-1)+j)-1],iNumOfInt)+','+FormatNum(aHigh[i-2,(k*(iNumOfCols-1)+j)-1],iNumOfInt) +'] \\end{scriptsize} '   
                        else:
                            sTable+=' & '
                sTable=sTable+' \\\\  \\midrule \n'  
                for j in range(0,iNumOfCols):        
                    if j==0:
                        sTable+='Inefficiency'
                    else:        
                        sTable+=' & '+FormatNum(aInef[(k*(iNumOfCols-1)+j)-1],2)
                sTable+=' \\\\ \n'
                for j in range(0,iNumOfCols):        
                    if j==0:
                        sTable+='ESS'
                    else:        
                        sTable+=' & '+FormatNum(aESS[(k*(iNumOfCols-1)+j)-1],2)
                if k==iNumOfTables-1:
                    sTable=sTable+' \\\\ \\bottomrule '
                else:
                    sTable=sTable+' \\\\ \\midrule \n '
            elif i==1:
                sTable=sTable+' \\\\ \\midrule \n '
            elif i==0:
                sTable=sTable+' \\\\ \n'
                for j in range(0,iNumOfCols):
                    if (j%4==1):
                       sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+4)+'} ' 
            else:            
                sTable=sTable+' \\\\ \n'
                for j in range(0,iNumOfCols):       
                    if j>0:
                        if(aMean[i-2,(k*(iNumOfCols-1)+j)-1]!=0):
                            sTable+=' & \\begin{scriptsize} ['+FormatNum(aLow[i-2,(k*(iNumOfCols-1)+j)-1],iNumOfInt)+','+FormatNum(aHigh[i-2,(k*(iNumOfCols-1)+j)-1],iNumOfInt) +'] \\end{scriptsize} '   
                        else:
                            sTable+=' & '
                sTable=sTable+' \\\\ \n'     
        
    
    
    sTable=sTable+'\n\\end{tabular}'

    f = open("../Paper/Tables/ResultTable"+Date +".tex", "w")
    f.write(sTable)
    f.close()           
    
    
def TracePlot(aData,sXAxisLabel,sFilename):
    iNumOfSample=len(aData)
    fig, ax1 = plt.subplots(subplot_kw={'axisbg':'#EEEEEE','axisbelow':True})
    ax1.grid(color='white', linestyle='-', linewidth=2)    
    lns1=ax1.plot(range(1,iNumOfSample+1), aData, linewidth=3, alpha=0.5, label=sXAxisLabel)        
    lns=lns1
    labs=[l.get_label() for l in lns]        
    ax1.legend(lns,labs, loc=0)
    ax1.set_xlabel('Iteration',fontsize=10)
    ax1.set_xlim(1,iNumOfSample)
    ax1.set_ylabel(sXAxisLabel, fontsize=10)
    for spine in ax1.spines.values():
        spine.set_color('#BBBBBB')
    plt.title('Parameter '+sXAxisLabel)  
    fig.savefig('trace_'+sFilename+'.png', bbox_inches=0) 
    fig.savefig('trace_'+sFilename+'.pdf', bbox_inches=0)    

def TracePlotAll(File, aData,Mean, iNumOfParam, iNumOfZeros, sXAxisLabel):
    iNumOfSample=len(aData)
    startPosition=np.ceil(iNumOfParam/2.0)*100+20+1
    fig = plt.figure()
    poz=0
    for i in range(0, iNumOfParam+iNumOfZeros):
        if Mean[i]!=0:
            print "plot "+str(i)  
            sub=startPosition+poz
            poz=poz+1
            ax1 = fig.add_subplot(sub,axisbg='#EEEEEE',axisbelow=True)
            ax1.grid(color='white', linestyle='-', linewidth=2)    
            lns1=ax1.plot(range(1,iNumOfSample+1), aData[:,i], linewidth=3, alpha=0.5, label=sXAxisLabel[i])        
            lns=lns1
    #        labs=[l.get_label() for l in lns]        
    #        ax1.legend(lns,labs, loc=0)
            ax1.set_xlim(1,iNumOfSample)
            ax1.set_xlabel(sXAxisLabel[i],fontsize=10)
            ax1.tick_params(axis='both', which='major',labelsize=7)
    #        ax1.set_ylabel(sXAxisLabel[i], fontsize=10)
            for spine in ax1.spines.values():
                spine.set_color('#BBBBBB')
        
    fig.suptitle('Trace plots of the parameters') 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92)     
    fig.savefig('../Paper/Pictures/'+File.rstrip('_EstimationResults.csv')+'trace_all.png', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+File.rstrip('_EstimationResults.csv')+'trace_all.pdf', bbox_inches=0)          
   
def DensityPlot(aData, Mean,Low, High,sXAxisLabel,sFilename):
        iNumOfSample=len(aData)
        x_grid = np.linspace(min(aData), max(aData), 1000)
        hpd_grid = np.linspace(Low, High, 1000)
#        grid = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0.1, 1.0, 30)}, cv=20) # 20-fold cross-validation
#        grid.fit(aData[:,i].tolist() )
#        kde = grid.best_estimator_
#        pdf = np.exp(kde.score_samples(x_grid[:, None]))
        bw=1.06*np.std(aData)*pow(iNumOfSample,-1.0/5) 
        pdf=kde_sklearn(aData, x_grid, bw)
        hpd=kde_sklearn(aData, hpd_grid, bw)
        pdfmean=kde_sklearn(aData, Mean, bw)
        #print pdfmean
        fig, ax1 = plt.subplots(subplot_kw={'axisbg':'#EEEEEE','axisbelow':True})
        ax1.grid(color='white', linestyle='-', linewidth=2)
        #SigmaG=np.sqrt(np.divide( float(6*(iNumOfSample-2)),((iNumOfSample+1)*(iNumOfSample+3))))
        #iNumOfBins=np.round(1+np.log2(iNumOfSample)+np.log2(1+np.abs(skew(aData[:,i]))/SigmaG))
        #iNumOfBins=np.round(np.sqrt(iNumOfSample))
        iNumOfBins=np.round(1+np.log2(iNumOfSample))
        #iNumOfBins=np.round(2*pow(iNumOfSample,1.0/3))
        ax1.hist(aData, iNumOfBins, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
        lns1=ax1.plot(x_grid, pdf, linewidth=1, alpha=0.5,color='b', label='kernel density estimate')
        lns3=ax1.axvline(x=aMean[i], ymin=0, ymax=pdfmean/ax1.get_ylim()[1], color='black')
        ymin=ax1.get_ylim()[0]
        ymax=ax1.get_ylim()[1]
        lns2=ax1.fill_between(hpd_grid, 0, hpd, hpd > 0, color='b', alpha=.25,label='hpd region')
        ax1.set_ylim([ymin,ymax]) 
        lns=lns1
        labs=[l.get_label() for l in lns]        
        ax1.legend(lns,labs, loc=0,prop={'size':10})
        plt.title('Posterior density of parameter '+sXAxisLabel)  
        ax1.set_xlabel(sXAxisLabel,fontsize=14)
        #ax.set_rasterized(True)
        ax1.tick_params(axis='both', which='major',labelsize=10)
        for spine in ax1.spines.values():
            spine.set_color('#BBBBBB')
        fig.savefig('dens_'+sFilename+'.png', bbox_inches=0)
        fig.savefig('dens_'+sFilename+'.pdf', bbox_inches=0)   
        
        
def DensityPlotAll(File, aData, iNumOfParam,iNumOfZeros, Mean,Low, High,sXAxisLabel):
    iNumOfSample=len(aData)
    startPosition=np.ceil(iNumOfParam/2.0)*100+20+1
    fig = plt.figure()

    poz=0
    for i in range(0, iNumOfParam+iNumOfZeros): 
        if Mean[i]!=0:
            print "plot "+str(i)          
            x_grid = np.linspace(min(aData[:,i]), max(aData[:,i]), 1000)
            hpd_grid = np.linspace(Low[i],High[i], 1000)
    #        grid = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0.1, 1.0, 30)}, cv=20) # 20-fold cross-validation
    #        grid.fit(aData[:,i].tolist() )
    #        kde = grid.best_estimator_
    #        pdf = np.exp(kde.score_samples(x_grid[:, None]))
            bw=1.06*np.std(aData[:,i])*pow(iNumOfSample,-1.0/5) 
            if bw>0:
                pdf=kde_sklearn(aData[:,i], x_grid, bw)
                hpd=kde_sklearn(aData[:,i], hpd_grid, bw)
                pdfmean=kde_sklearn(aData[:,i], np.array([Mean[i]]), bw)
                
                sub=startPosition+poz
                poz=poz+1
            
#                ax1 = fig.add_subplot(sub,axisbg='#EEEEEE',axisbelow=True)
#                ax1.grid(color='white', linestyle='-', linewidth=1)
                ax1 = fig.add_subplot(sub)
                #SigmaG=np.sqrt(np.divide( float(6*(iNumOfSample-2)),((iNumOfSample+1)*(iNumOfSample+3))))
                #iNumOfBins=np.round(1+np.log2(iNumOfSample)+np.log2(1+np.abs(skew(aData[:,i]))/SigmaG))
                #iNumOfBins=np.round(np.sqrt(iNumOfSample))
                iNumOfBins=np.round(1+np.log2(iNumOfSample))
                #iNumOfBins=np.round(2*pow(iNumOfSample,1.0/3))
                ax1.hist(aData[:,i], iNumOfBins, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
                lns1=ax1.plot(x_grid, pdf, linewidth=1, alpha=0.3,color='b', label='kernel density estimate')        
                lns3=ax1.axvline(x=Mean[i], ymin=0, ymax=pdfmean/ax1.get_ylim()[1], color='black')
                ymin=ax1.get_ylim()[0]
                ymax=ax1.get_ylim()[1]
                lns2=ax1.fill_between(hpd_grid, 0, hpd, hpd > 0, color='b', alpha=.25,label='hpd region')
                ax1.set_ylim([ymin,ymax])        
                lns=lns1
                labs=[l.get_label() for l in lns]        
                #ax1.legend(lns,labs, loc=2,prop={'size':8})
                 
                ax1.set_xlabel(sXAxisLabel[i],fontsize=10)
                ax1.tick_params(axis='both', which='major',labelsize=7)
#                for spine in ax1.spines.values():
#                    spine.set_color('#BBBBBB')
                
    fig.suptitle('Posterior densities of the parameters') 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92)     
    fig.savefig('../Paper/Pictures/'+File.rstrip('_EstimationResults.csv')+'dens_all.png', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+File.rstrip('_EstimationResults.csv')+'dens_all.pdf', bbox_inches=0)           
    


def ProcessLLOutPut(SkellamFiles,DNBFiles, OrdNFiles,OrdTFiles,OutSampleFiles):
    iNumOfFiles=len(SkellamFiles)
    startPosition=np.ceil(iNumOfFiles/2.0)*100+20+1
#    startPosition=111    
    fig = plt.figure()

    for i in range(0, iNumOfFiles):  
        print "plot "+str(i)

        Temp=[]        
        with open(OutSampleFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        Time=np.array(Temp, dtype=np.float64)  
        Return=Time[:,1]
        Time=Time[:,0]
        delta = [ datetime.timedelta(seconds=(s+3600*9.5) ) for s in Time] 
        
        year=int(OutSampleFiles[i][5:9])
        month=int(OutSampleFiles[i][9:11])
        if(len(OutSampleFiles[i].split('_')[1])==8):  
            day=int(OutSampleFiles[i][11:13])
        else:
            day=int(OutSampleFiles[i][11:12])
        dt= [datetime.datetime(year, month, day,0,0,0)+ k for k in delta ]
        t = np.array([ dates.date2num(l) for l in dt ] )

        Temp=[]        
        with open(SkellamFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        Skellam=np.array(Temp, dtype=np.float64)                
        Temp=[]        
        with open(DNBFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        DNB=np.array(Temp, dtype=np.float64) 
        Temp=[] 
        with open(OrdNFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        OrdN=np.array(Temp, dtype=np.float64)
        Temp=[]         
        with open(OrdTFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        OrdT=np.array(Temp, dtype=np.float64)
        
        BFSkellam=2*(np.cumsum(Skellam)-np.cumsum(OrdN))
        BFDNB=2*(np.cumsum(DNB)-np.cumsum(OrdN))
        BFOrdT=2*(np.cumsum(OrdT)-np.cumsum(OrdN))
        hfmt = dates.DateFormatter(' %H:%M \n %d/%m/%y ')
        sub=startPosition+i   
        ax1 = fig.add_subplot(int(sub),axisbg='#EEEEEE',axisbelow=True)
#        ax1 = fig.add_subplot(int(sub),axisbelow=True)
        ax1.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax1.xaxis.set_major_formatter(hfmt)
        lns1=ax1.plot(t, BFOrdT, linewidth=1,  label='Ord t',color='red', zorder=2) 
        lns2=ax1.plot(t, BFSkellam, linewidth=1,  label='Skellam',color='blue',zorder=2)     
        lns3=ax1.plot(t, BFDNB, linewidth=1,  label=' $\Delta$NB',color='black', zorder=2)     
        lns4=ax1.plot(t, 2*np.ones(np.size(t)), linewidth=1,color='r',ls='--',  label=SkellamFiles[i].split('_')[1],zorder=3)   
        lns5=ax1.plot(t,-2*np.ones(np.size(t)), linewidth=1,color='r',ls='--',  label=SkellamFiles[i].split('_')[1],zorder=3)   
        lns=lns1+lns2+lns3
        labs=[l.get_label() for l in lns]        
        legend=ax1.legend(lns,labs, loc=0,prop={'size':6})
        ax1.set_xlim(min(t), max(t))
        ax2=ax1.twinx()
        ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax2.xaxis.set_major_formatter(hfmt)
        ax2.xaxis.grid(color='black', linestyle=':', linewidth=1,which='major')
#        ax1.grid(color='white', linestyle='-', linewidth=1,which='both')
#        ax2.yaxis.grid(color='white', linestyle='-', linewidth=1)
     
        
        ln4=ax2.plot(t,Return, linewidth=1,color='g', alpha=0.7, label='tick returns',zorder=1)  
        plt.title(SkellamFiles[i].split('_')[1], fontsize=9)

        ax1.set_zorder(1)
        ax2.set_zorder(0)
        ax1.set_frame_on(False)
        ax2.set_frame_on(True)
        ax2.set_axisbelow(True)
#        ax2.set_axis_bgcolor('#EEEEEE')
      
        ax2.set_xlim(min(t), max(t))   
        ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax2.xaxis.set_major_formatter(hfmt)
        ax2.xaxis.grid(color='black', linestyle=':', linewidth=1,which='major')
#        ax1.set_axis_bgcolor('#EEEEEE')
        ax1.tick_params(axis='both', which='major',labelsize=7, )
        ax1.tick_params(axis='y', which='major',labelsize=7, colors='blue')       
        ax1.spines['top'].set_color('#BBBBBB')
        ax1.spines['bottom'].set_color('#BBBBBB')
        ax2.tick_params(axis='both', which='major',labelsize=7)
        ax2.tick_params(axis='y', which='major',labelsize=7, colors='green')
        ax2.spines['top'].set_color('#BBBBBB')
        ax2.spines['bottom'].set_color('#BBBBBB')
        ax1.spines['left'].set_color('blue')
        ax1.spines['right'].set_color('green')
        ax2.spines['right'].set_color('green')
        ax2.spines['left'].set_color('blue')
        ax1.set_ylabel('2 $\log$ BF', fontsize=7,color='blue')
        ax2.set_ylabel('Tick returns', fontsize=7,color='green')
#        ax1.set_ylabel(sXAxisLabel[i], fontsize=10)
#        for spine in ax1.spines.values():
#            spine.set_color('#BBBBBB')
    if(month==4):
        smonth='April'
    else:
        smonth='October'       
    fig.suptitle('Out of sample predicitve likelihood comparison in '+smonth +' '+str(year)) 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92, wspace=0.35, right=0.93,left=0.07, hspace=0.35)       
    fig.savefig('../Paper/Pictures/'+SkellamFiles[0][0:4]+'_OutLL.png', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+SkellamFiles[0][0:4]+'_OutLL.pdf', bbox_inches=0)          
#    fig.savefig('../Poster/'+SkellamFiles[0][0:4]+'_OutLL.png', bbox_inches=0)
#    fig.savefig('../Poster/'+SkellamFiles[0][0:4]+'_OutLL.pdf', bbox_inches=0)   


def ProcessBICOutPut(SkellamFiles,DNBFiles,NormalFiles,TFiles,EstimatesSk, EstimatesDNB,EstimatesNormal, EstimatesT,InSampleFile):
    
    iNumOfFiles=len(SkellamFiles)
    startPosition=np.ceil(iNumOfFiles/2.0)*100+20+1
#   startPosition=111  
    fig = plt.figure()

    for i in range(0, iNumOfFiles):  
        print "plot "+str(i)
        
        Temp=[] 
        with open(EstimatesSk[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        aData=np.array(Temp, dtype=np.float64)   
        aMean=np.mean(aData,axis=0)
        iNumOfParamSkellam=np.count_nonzero(aMean)
        print iNumOfParamSkellam
        Temp=[] 
        with open(EstimatesDNB[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        aData=np.array(Temp, dtype=np.float64)   
        aMean=np.mean(aData,axis=0)
        iNumOfParamDNB=np.count_nonzero(aMean)
        print iNumOfParamDNB
        Temp=[] 
        with open(EstimatesNormal[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        aData=np.array(Temp, dtype=np.float64)   
        aMean=np.mean(aData,axis=0)
        iNumOfParamNormal=np.count_nonzero(aMean)
        print iNumOfParamNormal
        Temp=[] 
        with open(EstimatesT[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        aData=np.array(Temp, dtype=np.float64)   
        aMean=np.mean(aData,axis=0)
        iNumOfParamT=np.count_nonzero(aMean)
        print iNumOfParamT
        Temp=[]      
        with open(InSampleFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        Time=np.array(Temp, dtype=np.float64)  
        
        Return=Time[:,1]
        Time=Time[:,0]
        
        delta = [ datetime.timedelta(seconds=(s+3600*9.5) ) for s in Time] 
        
        year=int(InSampleFiles[i][5:9])
        month=int(InSampleFiles[i][9:11])
        if(len(InSampleFiles[i].split('_')[1])==8):  
            day=int(InSampleFiles[i][11:13])
        else:
            day=int(InSampleFiles[i][11:12])
            
            
        dt=[]
        index=[]
        index.append(0)
        days=[]
        days.append(datetime.date(year, month, day).strftime('%d/%m'))
        PrevTime=0
        jump=0
        for k in range(0,np.size(Time)):
            if Time[k]<PrevTime:
                print day
                if jump==0:
                    day=day+3
                else:
                    day=day+1
                jump=jump+1
                days.append(datetime.date(year, month, day).strftime('%d/%m'))
                index.append(k)
                print day
                print 'XXXXX'
                
         
            dt.append(datetime.datetime(year, month, day,0,0,0)+delta[k])
            PrevTime=Time[k]
        
        t = np.array([ dates.date2num(l) for l in dt ] )
        iNumOfObs=np.count_nonzero(t)
        print iNumOfObs

        Temp=[]        
        with open(SkellamFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        Skellam=np.array(Temp, dtype=np.float64)                
        Temp=[]        
        with open(DNBFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        DNB=np.array(Temp, dtype=np.float64)          
        Temp=[]        
        with open(NormalFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                if(line[size-1]=='-inf'):
                    line[size-1]=-100
                Temp.append(line)
        Normal=np.array(Temp, dtype=np.float64)     
        Temp=[]        
        with open(TFiles[i], 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        T=np.array(Temp, dtype=np.float64)     
        
        BFOrdT=-((-2*np.cumsum(T)+iNumOfParamT*np.log(np.arange(1,iNumOfObs+1)))-(-2*np.cumsum(Normal)+iNumOfParamNormal*np.log(np.arange(1,iNumOfObs+1))))
        BFSkellam=-((-2*np.cumsum(Skellam)+iNumOfParamSkellam*np.log(np.arange(1,iNumOfObs+1)))-(-2*np.cumsum(Normal)+iNumOfParamNormal*np.log(np.arange(1,iNumOfObs+1))))
        BFDNB=-((-2*np.cumsum(DNB)+iNumOfParamDNB*np.log(np.arange(1,iNumOfObs+1)))-(-2*np.cumsum(Normal)+iNumOfParamNormal*np.log(np.arange(1,iNumOfObs+1))))
   
        hfmt = dates.DateFormatter(' %d/%m/%y ')
        sub=startPosition+i   
        ax1 = fig.add_subplot(int(sub),axisbg='#EEEEEE',axisbelow=True)
#        ax1 = plt.subplot(sub,axisbg='#EEEEEE',axisbelow=True)
#        ax1 = fig.add_subplot(sub)
       
       
        ax2 = ax1.twinx() 
#        ax1.grid(color='white', linestyle='-', linewidth=2)    
#        ax1.xaxis.set_major_locator(dates.HourLocator(interval=24))
       
        ax1.set_xticks(index)
        ax1.set_xticklabels(days)
        
#        ax1.set_xticklabels(days,rotation=45, horizontalalignment='right')
       
       
        lns1=ax1.plot(np.arange(np.size(Time)), BFOrdT, linewidth=1,  label='Ord t',color='red', zorder=2) 
        lns2=ax1.plot(np.arange(np.size(Time)), BFSkellam, linewidth=1,  label='Skellam',color='blue',zorder=2)     
        lns3=ax1.plot(np.arange(np.size(Time)), BFDNB, linewidth=1,  label=' $\Delta$NB',color='black', zorder=2)     
        lns4=ax1.plot(np.arange(np.size(Time)), 2*np.ones(np.size(t)), linewidth=1,color='r',ls='--',  label=SkellamFiles[i].split('_')[1])   
        lns5=ax1.plot(np.arange(np.size(Time)),-2*np.ones(np.size(t)), linewidth=1,color='r',ls='--',  label=SkellamFiles[i].split('_')[1])  
        
        ln6=ax2.plot(np.arange(np.size(Time)),Return, alpha=0.7, linewidth=1,color='g',  label='tick returns')  
        lns=lns1+lns2+lns3
        labs=[l.get_label() for l in lns]        
        legend=ax1.legend(lns,labs, loc=0,prop={'size':6})
        ax1.set_zorder(1)
        ax2.set_zorder(0)
        ax1.set_frame_on(False)
        ax2.set_frame_on(True)
        ax2.set_axisbelow(True)
#        labs=[l.get_label() for l in lns]        
#        ax1.legend(lns,labs, loc=0)
        plt.title(SkellamFiles[i].split('_')[1], fontsize=9)
#        ax1.set_xlim(min(t), max(t))
#        ax1.set_xlabel(sXAxisLabel[i],fontsize=10)
         
        ax1.tick_params(axis='both', which='major',labelsize=7, )
        ax1.tick_params(axis='y', which='major',labelsize=7, colors='blue')
        ax1.spines['left'].set_color('blue')
        ax1.spines['right'].set_color('green')
        ax1.spines['top'].set_color('#BBBBBB')
        ax1.spines['bottom'].set_color('#BBBBBB')
        ax2.tick_params(axis='both', which='major',labelsize=7)
        ax2.tick_params(axis='y', which='major',labelsize=7, colors='green')
        ax2.spines['right'].set_color('green')
        ax2.spines['top'].set_color('#BBBBBB')
        ax2.spines['bottom'].set_color('#BBBBBB')
        ax1.set_ylabel('2 $\log$ BF', fontsize=7,color='blue')
        ax2.set_ylabel('Tick returns', fontsize=7,color='green')
#        for spine in ax1.spines.values():
#            spine.set_color('#BBBBBB')
            
    if(month==4):
        smonth='April'
    else:
        smonth='October'
    fig.suptitle('In-sample BIC comparison in '+smonth +' '+str(year)) 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92, wspace=0.35, right=0.93,left=0.07, hspace=0.25)     

    fig.savefig('../Paper/Pictures/'+SkellamFiles[0][0:4]+'_InLL.png', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+SkellamFiles[0][0:4]+'_InLL.pdf', bbox_inches=0)   

#    fig.savefig('../Poster/'+SkellamFiles[0][0:4]+'_InLL.png', bbox_inches=0)
#    fig.savefig('../Poster/'+SkellamFiles[0][0:4]+'_InLL.pdf', bbox_inches=0)   

    
def ProcessMCMCOutPut(Files,BurnIn):
    iNumOfFiles=len(Files)
    iNumOfP=7    
    aMean=np.zeros((iNumOfP,iNumOfFiles)) 
    aHPDLow=np.zeros((iNumOfP,iNumOfFiles))
    aHPDHigh=np.zeros((iNumOfP,iNumOfFiles))     
    aIneff=np.zeros((iNumOfP,iNumOfFiles))  
    aESS=np.zeros((iNumOfP,iNumOfFiles))  
    numFiles=0
    aMaxInef=[]
    aMinESS=[]
    ColNames=[]
    for f in Files:
        print 'XXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXXXXXXXXX'
        print 'XXXXXXXXX '+f+' XXXXXXXXXXX'
        print 'XXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxXXXXXXXXXX'
        Temp=[]
        with open(f, 'r') as file:
            for record in file:
                line=record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Temp.append(line)
        aData=np.array(Temp, dtype=np.float64)        
       
        ''' Initialize variables'''    
        aData=aData[BurnIn:,:]    
        
        iNumOfParam=np.size(aData,1)
        iNumOfSample=np.size(aData,0)
        print 'Number of sample: '+str(iNumOfSample)
        if iNumOfSample <=0:
            print 'Burn in sample is longer than the actual number of draws!!! '
            
        iMaxLag=3000;
        dCoverage=0.95;  
        
        if(f.find('DNB')!=-1 or f.find('OrdT')!=-1):
#            ParamNamesLatex=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$\gamma$', r'$ \beta_{1} $', r'$ \beta_{2} $',r'$ \beta_{3} $',r'$ \beta_{4} $' ,'$\\nu$' ]
#            ParamNames=['mu', 'phi' ,'sigma2','gamma', 'beta1', 'beta2','beta3 ','beta4' , 'nu'  ]
            ParamNamesLatex=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$\gamma$', r'$ \beta_{1} $' ,r'$ \beta_{2} $' ,'$\\nu$' ]
            ParamNames=['mu', 'phi' ,'sigma2','gamma', 'beta1','beta2',  'nu'  ]
        else:    
#            ParamNamesLatex=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$\gamma$', r'$ \beta_{1} $', r'$ \beta_{2} $',r'$ \beta_{3} $',r'$ \beta_{4} $'  ]
#            ParamNames=['mu', 'phi' ,'sigma2','gamma', 'beta1', 'beta2','beta3 ','beta4'  ]
            ParamNamesLatex=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$ \gamma $', r'$ \beta_1 $',r'$ \beta_2 $'  ]
            ParamNames=['mu', 'phi' ,'sigma2','gamma', 'beta1', 'beta2' ]
        
        
        print 'Number of parameters: '+str(iNumOfParam)
        iNumOfZero=0
        for i in range(0,iNumOfParam):   
            ''' Mean HPD region '''        
            aMean[i,numFiles]=np.mean(aData[:,i])        
            
            '''  HPD region '''                
            TempDataSorted=np.sort(aData[:,i],axis=0)
            iNumOfInterval=int(np.floor(iNumOfSample*dCoverage))
            for j in range(0,iNumOfSample-iNumOfInterval):
                if j==0:
                    Low=TempDataSorted[j]
                    High=TempDataSorted[j+iNumOfInterval]
                    Range=TempDataSorted[j+iNumOfInterval]-TempDataSorted[j]
                else:
                    TempRange=TempDataSorted[j+iNumOfInterval]-TempDataSorted[j]
                    if(TempRange<Range):
                        Range=TempRange
                        Low=TempDataSorted[j]
                        High=TempDataSorted[j+iNumOfInterval]
            
            print ParamNames[i]
#            if(Low<0 and High<0):
#                aHPDLow[i,numFiles]=High
#                aHPDHigh[i,numFiles]=Low
#            else:
#                aHPDLow[i,numFiles]=Low
#                aHPDHigh[i,numFiles]=High
             
            aHPDLow[i,numFiles]=Low
            aHPDHigh[i,numFiles]=High
                    
            ''' Inefficiency '''        
#            dTempIneff=1  
#            iLag=1
#            dTempAcf=acf(aData[:,i],iLag,False)
#            print 'acf: '+str(np.abs(dTempAcf[iLag]))
#            print 'threshold: '+str(2.0/np.sqrt(iNumOfSample-iLag))
#            while(np.abs(dTempAcf[iLag]) > (2/np.sqrt(iNumOfSample-iLag)) and iLag<=iMaxLag):
#                dTempIneff=dTempIneff+2*dTempAcf[iLag]
#                iLag=iLag+1
#                dTempAcf=acf(aData[:,i],iLag,False)
#                if iLag==iNumOfSample:
#                    break
#                print '+++++++++++++++++++'  
#                print i
#                print 'Lag: '+str(iLag)
#                print 'acf: '+str(np.abs(dTempAcf[iLag]))
#                print 'threshold: '+str(2.0/np.sqrt(iNumOfSample-iLag))
#                print 'dTempIneff: '+str(dTempIneff)
            if aMean[i,numFiles]!=0:  
                
                aIneff[i,numFiles]=CalculateIneff(aData[:,i])
                aESS[i,numFiles]=iNumOfSample/aIneff[i,numFiles]
            else:
                iNumOfZero=iNumOfZero+1
                aESS[i,numFiles]=iNumOfSample
#            ''' Trace plot '''
#            TracePlot(aData[:,i],ParamNamesLatex[i],ParamNames[i])
#            
#            
#            ''' Density plot '''
#            DensityPlot(aData[:,i], aMean[i],aHPDLow[i,0],aHPDHigh[i,0],ParamNamesLatex[i],ParamNames[i])  
        
       
        print 'iNumOfZero'
        print iNumOfZero
        DensityPlotAll(f, aData, iNumOfParam-iNumOfZero,iNumOfZero, aMean[:,numFiles], aHPDLow[:,numFiles],   aHPDHigh[:,numFiles],ParamNamesLatex)   
        TracePlotAll(f, aData,  aMean[:,numFiles],iNumOfParam-iNumOfZero,iNumOfZero, ParamNamesLatex)
        
          
        aMaxInef.append(max(aIneff[:iNumOfParam,numFiles]))
        aMinESS.append(min(aESS[:iNumOfParam,numFiles]))
        
        ColNames.append(f.split('_')[1])
        print ColNames
        numFiles+=1   
        
    
    ColNamesU=[ ColNames[x] for x in range(0,len(ColNames)) if (x % 4) == 0] 
    with open('results_mean_'+Files[0].split('_')[1]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(aMean)  
    with open('results_low_'+Files[0].split('_')[1]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(aHPDLow)  
    with open('results_high_'+Files[0].split('_')[1]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(aHPDHigh)  
    with open('results_ineff_'+Files[0].split('_')[1]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(aIneff)  
    with open('results_ess_'+Files[0].split('_')[1]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(aESS)  
   
   
    RowNames=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$\gamma$', r'$ \beta_{1} $', r'$ \beta_{2} $','$\\nu$' ]
#    RowNames=['$ \mu $', '$ \phi $' ,'$ \sigma^{2} $','$\gamma$', r'$ \beta_{1} $','$\\nu$' ]   
#    ResultTable(Files[0].split('_')[0], RowNames,ColNamesU,aMean,aHPDLow, aHPDHigh,aMaxInef,aMinESS)     
    ResultTable2(Files[0].split('_')[0], RowNames,ColNamesU,aMean,aHPDLow, aHPDHigh,aMaxInef,aMinESS,2)     
        
        
def Sesonal(DataIn,XIn, SIn, IntIn):
    Temp=[] 
    with open(DataIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    Return=aData[:,1]
    Time=aData[:,0]
    Temp=[] 
    with open(XIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    X=aData[:]    
    
    Temp=[] 
    with open(SIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    S=aData[:]
    
    Temp=[] 
    with open(IntIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    LogInt=aData[:]
    
        
    delta = [ datetime.timedelta(seconds=(s+3600*9.5) ) for s in Time] 
        
    year=int(DataIn[5:9])
    month=int(DataIn[9:11])
    if(len(DataIn.split('_')[1])==8):  
        day=int(DataIn[11:13])
    else:
        day=int(DataIn[11:12])
            
            
    dt=[]
    index=[]
    index.append(0)
    days=[]
    days.append(datetime.date(year, month, day).strftime('%d/%m'))
    PrevTime=0
    jump=0
    for k in range(0,np.size(Time)):
        if Time[k]<PrevTime:
            if jump==0:
                day=day+3
            else:
                day=day+1
            jump=jump+1
            days.append(datetime.date(year, month, day).strftime('%d/%m'))
            index.append(k)
           
                
         
        dt.append(datetime.datetime(year, month, day,0,0,0)+delta[k])
        PrevTime=Time[k]
        
    t = np.array([ dates.date2num(l) for l in dt ] )   
    fig = plt.figure()
    ax1 = fig.add_subplot(411)
    lns1=ax1.plot(np.arange(np.size(Time)), Return, linewidth=1,  label='tick retruns', color='green')
    ax1.set_xticks(index)
    ax1.set_xticklabels(days)
    ax1.tick_params(axis='both', which='major',labelsize=7 )
    ax1.set_xlim(1, np.size(Time))
    ax1.grid()
    lns=lns1
    labs=[l.get_label() for l in lns]        
    ax1.legend(lns,labs, loc=0,prop={'size':8})
    ax2 = fig.add_subplot(412)
    lns2=ax2.plot(np.arange(np.size(Time)), X, linewidth=1,  label='$x_t$')
    ax2.set_xticks(index)
    ax2.set_xticklabels(days) 
    ax2.tick_params(axis='both', which='major',labelsize=7 )
    ax2.set_xlim(1, np.size(Time))
    ax2.grid()
    lns=lns2
    labs=[l.get_label() for l in lns]        
    ax2.legend(lns,labs, loc=0,prop={'size':8})
    ax3 = fig.add_subplot(413)
    lns3=ax3.plot(np.arange(np.size(Time)), S, linewidth=1,  label='$s_t$')
    ax3.set_xticks(index)
    ax3.set_xticklabels(days)    
    ax3.tick_params(axis='both', which='major',labelsize=7 )   
    ax3.set_xlim(1, np.size(Time))
    ax3.grid()
    lns=lns3
    labs=[l.get_label() for l in lns]        
    ax3.legend(lns,labs, loc=0,prop={'size':8})    
    ax4 = fig.add_subplot(414)
    lns4=ax4.plot(np.arange(np.size(Time)), LogInt, linewidth=1,  label='$\log \\lambda_t$')
    ax4.set_xticks(index)
    ax4.set_xticklabels(days) 
    ax4.tick_params(axis='both', which='major',labelsize=7 )
    ax4.set_xlim(1, np.size(Time))
    ax4.grid()
    lns=lns4
    labs=[l.get_label() for l in lns]        
    ax4.legend(lns,labs, loc=0,prop={'size':8})
    fig.suptitle('Volatility decompostion of '+ DataIn.split('_')[3].strip('.csv')+' tick returns from 23rd to 29th April 2010') 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92, wspace=0.35, right=0.93,left=0.07, hspace=0.2)   
    fig.savefig('../Paper/Pictures/'+DataIn.split('_')[3].strip('.csv')+'_season.pdf', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+DataIn.split('_')[3].strip('.csv')+'_season.png', bbox_inches=0)   


def Sesonal2(DataIn,XIn, XLow, XHigh,SIn, SLow, SHigh):
    Temp=[] 
    with open(DataIn, 'r') as file:
        for record in file:
#            line=record.split(',')
            line=re.split(',|\t',record)            
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    Return=aData[:,1]
    Time=aData[:,0]
    Temp=[] 
    with open(XIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    X=aData[:]    
    
    Temp=[] 
    with open(XHigh, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line[2])
    aData=np.array(Temp, dtype=np.float64)       
    XH=aData[:]    

    Temp=[] 
    with open(XLow, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line[2])
    aData=np.array(Temp, dtype=np.float64)       
    XL=aData[:]        
    
    Temp=[] 
    with open(SIn, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line)
    aData=np.array(Temp, dtype=np.float64)       
    S=aData[:]
    
    Temp=[] 
    with open(SLow, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line[2])
    aData=np.array(Temp, dtype=np.float64)       
    SL=aData[:]
  
    Temp=[] 
    with open(SHigh, 'r') as file:
        for record in file:
            line=record.split(',')
            size=len(line)
            line[size-1]=line[size-1].strip()
            Temp.append(line[2])
    aData=np.array(Temp, dtype=np.float64)       
    SH=aData[:]  
        
    delta = [ datetime.timedelta(seconds=(s+3600*9.5) ) for s in Time] 
        
    year=int(DataIn[5:9])
    month=int(DataIn[9:11])
    if(len(DataIn.split('_')[1])==8):  
        day=int(DataIn[11:13])
    else:
        day=int(DataIn[11:12])
            
            
    dt=[]
    index=[]
    index.append(0)
    days=[]
    days.append(datetime.date(year, month, day).strftime('%d/%m'))
    PrevTime=0
    jump=0
    for k in range(0,np.size(Time)):
        if Time[k]<PrevTime:
            if jump==0:
                day=day+3
            else:
                day=day+1
            jump=jump+1
            days.append(datetime.date(year, month, day).strftime('%d/%m'))
            index.append(k)
           
                
         
        dt.append(datetime.datetime(year, month, day,0,0,0)+delta[k])
        PrevTime=Time[k]
        
    t = np.array([ dates.date2num(l) for l in dt ] )   
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    lns1=ax1.plot(np.arange(np.size(Time)), Return, linewidth=1,  label='tick retruns', color='green')
    ax1.set_xticks(index)
    ax1.set_xticklabels(days)
    ax1.tick_params(axis='both', which='major',labelsize=7 )
    ax1.set_xlim(1, np.size(Time))
    ax1.grid()
    lns=lns1
    labs=[l.get_label() for l in lns]        
    ax1.legend(lns,labs, loc=0,prop={'size':8})
    ax2 = fig.add_subplot(312)
    lns2=ax2.plot(np.arange(np.size(Time)), X, linewidth=0.5,  label='$x_t$')
#    lns2l=ax2.plot(np.arange(np.size(Time)), XL, linewidth=1, color='red' , alpha=0, label='$x_t$')
#    lns2h=ax2.plot(np.arange(np.size(Time)), XH, linewidth=1, color='red' ,  alpha=0, label='$x_t$')
    fbk = {'lw':0.0, 'edgecolor':None}
    ax2.fill_between(np.arange(np.size(Time)), XL,XH, facecolor='blue', alpha=0.5,**fbk)
    ax2.set_xticks(index)
    ax2.set_xticklabels(days) 
    ax2.tick_params(axis='both', which='major',labelsize=7 )
    ax2.set_xlim(1, np.size(Time))
    ax2.grid()
    lns=lns2
    labs=[l.get_label() for l in lns]        
    ax2.legend(lns,labs, loc=0,prop={'size':8})
    ax3 = fig.add_subplot(313)
    lns3=ax3.plot(np.arange(np.size(Time)), S, linewidth=1,  label='$s_t$')
#    lns3l=ax3.plot(np.arange(np.size(Time)), SL, linewidth=1, color='red' ,  label='$s_t$')
#    lns3h=ax3.plot(np.arange(np.size(Time)), SH, linewidth=1, color='red' ,   label='$s_t$')
    ax3.fill_between(np.arange(np.size(Time)), SL,SH, facecolor='blue', alpha=0.3,**fbk)
    ax3.set_xticks(index)
    ax3.set_xticklabels(days)    
    ax3.tick_params(axis='both', which='major',labelsize=7 )   
    ax3.set_xlim(1, np.size(Time))
    ax3.grid()
    lns=lns3
    labs=[l.get_label() for l in lns]        
    ax3.legend(lns,labs, loc=0,prop={'size':8})    
    fig.suptitle('Volatility decompostion of '+ DataIn.split('_')[3].strip('.csv')+' tick returns from 23rd to 29th April 2010') 
    
    fig.tight_layout()  
    fig.subplots_adjust(top=0.92, wspace=0.35, right=0.93,left=0.07, hspace=0.2)   
    fig.savefig('../Paper/Pictures/'+DataIn.split('_')[3].strip('.csv')+'_season2.pdf', bbox_inches=0)
    fig.savefig('../Paper/Pictures/'+DataIn.split('_')[1][0:4]+'_'+DataIn.split('_')[3].strip('.csv')+'_season2.png', bbox_inches=0)   
    
''' XXXXXXXXXXXXXXXXXXX 2010 XXXXXXXXXXXXXXXXXXXXXXX '''        

#''' XXXXXXXXXXXXXXXXXXX 2010 Parameters XXXXXXXXXXXXXXXXXXXXXXX '''  
#Files=['2010_AA_Sk_EstimationResults.csv', '2010_AA_OrdNormal_EstimationResults.csv','2010_AA_DNB_EstimationResults.csv','2010_AA_OrdT_EstimationResults.csv',
#       '2010_F_Sk_EstimationResults.csv', '2010_F_OrdNormal_EstimationResults.csv','2010_F_DNB_EstimationResults.csv','2010_F_OrdT_EstimationResults.csv',
#       '2010_IBM_Sk_EstimationResults.csv', '2010_IBM_OrdNormal_EstimationResults.csv','2010_IBM_DNB_EstimationResults.csv','2010_IBM_OrdT_EstimationResults.csv',
#       '2010_JPM_Sk_EstimationResults.csv', '2010_JPM_OrdNormal_EstimationResults.csv','2010_JPM_DNB_EstimationResults.csv','2010_JPM_OrdT_EstimationResults.csv',
#       '2010_KO_Sk_EstimationResults.csv', '2010_KO_OrdNormal_EstimationResults.csv','2010_KO_DNB_EstimationResults.csv','2010_KO_OrdT_EstimationResults.csv',
#       '2010_XRX_Sk_EstimationResults.csv', '2010_XRX_OrdNormal_EstimationResults.csv','2010_XRX_DNB_EstimationResults.csv','2010_XRX_OrdT_EstimationResults.csv']  
#       
#
#ProcessMCMCOutPut(Files,20000)
#
#
#''' XXXXXXXXXXXXXXXXXXX 2010 Out of sample XXXXXXXXXXXXXXXXXXXXXXX '''  
#
#SkellamFiles=['2010_AA_Sk_vLogLikeSkellam.csv',
#              '2010_F_Sk_vLogLikeSkellam.csv',
#              '2010_IBM_Sk_vLogLikeSkellam.csv',
#              '2010_JPM_Sk_vLogLikeSkellam.csv',
#              '2010_KO_Sk_vLogLikeSkellam.csv',
#              '2010_XRX_Sk_vLogLikeSkellam.csv']
#DNBFiles=['2010_AA_DNB_vLogLikeDNB.csv',
#          '2010_F_DNB_vLogLikeDNB.csv',
#          '2010_IBM_DNB_vLogLikeDNB.csv',
#          '2010_JPM_DNB_vLogLikeDNB.csv',
#          '2010_KO_DNB_vLogLikeDNB.csv',
#          '2010_XRX_DNB_vLogLikeDNB.csv']
#OrdNFiles=['2010_AA_OrdNormal_vLogLikeNormal.csv',
#          '2010_F_OrdNormal_vLogLikeNormal.csv',
#          '2010_IBM_OrdNormal_vLogLikeNormal.csv',
#          '2010_JPM_OrdNormal_vLogLikeNormal.csv',
#          '2010_KO_OrdNormal_vLogLikeNormal.csv',
#          '2010_XRX_OrdNormal_vLogLikeNormal.csv']
#OrdTFiles=['2010_AA_OrdT_vLogLikeT.csv',
#          '2010_F_OrdT_vLogLikeT.csv',
#          '2010_IBM_OrdT_vLogLikeT.csv',
#          '2010_JPM_OrdT_vLogLikeT.csv',
#          '2010_KO_OrdT_vLogLikeT.csv',
#          '2010_XRX_OrdT_vLogLikeT.csv']          
#          
#OutSampleFiles=['Data_20100430_30_AA.csv',
#                'Data_20100430_30_F.csv',
#                'Data_20100430_30_IBM.csv',
#                'Data_20100430_30_JPM.csv',
#                'Data_20100430_30_KO.csv',
#                'Data_20100430_30_XRX.csv' ]
#
#ProcessLLOutPut(SkellamFiles,DNBFiles,OrdNFiles,OrdTFiles,OutSampleFiles)



''' XXXXXXXXXXXXXXXXXXX 2010 in sample XXXXXXXXXXXXXXXXXXXXXXX '''  

#SkellamFiles=['2010_AA_Sk_vInSkellamBIC.csv',
#              '2010_F_Sk_vInSkellamBIC.csv',
#              '2010_IBM_Sk_vInSkellamBIC.csv',
#              '2010_JPM_Sk_vInSkellamBIC.csv',
#              '2010_KO_Sk_vInSkellamBIC.csv',
#              '2010_XRX_Sk_vInSkellamBIC.csv']
#DNBFiles=['2010_AA_DNB_vInDNBBIC.csv',
#          '2010_F_DNB_vInDNBBIC.csv',
#          '2010_IBM_DNB_vInDNBBIC.csv',
#          '2010_JPM_DNB_vInDNBBIC.csv',
#          '2010_KO_DNB_vInDNBBIC.csv',
#          '2010_XRX_DNB_vInDNBBIC.csv']
#
#NormalFiles=['2010_AA_OrdNormal_vInNormalBIC.csv',
#              '2010_F_OrdNormal_vInNormalBIC.csv',
#              '2010_IBM_OrdNormal_vInNormalBIC.csv',
#              '2010_JPM_OrdNormal_vInNormalBIC.csv',
#              '2010_KO_OrdNormal_vInNormalBIC.csv',
#              '2010_XRX_OrdNormal_vInNormalBIC.csv']
#TFiles=['2010_AA_OrdT_vTBIC.csv',
#          '2010_F_OrdT_vTBIC.csv',
#          '2010_IBM_OrdT_vTBIC.csv',
#          '2010_JPM_OrdT_vTBIC.csv',
#          '2010_KO_OrdT_vTBIC.csv',
#          '2010_XRX_OrdT_vTBIC.csv']          
#          
#
#EstimatesSk=['2010_AA_Sk_EstimationResults.csv',
#       '2010_F_Sk_EstimationResults.csv',
#       '2010_IBM_Sk_EstimationResults.csv',
#       '2010_JPM_Sk_EstimationResults.csv', 
#       '2010_KO_Sk_EstimationResults.csv', 
#       '2010_XRX_Sk_EstimationResults.csv']  
#EstimatesDNB=[ '2010_AA_DNB_EstimationResults.csv',
#        '2010_F_DNB_EstimationResults.csv',
#        '2010_IBM_DNB_EstimationResults.csv',
#        '2010_JPM_DNB_EstimationResults.csv',
#        '2010_KO_DNB_EstimationResults.csv',
#        '2010_XRX_DNB_EstimationResults.csv']         
#
#EstimatesNormal=['2010_AA_OrdNormal_EstimationResults.csv',
#       '2010_F_OrdNormal_EstimationResults.csv',
#       '2010_IBM_OrdNormal_EstimationResults.csv',
#       '2010_JPM_OrdNormal_EstimationResults.csv',
#       '2010_KO_OrdNormal_EstimationResults.csv',
#       '2010_XRX_OrdNormal_EstimationResults.csv']  
#
#EstimatesT=['2010_AA_OrdT_EstimationResults.csv',
#       '2010_F_OrdT_EstimationResults.csv',
#       '2010_IBM_OrdT_EstimationResults.csv',
#       '2010_JPM_OrdT_EstimationResults.csv',
#       '2010_KO_OrdT_EstimationResults.csv',
#       '2010_XRX_OrdT_EstimationResults.csv']  
#
#InSampleFiles=['Data_20100423_29_AA.csv',
#                'Data_20100423_29_F.csv',
#                'Data_20100423_29_IBM.csv',
#                'Data_20100423_29_JPM.csv',
#                'Data_20100423_29_KO.csv',
#                'Data_20100423_29_XRX.csv' ]
#
#
#ProcessBICOutPut(SkellamFiles,DNBFiles,NormalFiles,TFiles,EstimatesSk, EstimatesDNB,EstimatesNormal, EstimatesT,InSampleFiles)
'''
========================================================================================
========================================================================================
========================================================================================'''

#''' XXXXXXXXXXXXXXXXXXX 2008 XXXXXXXXXXXXXXXXXXXXXXX '''        
#
#''' XXXXXXXXXXXXXXXXXXX 2008 Parameters XXXXXXXXXXXXXXXXXXXXXXX '''  
#Files=['2008_AA_Sk_EstimationResults.csv', '2008_AA_OrdNormal_EstimationResults.csv','2008_AA_DNB_EstimationResults.csv','2008_AA_OrdT_EstimationResults.csv',
#       '2008_F_Sk_EstimationResults.csv', '2008_F_OrdNormal_EstimationResults.csv','2008_F_DNB_EstimationResults.csv','2008_F_OrdT_EstimationResults.csv',
#       '2008_IBM_Sk_EstimationResults.csv', '2008_IBM_OrdNormal_EstimationResults.csv','2008_IBM_DNB_EstimationResults.csv','2008_IBM_OrdT_EstimationResults.csv',
#       '2008_JPM_Sk_EstimationResults.csv', '2008_JPM_OrdNormal_EstimationResults.csv','2008_JPM_DNB_EstimationResults.csv','2008_JPM_OrdT_EstimationResults.csv',
#       '2008_KO_Sk_EstimationResults.csv', '2008_KO_OrdNormal_EstimationResults.csv','2008_KO_DNB_EstimationResults.csv','2008_KO_OrdT_EstimationResults.csv',
#       '2008_XRX_Sk_EstimationResults.csv', '2008_XRX_OrdNormal_EstimationResults.csv','2008_XRX_DNB_EstimationResults.csv','2008_XRX_OrdT_EstimationResults.csv']  
#       
#
#ProcessMCMCOutPut(Files,20000)


#''' XXXXXXXXXXXXXXXXXXX 2008 Out of sample XXXXXXXXXXXXXXXXXXXXXXX '''  
#
#SkellamFiles=['2008_AA_Sk_vLogLikeSkellam.csv',
#              '2008_F_Sk_vLogLikeSkellam.csv',
#              '2008_IBM_Sk_vLogLikeSkellam.csv',
#              '2008_JPM_Sk_vLogLikeSkellam.csv',
#              '2008_KO_Sk_vLogLikeSkellam.csv',
#              '2008_XRX_Sk_vLogLikeSkellam.csv']
#DNBFiles=['2008_AA_DNB_vLogLikeDNB.csv',
#          '2008_F_DNB_vLogLikeDNB.csv',
#          '2008_IBM_DNB_vLogLikeDNB.csv',
#          '2008_JPM_DNB_vLogLikeDNB.csv',
#          '2008_KO_DNB_vLogLikeDNB.csv',
#          '2008_XRX_DNB_vLogLikeDNB.csv']
#OrdNFiles=['2008_AA_OrdNormal_vLogLikeNormal.csv',
#          '2008_F_OrdNormal_vLogLikeNormal.csv',
#          '2008_IBM_OrdNormal_vLogLikeNormal.csv',
#          '2008_JPM_OrdNormal_vLogLikeNormal.csv',
#          '2008_KO_OrdNormal_vLogLikeNormal.csv',
#          '2008_XRX_OrdNormal_vLogLikeNormal.csv']
#OrdTFiles=['2008_AA_OrdT_vLogLikeT.csv',
#          '2008_F_OrdT_vLogLikeT.csv',
#          '2008_IBM_OrdT_vLogLikeT.csv',
#          '2008_JPM_OrdT_vLogLikeT.csv',
#          '2008_KO_OrdT_vLogLikeT.csv',
#          '2008_XRX_OrdT_vLogLikeT.csv']          
#          
#OutSampleFiles=['Data_20081010_10_AA.csv',
#                'Data_20081010_10_F.csv',
#                'Data_20081010_10_IBM.csv',
#                'Data_20081010_10_JPM.csv',
#                'Data_20081010_10_KO.csv',
#                'Data_20081010_10_XRX.csv' ]
#
#ProcessLLOutPut(SkellamFiles,DNBFiles,OrdNFiles,OrdTFiles,OutSampleFiles)
#
#
#
#''' XXXXXXXXXXXXXXXXXXX 2008 in sample XXXXXXXXXXXXXXXXXXXXXXX '''  
#
#SkellamFiles=['2008_AA_Sk_vInSkellamBIC.csv',
#              '2008_F_Sk_vInSkellamBIC.csv',
#              '2008_IBM_Sk_vInSkellamBIC.csv',
#              '2008_JPM_Sk_vInSkellamBIC.csv',
#              '2008_KO_Sk_vInSkellamBIC.csv',
#              '2008_XRX_Sk_vInSkellamBIC.csv']
#DNBFiles=['2008_AA_DNB_vInDNBBIC.csv',
#          '2008_F_DNB_vInDNBBIC.csv',
#          '2008_IBM_DNB_vInDNBBIC.csv',
#          '2008_JPM_DNB_vInDNBBIC.csv',
#          '2008_KO_DNB_vInDNBBIC.csv',
#          '2008_XRX_DNB_vInDNBBIC.csv']
#
#NormalFiles=['2008_AA_OrdNormal_vInNormalBIC.csv',
#              '2008_F_OrdNormal_vInNormalBIC.csv',
#              '2008_IBM_OrdNormal_vInNormalBIC.csv',
#              '2008_JPM_OrdNormal_vInNormalBIC.csv',
#              '2008_KO_OrdNormal_vInNormalBIC.csv',
#              '2008_XRX_OrdNormal_vInNormalBIC.csv']
#TFiles=['2008_AA_OrdT_vTBIC.csv',
#          '2008_F_OrdT_vTBIC.csv',
#          '2008_IBM_OrdT_vTBIC.csv',
#          '2008_JPM_OrdT_vTBIC.csv',
#          '2008_KO_OrdT_vTBIC.csv',
#          '2008_XRX_OrdT_vTBIC.csv']          
#          
#
#EstimatesSk=['2008_AA_Sk_EstimationResults.csv',
#       '2008_F_Sk_EstimationResults.csv',
#       '2008_IBM_Sk_EstimationResults.csv',
#       '2008_JPM_Sk_EstimationResults.csv', 
#       '2008_KO_Sk_EstimationResults.csv', 
#       '2008_XRX_Sk_EstimationResults.csv']  
#EstimatesDNB=[ '2008_AA_DNB_EstimationResults.csv',
#        '2008_F_DNB_EstimationResults.csv',
#        '2008_IBM_DNB_EstimationResults.csv',
#        '2008_JPM_DNB_EstimationResults.csv',
#        '2008_KO_DNB_EstimationResults.csv',
#        '2008_XRX_DNB_EstimationResults.csv']         
#
#EstimatesNormal=['2008_AA_OrdNormal_EstimationResults.csv',
#       '2008_F_OrdNormal_EstimationResults.csv',
#       '2008_IBM_OrdNormal_EstimationResults.csv',
#       '2008_JPM_OrdNormal_EstimationResults.csv',
#       '2008_KO_OrdNormal_EstimationResults.csv',
#       '2008_XRX_OrdNormal_EstimationResults.csv']  
#
#EstimatesT=['2008_AA_OrdT_EstimationResults.csv',
#       '2008_F_OrdT_EstimationResults.csv',
#       '2008_IBM_OrdT_EstimationResults.csv',
#       '2008_JPM_OrdT_EstimationResults.csv',
#       '2008_KO_OrdT_EstimationResults.csv',
#       '2008_XRX_OrdT_EstimationResults.csv']  
#
#InSampleFiles=['Data_2008103_9_AA.csv',
#                'Data_2008103_9_F.csv',
#                'Data_2008103_9_IBM.csv',
#                'Data_2008103_9_JPM.csv',
#                'Data_2008103_9_KO.csv',
#                'Data_2008103_9_XRX.csv' ]
#
#
#ProcessBICOutPut(SkellamFiles,DNBFiles,NormalFiles,TFiles,EstimatesSk, EstimatesDNB,EstimatesNormal, EstimatesT,InSampleFiles)

''' XXXXXXXXXXXXXXXXXXX Sesonal XXXXXXXXXXXXXXXXXXXXXXX '''  
#Sesonal('Data_20100423_29_F.csv','2010_F_OrdT_vXEst.csv', '2010_F_OrdT_vSeasonEst.csv','2010_F_OrdT_vLogIntEst.csv')

#Sesonal2('Data_20100423_29_IBM.csv','2010_IBM_OrdT_vXEst.csv','2010_IBM_OrdT_vXLow.csv','2010_IBM_OrdT_vXHigh.csv', '2010_IBM_OrdT_vSeasonEst.csv','2010_IBM_OrdT_vSLow.csv','2010_IBM_OrdT_vSHigh.csv')
#Sesonal2('Data_2008103_9_IBM.csv','2008_IBM_DNB_vXEst.csv','2008_IBM_DNB_vXLow.csv','2008_IBM_DNB_vXHigh.csv', '2008_IBM_DNB_vSeasonEst.csv','2008_IBM_DNB_vSLow.csv','2008_IBM_DNB_vSHigh.csv')
#ProcessLLOutPut(['2010_IBM_Sk_vLogLikeSkellam.csv'],['2010_IBM_DNB_vLogLikeDNB.csv'], ['Data_20100430_30_IBM.csv'] )
#ProcessBICOutPut(['2010_IBM_Sk_vInSkellamBIC.csv'],['2010_IBM_DNB_vInDNBBIC.csv'],['2010_IBM_Sk_EstimationResults.csv'], ['2010_IBM_DNB_EstimationResults.csv'],['Data_20100423_29_IBM.csv'])
