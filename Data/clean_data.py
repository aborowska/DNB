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
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from scipy.stats import skew

plt.rc('font', family='serif')
#
#def convert_time(X):
#       h,m,s=re.split(':',X) # split the string to get time
#       return datetime.timedelta(hours=int(h), minutes=int(m),seconds=float(s)).total_seconds() # get time expressed in "seconds"
#
#def offset(time, offset):
#     h,m,s=re.split(':',time) # split the string to get time 
#     sec,microsec=s.split('.')
#     return datetime.time(int(h)+int(offset),int(m), int(sec) ,int(microsec))
#
#
#def offset_and_convert_time(time, offset):
#     h,m,s=re.split(':',time) # split the string to get time    
#     return datetime.timedelta(hours=int(h)+int(offset), minutes=int(m),seconds=float(s)).total_seconds()


def kde_sklearn(x, x_grid, bandwidth, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

def datetime_offset(date,time, offset):
     y=str(date)[0:4]
     mon=str(date)[4:6]
     day=str(date)[6:8]
     h,m,s=re.split(':',time) # split the string to get time 
     sec,microsec=s.split('.')
     current=datetime.datetime(int(y),int(mon),int(day) , int(h),int(m), int(sec) ,int(microsec))        
     delta=datetime.timedelta(hours=offset)
     return delta+current
     
def datetime_date(date):
     y=str(date)[0:4]
     mon=str(date)[4:6]
     day=str(date)[6:8]
     current=datetime.datetime(int(y),int(mon),int(day) , 0,0, 0 ,0)        
     return current
     
def datetime_offset_milisec(date,time, offset):
     y=str(date)[0:4]
     mon=str(date)[4:6]
     day=str(date)[6:8]
     h,m,s=re.split(':',time) # split the string to get time 
     sec,milisec=s.split('.')
     current=datetime.datetime(int(y),int(mon),int(day) , int(h),int(m), int(sec) ,int(milisec)*1000)        
     delta=datetime.timedelta(hours=offset)
     return delta+current     

def time_from_datetime(time): 
    return datetime.time(time.hour,time.minute,time.second ,time.microsecond)
   
def time_in_sec(time):
    return datetime.timedelta(hours=time.hour,minutes=time.minute,seconds=time.second ,microseconds=time.microsecond).total_seconds()  
    
def time_in_sec_from_start(time):
    time=datetime.timedelta(hours=time.hour,minutes=time.minute,seconds=time.second ,microseconds=time.microsecond).total_seconds()
    start=datetime.timedelta(hours=9,minutes=30,seconds=0 ,microseconds=0).total_seconds()
    return time-start

def get_day(time):
    return time.day


def autocorr(x, t=1):
    return np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))[1,0]


def CheckData(DataFiles):
    #loading data
    for f in DataFiles: 
        print 'XXXXXXXXX '+f+' XXXXXXXXXXX'
        #Reading in the file
        df=DataFrame.from_csv('RawData/'+f).reset_index()
        
    
        #Dropping unnecessary columns 
        #df=df[['#RIC' , 'Date[G]', 'Time[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Exch Time']]
        df=df[['#RIC' , 'Date[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Time[G]','Exch Time']]
        #df=df[['#RIC' , 'Date[G]', 'Time[G]','GMT Offset', 'Type','Price', 'Volume']]
        # Forward filling missing values fro Bid Price and Ask Price 

        print 'Number of rows: '+str(len(df))
       
        #Filtering on trades
        trades=df[df['Type']=='Trade']
        print 'Number of trades: '+str(len(trades))
        
        trades=trades[(np.isnan(trades['Price'])==0) & ( np.isnan(trades['Volume'])==0) ]
        print 'Number of trades without price: '+str(len(trades))

        
        #trades=trades[:5000]


        trades_timeg=trades.sort(['Date[G]','Time[G]'])
        trades_timee=trades.sort(['Date[G]','Exch Time'])       
        
        price_timeg=np.array(trades_timeg['Price'])
        price_timee=np.array(trades_timee['Price'])
        
        price_diff=price_timeg-price_timee
        sum_diff=np.sum(price_diff)
        print sum_diff
       
       
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

def DataTable(Table, RowNames, Ticks,FileName):
    iNumOfCols=2*len(Ticks)+1
    print 'Col'
    print iNumOfCols
    iNumOfRows=len(RowNames)+2
    print 'Row'
    print iNumOfRows
    iNumOfDig=3      
    
    sTable='\\begin{tabular}{'
    for i in range(0, iNumOfCols):
        sTable=sTable+'r'
    sTable=sTable+'} \\toprule \n ' 

    for i in range(0,iNumOfRows):
        for j in range(0,iNumOfCols):  
            #first row            
            if(i==0):
                if( (j>0) & (j<=len(Ticks))):
                    sTable+='& \\multicolumn{2}{c}{'+Ticks[j-1].strip('.N')+'}'
            #second row        
            elif(i==1):
                if((j>0) & (j%2==1) ):
                    sTable+='&  \\multicolumn{1}{c}{$ \\# $}'  
                if((j>0) & (j%2==0) ):
                    sTable+='&  \\multicolumn{1}{c}{$ \\% $ dropped}'
            
            #other rows
            else:
                if(j==0):
                     sTable+=RowNames[i-2]
                else:
                    if(j%2==1):
                        sTable+='&'+FormatNum(Table[(j-1)/2][i-2],iNumOfDig)
                    else:
                        if(i>2):
                            sTable+='&'+FormatNum((1-float(Table[(j-1)/2][i-2])/float(Table[(j-1)/2][i-3]))*100,2)
                        else:
                            sTable+='& '
        
        if i==iNumOfRows-1:
            sTable=sTable+' \\\\ \\bottomrule '
        elif i==1:
            sTable=sTable+' \\\\ \\midrule \n '
        elif i==0:
            sTable=sTable+' \\\\ \n'
            for j in range(0,iNumOfCols):
                if (j%2==1):
                   sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+2)+'} ' 
        else:            
            sTable=sTable+' \\\\ \n'
                     
                 
    sTable=sTable+'\n\\end{tabular}'

    f = open('../Paper/Tables/data'+FileName+ '_table.tex', "w")
    f.write(sTable)
    f.close() 
       
def DescriptiveLatexTable2(DataIn, DataOut, RowNames, Ticks,FileName,TickPerRow):
    
    iNumOfTables=int(np.ceil(float(len(Ticks))/TickPerRow))
    print iNumOfTables
    print len(Ticks)
    print TickPerRow
    print min(len(Ticks),TickPerRow)
    iNumOfCols=2*min(len(Ticks),TickPerRow)+1
    print 'Col'
    print iNumOfCols
    
    iNumOfRows=len(RowNames)+2
    print 'Row'
    print iNumOfRows
    iNumOfDig=3      
    
    sTable='\\begin{tabular}{l'
    for i in range(0, iNumOfCols-1):
        sTable=sTable+'r'
    sTable=sTable+'} \\toprule \n ' 

    for k in range(0,iNumOfTables):
        for i in range(0,iNumOfRows):
            for j in range(0,iNumOfCols):  
                #first row            
                if(i==0):
                    if( (j>0) & (j<=TickPerRow)):
                        sTable+='& \\multicolumn{2}{c}{'+Ticks[k*TickPerRow +j-1]+'}'
                #second row        
                elif(i==1):
                    if((j>0) & ((k*(iNumOfCols)+j)%2==1) ):
                        sTable+='&  \\multicolumn{1}{c}{In}'  
                    if((j>0) & ((k*(iNumOfCols)+j)%2==0) ):
                        sTable+='&  \\multicolumn{1}{c}{ Out}'
                
                #other rows
                else:
                    if(j==0):
                         sTable+=RowNames[i-2]
                    else:
                        if(j%2==1):
                            if i==iNumOfRows-1:
                                sTable+='&'+FormatNum(DataOut[i-2,((k*(iNumOfCols-1)+j)-1)/2],2)
                            else:
                                sTable+='&'+FormatNum(DataIn[i-2,((k*(iNumOfCols-1)+j)-1)/2],iNumOfDig)
                        else:
                            if i==iNumOfRows-1:
                                sTable+='&'+FormatNum(DataOut[i-2,((k*(iNumOfCols-1)+j)-1)/2],2)
                            else:
                                sTable+='&'+FormatNum(DataOut[i-2,((k*(iNumOfCols-1)+j)-1)/2],iNumOfDig)
                            
            
            if i==iNumOfRows-1:
                if k==iNumOfTables-1:
                    sTable=sTable+' \\\\ \\bottomrule '
                else:
                    sTable=sTable+' \\\\ \\midrule \n '
            elif i==1:
                sTable=sTable+' \\\\ \\midrule \n '
            elif i==0:
                sTable=sTable+' \\\\ \n'
                for j in range(0,iNumOfCols):
                    if (j%2==1):
                       sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+2)+'} ' 
            else:            
                sTable=sTable+' \\\\ \n'
                     
                 
    sTable=sTable+'\n\\end{tabular}'

    f = open('../Paper/Tables/desc'+FileName+ '_table.tex', "w")
    f.write(sTable)
    f.close() 
   


def DataTable2(Table, RowNames, Ticks,FileName,TickPerRow):
    
    iNumOfTables=int(np.ceil(float(len(Ticks))/TickPerRow))
    iNumOfCols=2*min(len(Ticks),TickPerRow)+1
    print 'Col'
    print iNumOfCols
    iNumOfRows=len(RowNames)+2
    print 'Row'
    print iNumOfRows
    iNumOfDig=3      
    
    sTable='\\begin{tabular}{'
    for i in range(0, iNumOfCols):
        sTable=sTable+'r'
    sTable=sTable+'} \\toprule \n '
    
    for k in range(0,iNumOfTables):
        for i in range(0,iNumOfRows):
            for j in range(0,iNumOfCols):  
                #first row            
                if(i==0):
                    if( (j>0) & (j<=TickPerRow)):
                        sTable+='& \\multicolumn{2}{c}{'+Ticks[k*TickPerRow +j-1]+'}'
                #second row        
                elif(i==1):
                    if((j>0) & ((k*(iNumOfCols)+j)%2==1) ):
                        sTable+='&  \\multicolumn{1}{c}{$ \\# $}'  
                    if((j>0) & ((k*(iNumOfCols)+j)%2==0) ):
                        sTable+='&  \\multicolumn{1}{c}{$ \\% $ dropped}'
                
                #other rows
                else:
                    if(j==0):
                         sTable+=RowNames[i-2]
                    else:
                        if(j%2==1):
                            sTable+='&'+FormatNum(Table[((k*(iNumOfCols-1)+j)-1)/2][i-2],iNumOfDig)
                        else:
                            if(i>2):
                                sTable+='&'+FormatNum((1-float(Table[((k*(iNumOfCols-1)+j)-1)/2][i-2])/float(Table[((k*(iNumOfCols-1)+j)-1)/2][i-3]))*100,2)
                            else:
                                sTable+='& '
            
            if i==iNumOfRows-1:
                sTable=sTable+' \\\\ \\bottomrule '
            elif i==1:
                sTable=sTable+' \\\\ \\midrule \n '
            elif i==0:
                sTable=sTable+' \\\\ \n'
                for j in range(0,iNumOfCols):
                    if (j%2==1):
                       sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+2)+'} ' 
            else:            
                sTable=sTable+' \\\\ \n'
                     
                 
    sTable=sTable+'\n\\end{tabular}'

    f = open('../Paper/Tables/data'+FileName+ '_table.tex', "w")
    f.write(sTable)
    f.close() 
       

def CleanData(DataFiles,Ticks,DaysIn,DaysOut):    
    #loading data
    table=[]
    Days=DaysIn+DaysOut
    for tick in Ticks:
        data=[]
        nums=[]
        numfiles=1
        print 'XXXXXXXXXXxxxxxxxxxXXXXXXXXXX'
        print 'XXXXXXXXX '+tick+' XXXXXXXXXXX'
        print 'XXXXXXXXXXxxxxxxxxxXXXXXXXXXX'
        
        for f in DataFiles: 
            print 'XXXXXXXXX '+f+' XXXXXXXXXXX'
            #Reading in the file
            df=DataFrame.from_csv('RawData/'+f).reset_index()
            #Dropping unnecessary columns 
            #df=df[['#RIC' , 'Date[G]', 'Time[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Exch Time']]
            df=df[['#RIC' , 'Date[G]','Time[G]','GMT Offset', 'Type','Price', 'Volume','Bid Price', 'Ask Price','Exch Time']]
                        
            # Filtering on days and tick 
            df=df[df['#RIC']==tick]
            print 'Filter on Tick is done'
            #Calculating datetime and offset
            if(len(df)!=0):
                df['Date'] = df.apply(lambda x: datetime_offset(x['Date[G]'],x['Time[G]'], x['GMT Offset']), axis=1)
                # Filtering on days and tick           
                df=df[df['Date'].map(get_day).isin(Days)]  
                
            if(len(df)!=0):
                if numfiles==1:
                    data=df
                    numfiles+=1
                else:
                    data=data.append(df, ignore_index=True)
              
            
    
           
        #Calculate bid ask spread
        data['Bid Price']=data['Bid Price'].ffill()     
        data['Ask Price']=data['Ask Price'].ffill() 
            
        nums.append(len(data.index))    
        print 'Number of rows: '+str(len(data.index))
       
        
        #Filtering on trades
        trades=data[data['Type']=='Trade']
        nums.append(len(trades.index))     
        print 'Number of trades: '+str(len(trades.index))  
        
        
            
        #Delete rows with no price or volume
        trades=trades[(np.isnan(trades['Price'])==0) & ( np.isnan(trades['Volume'])==0) &(trades['Price']!=0)&(trades['Volume']!=0)]
        nums.append(len(trades.index)) 
        print 'Number of trades without price: '+str(len(trades.index))
      
        #print 'Trades per day'
        #print trades.groupby(['#RIC', 'Date[G]'],sort=True, as_index=False).agg({  'Exch Time': 'count'})
         
        #Deleting trades out side of the considered interval 
        #Calculating datetime and offset
        trades['DateTime'] = trades.apply(lambda x: datetime_offset_milisec(x['Date[G]'],x['Exch Time'], x['GMT Offset']), axis=1)
        
        #Calculating time
        trades['Time'] = trades['DateTime'].map(time_from_datetime)      
        OpeningTime=datetime.time(9,30,0)
        ClosingTime=datetime.time(16,0,0)    
        trades=trades[ (trades['Time'] >= OpeningTime) & (trades["Time"] <= ClosingTime)]   
        nums.append(len(trades.index))     
        print 'Number of trades in trading period: '+str(len(trades.index))         
        
            
        #Take the last observation in every multiple trade per microsecond      
        trades=trades.groupby(['#RIC', 'DateTime'],sort=True, as_index=False).agg({  'Time': 'last',  'Price': 'last', 'Bid Price': 'last', 'Ask Price': 'last', 'Volume': np.sum})        
        nums.append(len(trades.index))             
        print 'Number of trades after aggregation: '+str(len(trades.index))    
        

     
        trades['Spread']=trades['Ask Price']-trades['Bid Price']   
        #print trades[:20]
          
        # Delete rows with price < bid - bidask or price > ask + bidask        
        trades=trades[(trades['Price']>= trades['Bid Price']- trades['Spread'] ) & (trades['Price']<= trades['Ask Price']+ trades['Spread'] ) ]
        nums.append(len(trades.index))             
        print 'Number of trades in bid ask range: '+str(len(trades.index))
            
            
           
           
        trades=trades[['#RIC' , 'DateTime', 'Time',  'Price', 'Bid Price','Ask Price','Volume']]
        
        trades['Duration']=trades['Time'].map(time_in_sec).diff()
        trades['TimeInSec']=trades['Time'].map(time_in_sec_from_start)
        #Calculating tick returns
        trades['TickReturn']=(trades['Price']/0.01).diff()   
        trades['LogReturn']=(np.log(trades['Price'])).diff() 
        trades=trades[trades['Duration'] >=0 ]  
        nums.append(len(trades.index))             
        print 'Number of trades in without opening trades: '+str(len(trades.index))        
        
        #Float formatting 
        trades['TickReturn']= trades['TickReturn'].apply(lambda x: '%d' % x)    
        trades['Volume']= trades['Volume'].apply(lambda x: '%d' % x)
        trades['Price']= trades['Price'].apply(lambda x: '%.2f' % x)
        trades['Duration']= trades['Duration'].apply(lambda x: '%.6f' % x)
        trades['TimeInSec']= trades['TimeInSec'].apply(lambda x: '%.6f' % x)
        
        
                     
        print 'Befor writing out: '+str(len(trades.index))
       
        
        trades.to_csv('DataFull_'+f[5:11]+str(Days[0])+'_'+str(Days[len(Days)-1])+'_'+tick.rstrip('.N')+'.csv',index=False)
       
        tradesIn=trades[trades['DateTime'].map(get_day).isin(DaysIn)]
        tradesIn=tradesIn[['TimeInSec', 'TickReturn','Price','LogReturn']]
        tradesIn.to_csv('Data_'+f[5:11]+str(DaysIn[0])+'_'+str(DaysIn[len(DaysIn)-1])+'_'+tick.rstrip('.N')+'.csv' ,index=False,header=False)
        
        tradesOut=trades[trades['DateTime'].map(get_day).isin(DaysOut)]
        tradesOut=tradesOut[['TimeInSec', 'TickReturn','Price','LogReturn']]
        tradesOut.to_csv('Data_'+f[5:11]+str(DaysOut[0])+'_'+str(DaysOut[len(DaysOut)-1])+'_'+tick.rstrip('.N')+'.csv' ,index=False,header=False)
        
        table.append(nums)
        
    RowNames=['Raw quotes and trades','Trades', 'Non missing price and volume', 'Trades between 9:30 and 16:00', 'Aggrageted trades', 'Without outliers', 'Without opening trades' ]                
        
    with open('table'+DataFiles[0][5:11]+'.csv', "wb") as f:
        writer = csv.writer(f)
        writer.writerows(table)    
        
    DataTable(table, RowNames, Ticks,DataFiles[0][5:11])
        
   
def AddColumnName(DataFiles,ColNameFile):
    
    with open('RawData/'+ColNameFile, 'r') as file:
        for record in file:
            names=record.split(',')
            size=len(names)
            names[size-1]=names[size-1].strip()

    
    print names
    for f in DataFiles:
        Data=[]
        Data.append(names)
        with open('RawData/'+f, 'r') as file:
            for record in file: 
                line = record.split(',')
                size=len(line)
                line[size-1]=line[size-1].strip()
                Data.append(line)
                
        with open('RawData/'+'new_'+f, 'wb') as fout:
            writer = csv.writer(fout)
            writer.writerows(Data)             
               
        print f+' is done!'
        
def SkellamDens(x, l):
    return np.exp(-2*l)*special.iv(np.abs(x), 2*l)  

def ExponentialDens(x,l):
    return l*np.exp(-l*x) 

def GaussianDens(x,mu,sigma):
    u=(x-mu)/sigma
    y=(1/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*u*u)
    return y

def OrdNormDens(x,sigma):
    return norm.cdf((x+0.5)/sigma)-norm.cdf((x-0.5)/sigma) 

def SkellamVSGaussian(b, var):
    x=np.arange(-b,b+1,1)
    gaussian=[np.log(SkellamDens(i,var/2.0)) for i in x  ]
    skellam=[np.log(GaussianDens(i,0,np.sqrt(var)) )for i in x  ]  
    ordnorm=[np.log( OrdNormDens(i,np.sqrt(var)) )for i in x  ]
    print ordnorm
#    gaussian=[SkellamDens(i,var/2.0) for i in x  ]
#    skellam =[GaussianDens(i,0,np.sqrt(var)) for i in x  ]
    fig = plt.figure()
    ax1 = fig.add_subplot(111) 
    lns1=ax1.plot(x,skellam,label='Skellam', linestyle='-', marker='o',color='blue')
    lns2=ax1.plot(x,gaussian,label='Gaussian', linestyle='-', marker='x',color='red')
    lns3=ax1.plot(x,ordnorm,label='Ordered normal', linestyle='-', marker='s',color='green')
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    #ax1.legend(lns, labs, loc=2, prop={'size':8})
    ax1.legend(lns, labs, loc=2)
    plt.title('Discrete distributions', fontsize=9)

    fig.savefig('skellam_vs_gaussian.png', bbox_inches=0)
    fig.savefig('skellam_vs_gaussian.pdf', bbox_inches=0)      

def EmpiricalLogReturnHist(DataFiles,Bound):
    
    startPosition=np.ceil(len(DataFiles)/2.0)*100+20+1  
    fig = plt.figure()
    i=0
    for f in DataFiles:
        # years month and ticker
        Year=f[0:4]        
        Month=f[4:6]        
        Ticker=f[7:].strip('.csv')
        # load data to data frame
        df=DataFrame.from_csv(f).reset_index()
        x=np.array(df['LogReturn'])
        x=x[(x>-Bound[i] )& (x<Bound[i])]
        #Histogram
       
        xmin=min(x)
        xmax=max(x)
      
        ax1 = fig.add_subplot(startPosition)
        iNumOfBins=np.round(np.sqrt(len(x)))
        lns1=ax1.hist(x, iNumOfBins, fc='blue', histtype='stepfilled', normed=True)
        ax1.tick_params(axis='both', which='major', labelsize=7)  
        ax1.set_xlim([xmin,xmax])  
        ymin=ax1.get_ylim()[0]
        ymax=ax1.get_ylim()[1]
        ax1.set_ylim([ymin,0.1*ymax])  
        plt.title(Ticker, fontsize=9)
      
        
        startPosition+=1
        i+=1
        
    fig.suptitle('Empirical distribution of the log returns in '+Month+'/'+Year)     
    fig.tight_layout()  
    fig.subplots_adjust(top=0.9)           
    fig.savefig('hist_logreturn_'+Year+Month+'.png', bbox_inches=0)
    fig.savefig('hist_logreturn_'+Year+Month+'.pdf', bbox_inches=0)      

def EmpiricalTickReturnHist(DataFiles):
    startPosition=np.ceil(len(DataFiles)/2.0)*100+20+1  
   
    fig = plt.figure()
    i=0
    for f in DataFiles:
        # years month and ticker
        Year=f[0:4]        
        Month=f[4:6]        
        Ticker=f[7:].strip('.csv')
        # load data to data frame
        df=DataFrame.from_csv(f).reset_index()
        #Histogram
        hist=df.groupby('TickReturn').size().reset_index()
        hist=hist.applymap(int).sort('TickReturn')   
        x=np.array(hist.ix[:,0])
        y=np.array(hist.ix[:,1],dtype=float )/len(df)
        # Skellam
        x2=np.square(x)
        var= np.dot(x2, y)
        lamb=var/2;
        y_skellam=[SkellamDens(i, lamb) for i in x]
        ax1 = fig.add_subplot(startPosition)
        lns1=ax1.plot(x,y,label='empirical', linestyle='-', marker='o',color='blue')
        lns2=ax1.plot(x,y_skellam,label='skellam', linestyle='-.', marker='+',color='red')
        ax1.set_yscale('log')
        ymin=np.amin(y) 
        ymax=np.amax(y) 
        ax1.set_ylim([ymin,ymax])
        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        #ax1.legend(lns, labs, loc=2, prop={'size':8})
        ax1.legend(lns, labs, loc=2)
        ax1.tick_params(axis='both', which='major', labelsize=7)
        plt.ylabel('log density',fontsize=7)   
        ax1.legend(lns,labs, loc=2,prop={'size':7})
        #plt.xlabel('tick returns')    
        plt.title(Ticker, fontsize=9)
        startPosition+=1
        i+=1
        
    fig.suptitle('Empirical distribution of the tick returns in '+Month+'/'+Year)     
    fig.tight_layout()  
    fig.subplots_adjust(top=0.9)           
    fig.savefig('hist_tickreturn_'+Year+Month+'.png', bbox_inches=0)
    fig.savefig('hist_tickreturn_'+Year+Month+'.pdf', bbox_inches=0)      


def EmpiricalLogVsTick(DataFile,Bound):
    fig = plt.figure()

    df=DataFrame.from_csv(DataFile).reset_index()    
    x=np.array(df['LogReturn'])
    x=x[(x>-Bound )& (x<Bound)]
    #Histogram
       
    xmin=min(x)
    xmax=max(x)
      
    ax1 = fig.add_subplot(121)
    iNumOfBins=np.round(np.sqrt(len(x)))
    lns1=ax1.hist(x, iNumOfBins, fc='blue', histtype='stepfilled', normed=True)
    ax1.tick_params(axis='both', which='major', labelsize=7)  
    ax1.set_xlim([xmin,xmax])  
    ymin=ax1.get_ylim()[0]
    ymax=ax1.get_ylim()[1]
    ax1.set_ylim([ymin,0.1*ymax])  
    plt.title('Log retrun', fontsize=9)

    Year=DataFile[0:4]        
    Month=DataFile[4:6]        
    Ticker=DataFile[7:].strip('.csv')
    # load data to data frame
    df=DataFrame.from_csv(DataFile).reset_index()
    #Histogram
    hist=df.groupby('TickReturn').size().reset_index()
    hist=hist.applymap(int).sort('TickReturn')   
    x=np.array(hist.ix[:,0])
    y=np.array(hist.ix[:,1],dtype=float )/len(df)
    # Skellam
    x2=np.square(x)
    var= np.dot(x2, y)
    lamb=var/2;
    y_skellam=[SkellamDens(i, lamb) for i in x]
    ax1 = fig.add_subplot(120)
    lns1=ax1.plot(x,y,label='empirical', linestyle='-', marker='o',color='blue')
    lns2=ax1.plot(x,y_skellam,label='skellam', linestyle='-.', marker='+',color='red')
    ax1.set_yscale('log')
    ymin=np.amin(y) 
    ymax=np.amax(y) 
    ax1.set_ylim([ymin,ymax])
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    #ax1.legend(lns, labs, loc=2, prop={'size':8})
    ax1.legend(lns, labs, loc=2)
    ax1.tick_params(axis='both', which='major', labelsize=7)
    plt.ylabel('log density',fontsize=7)   
    ax1.legend(lns,labs, loc=2,prop={'size':7})
    #plt.xlabel('tick returns')    
    plt.title('Tick return', fontsize=9)    
    
    
    fig.suptitle('Empirical distribution of JPM returns in '+Month+'/'+Year)     
    fig.tight_layout()  
    fig.subplots_adjust(top=0.9)           
    fig.savefig('log_vs_tick_'+Year+Month+'.png', bbox_inches=0)
    fig.savefig('log_vs_tick_'+Year+Month+'.pdf', bbox_inches=0)          
    

def EmpiricalDurationHist(DataFiles):
    
    for f in DataFiles:
        # load data to data frame
        df=DataFrame.from_csv(f).reset_index()
        x=np.array(df['Duration'])
        #Exponential
        lamb=1/np.mean(x)
        x_grid = np.linspace(min(x), max(x), 2000)
        y_expo=[ExponentialDens(i,lamb) for i in x_grid]
        #Histogram
       
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        iNumOfBins=np.round(np.sqrt(len(df)))
        lns1=ax1.hist(x, iNumOfBins, fc='blue', histtype='stepfilled', alpha=0.3, normed=True)
        lns2=ax1.plot(x_grid,y_expo,label='exponential', linestyle='-', color='red')
        lns = lns2
        labs = [l.get_label() for l in lns]
        #ax1.legend(lns, labs, loc=2, prop={'size':8})
        ax1.legend(lns, labs, loc=1)        
        xmax=np.percentile(x,95)
        ax1.set_xlim([0,xmax])
        ax1.tick_params(axis='both', which='major', labelsize=8)   
        plt.xlabel('durations')    
        plt.title('Empirical distribution of the '+f.replace('200810_', '').strip('.csv')  +' durations')
        fig.savefig('hist_dur'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('hist_dur'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0)   

def EmpiricalVolumeHist(DataFiles):
    
    for f in DataFiles:
        # load data to data frame
        df=DataFrame.from_csv(f).reset_index()
        x=np.array(df['Volume'])
        #Exponential
        lamb=1/np.mean(x)
        x_grid = np.linspace(min(x), max(x), 2000)
        y_expo=[ExponentialDens(i,lamb) for i in x_grid]
        #Histogram
       
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        iNumOfBins=np.round(np.sqrt(len(df)))
        lns1=ax1.hist(x, iNumOfBins, fc='blue', histtype='stepfilled', alpha=0.3, normed=True)
        lns2=ax1.plot(x_grid,y_expo,label='exponential', linestyle='-', color='red')
        lns = lns2
        labs = [l.get_label() for l in lns]
        #ax1.legend(lns, labs, loc=2, prop={'size':8})
        ax1.legend(lns, labs, loc=1)        
        xmax=np.percentile(x,99)
        ax1.set_xlim([0,xmax])
        ax1.tick_params(axis='both', which='major', labelsize=8)   
        plt.xlabel('volume')    
        plt.title('Empirical distribution of the '+f.replace('200810_', '').strip('.csv')  +' volume')
        fig.savefig('hist_vol'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('hist_vol'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0) 


def EmpiricalTickReturnSeries(DataFiles, Day):  
    for f in DataFiles:
        '''load data to data frame '''
        df=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()
        dtdf=to_datetime(df['DateTime'])
        dt = np.array(dtdf[dtdf.map(get_day)==Day].map(dates.date2num))
        p=df['TickReturn']
        p =  np.array(p[dtdf.map(get_day)==Day])
        p2=np.abs(p)
        
        
        '''matplotlib date format object'''
        hfmt = dates.DateFormatter(' %H:%M \n %d/%m/%y ')
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        lns1=ax1.plot(dt,p, label='Tick return', linestyle='-', color='blue')
        ax1.set_ylabel('tick return', fontsize=8 )
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax2 = ax1.twinx()  
        lns2=ax2.plot(dt,p2, label='Abs tick return', linestyle='-', color='green')
        ax2.set_ylabel('abs tick return', fontsize=8 )
        ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax2.xaxis.set_major_formatter(hfmt)
        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc=2, prop={'size':8})
        ax2.tick_params(axis='both', which='major', labelsize=8)
        plt.title(f.replace('200810_', '').strip('.csv')  +' tick returns')
        plt.xlabel('Time')
        plt.xlim((min(dt), max(dt)))
        fig.savefig('tick_'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('tick_'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0)
        

def VolAndDuration(DataFiles, Day):  
    for f in DataFiles:
        '''load data to data frame '''
        df=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()
        dtdf=to_datetime(df['DateTime'])
        dt = np.array(dtdf[dtdf.map(get_day)==Day].map(dates.date2num))
        p=df['TickReturn']
        p =  np.array(p[dtdf.map(get_day)==Day])
        p2=np.abs(p)
        d=df['Duration']
        d =  np.array(d[dtdf.map(get_day)==Day])
       
        
        
        '''matplotlib date format object'''
        hfmt = dates.DateFormatter(' %H:%M \n %d/%m/%y ')
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        lns1=ax1.plot(dt,d, label='Duration', linestyle='-', color='blue')
        ax1.set_ylabel('duration', fontsize=8 )
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax2 = ax1.twinx()  
        lns2=ax2.plot(dt,p2, label='Abs tick return', linestyle='-', color='green')
        ax2.set_ylabel('abs tick return', fontsize=8 )
        ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax2.xaxis.set_major_formatter(hfmt)
        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc=2, prop={'size':8})
        ax2.tick_params(axis='both', which='major', labelsize=8)
        plt.title(f.replace('200810_', '').strip('.csv')  +' duaration vs abs tick returns')
        plt.xlabel('Time')
        plt.xlim((min(dt), max(dt)))
        fig.savefig('dur'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('dur_'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0)        
        

def VolAndVolume(DataFiles, Day):  
    for f in DataFiles:
        '''load data to data frame '''
        df=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()
        dtdf=to_datetime(df['DateTime'])
        dt = np.array(dtdf[dtdf.map(get_day)==Day].map(dates.date2num))
        p=df['TickReturn']
        p =  np.array(p[dtdf.map(get_day)==Day])
        p2=np.abs(p)
        d=df['Volume']
        d =  np.array(d[dtdf.map(get_day)==Day])
       
        
        
        '''matplotlib date format object'''
        hfmt = dates.DateFormatter(' %H:%M \n %d/%m/%y ')
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        lns1=ax1.plot(dt,d, label='volume', linestyle='-', color='blue',alpha=0.5)
        ax1.set_ylabel('volume', fontsize=8 )
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax2 = ax1.twinx()  
        lns2=ax2.plot(dt,p2, label='Abs tick return', linestyle='-', color='green',alpha=0.5)
        ax2.set_ylabel('abs tick return', fontsize=8 )
        ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax2.xaxis.set_major_formatter(hfmt)
        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc=2, prop={'size':8})
        ax2.tick_params(axis='both', which='major', labelsize=8)
        plt.title(f.replace('200810_', '').strip('.csv')  +' volume vs abs tick returns')
        plt.xlabel('Time')
        plt.xlim((min(dt), max(dt)))
        fig.savefig('volume'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('volume_'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0)        
 

def VolPersistence(DataFiles, Day):  
    for f in DataFiles:
        '''load data to data frame '''
        df=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()
        dtdf=to_datetime(df['DateTime'])
        dt = np.array(dtdf[dtdf.map(get_day)==Day].map(dates.date2num))
        p=df['TickReturn']
        p =  np.array(p[dtdf.map(get_day)==Day])
#        p =  np.array(p)
        p2=np.square(p)
        x =range(0,101)
        
        auto=[autocorr(p2, i) for i in x]  
        print p2
        print x
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        lns1=ax1.plot(x,auto,  linestyle='-', marker='o' ,color='blue',alpha=0.75)
        ax1.set_ylabel('auto correlation', fontsize=8 )
        ax1.tick_params(axis='both', which='major', labelsize=8)
#        lns = lns1+lns2
#        labs = [l.get_label() for l in lns]
#        ax2.legend(lns, labs, loc=2, prop={'size':8})
   
        plt.title(f.replace('200810_', '').strip('.csv')  +' auto correlation of the squared tick returns')
        plt.xlabel('lag')
        fig.savefig('volper'+f.replace('200810_', '').strip('.csv') +'.png', bbox_inches=0)
        fig.savefig('volper_'+f.replace('200810_', '').strip('.csv') +'.pdf', bbox_inches=0)          
        
def PrepareData(DataFiles, Days):
    for f in DataFiles: 
        print 'XXXXXXXXX '+f+' XXXXXXXXXXX'
        df=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()
        df=df[to_datetime(df['DateTime']).map(get_day).isin(Days)]
        df=df[['TimeInSec', 'TickReturn']]
        print 'Number of observations: '+str(len(df))
        df.to_csv('Data_'+f[0:6]+str(Days[0])+'_'+str(Days[len(Days)-1])+'_'+f[7:] ,index=False,header=False)
  

def DescriptiveLatexTable(DataIn, DataOut, RowNames, Ticks,FileName):
    iNumOfCols=2*len(Ticks)+1
    print 'Col'
    print iNumOfCols
    iNumOfRows=len(RowNames)+2
    print 'Row'
    print iNumOfRows
    iNumOfDig=3      
    
    sTable='\\begin{tabular}{l'
    for i in range(0, iNumOfCols-1):
        sTable=sTable+'r'
    sTable=sTable+'} \\toprule \n ' 

    for i in range(0,iNumOfRows):
        for j in range(0,iNumOfCols):  
            #first row            
            if(i==0):
                if( (j>0) & (j<=len(Ticks))):
                    sTable+='& \\multicolumn{2}{c}{'+Ticks[j-1]+'}'
            #second row        
            elif(i==1):
                if((j>0) & (j%2==1) ):
                    sTable+='&  \\multicolumn{1}{c}{In}'  
                if((j>0) & (j%2==0) ):
                    sTable+='&  \\multicolumn{1}{c}{ Out}'
            
            #other rows
            else:
                if(j==0):
                     sTable+=RowNames[i-2]
                else:
                    if(j%2==1):
                        if i==iNumOfRows-1:
                            sTable+='&'+FormatNum(DataOut[i-2,(j-1)/2],2)
                        else:
                            sTable+='&'+FormatNum(DataIn[i-2,(j-1)/2],iNumOfDig)
                    else:
                        if i==iNumOfRows-1:
                            sTable+='&'+FormatNum(DataOut[i-2,(j-1)/2],2)
                        else:
                            sTable+='&'+FormatNum(DataOut[i-2,(j-1)/2],iNumOfDig)
                        
        
        if i==iNumOfRows-1:
            sTable=sTable+' \\\\ \\bottomrule '
        elif i==1:
            sTable=sTable+' \\\\ \\midrule \n '
        elif i==0:
            sTable=sTable+' \\\\ \n'
            for j in range(0,iNumOfCols):
                if (j%2==1):
                   sTable=sTable+'\cmidrule(r){'+str(j+1)+'-'+str(j+2)+'} ' 
        else:            
            sTable=sTable+' \\\\ \n'
                     
                 
    sTable=sTable+'\n\\end{tabular}'

    f = open('../Paper/Tables/desc'+FileName+ '_table.tex', "w")
    f.write(sTable)
    f.close() 
   



     
def DescriptiveTable(InFiles,DaysIn,DaysOut):
     iNumOfFiles=len(InFiles)
     DataIn=np.zeros((7,iNumOfFiles))
     DataOut=np.zeros((7,iNumOfFiles))
     col=0     
     Ticks=[]
     for f in InFiles:
         print f
         ''' Loading in sample and out of sample data '''    
         Data=DataFrame.from_csv(f,parse_dates=True, infer_datetime_format=True).reset_index()    
         
         In=Data[to_datetime(Data['DateTime']).map(get_day).isin(DaysIn)]
         Out=Data[to_datetime(Data['DateTime']).map(get_day).isin(DaysOut)]    
         
         DataIn[0:2,col]=np.array(In['Price'].describe())[0:2]
         DataIn[2:5,col]=np.array(In['TickReturn'].describe())[1:4]
         DataIn[5,col]=np.array(In['TickReturn'].describe())[7]
         DataIn[6,col]=100*np.array(In[In['TickReturn']==0].describe())[0,7]/DataIn[0,col]
         
         DataOut[0:2,col]=np.array(Out['Price'].describe())[0:2]
         DataOut[2:5,col]=np.array(Out['TickReturn'].describe())[1:4]
         DataOut[5,col]=np.array(Out['TickReturn'].describe())[7]
         DataOut[6,col]=100*np.array(Out[Out['TickReturn']==0].describe())[0,7]/DataOut[0,col]
         
         Ticks.append(f.split('_')[3].rstrip('.csv'))
         col=col+1
       
       
     print DataIn
     print DataOut
     print Ticks     
     RowNames=['Num. obs', 'Avg. price', 'Mean', 'Std', 'Min','Max', '\\% Zeros']
  #   DescriptiveLatexTable(DataIn, DataOut, RowNames, Ticks,InFiles[0][9:13])     
     DescriptiveLatexTable2(DataIn, DataOut, RowNames, Ticks,InFiles[0][9:13],3)
     
     
''' Adding column headers to the raw files '''     
#ColNameFile='column_names.csv'
#DataFiles=[ 'NYSE_200810_data_42.csv','NYSE_200810_data_43.csv']  
#AddColumnName(DataFiles,ColNameFile)     
#DataFiles=['NYSE_201004_data_12.csv', 
#           'NYSE_201004_data_22.csv', 
#           'NYSE_201004_data_23.csv', 
#           'NYSE_201004_data_24.csv', 
#           'NYSE_201004_data_32.csv', 
#           'NYSE_201004_data_33.csv', 
#           'NYSE_201004_data_34.csv', 
#           'NYSE_201004_data_35.csv', 
#           'NYSE_201004_data_36.csv', 
#           'NYSE_201004_data_52.csv']
#AddColumnName(DataFiles,ColNameFile)    

''' Cleaning monthly data files ''' 

#''' 2008 October ''' 
#DataFiles=['NYSE_200810_data_1.csv', 
#           'NYSE_200810_data_2.csv',
#           'NYSE_200810_data_3.csv',
#           'NYSE_200810_data_41.csv', 
#           'NYSE_200810_data_42.csv',
#           'NYSE_200810_data_43.csv',
#           'NYSE_200810_data_51.csv']
#InSampleDays=[3,6,7,8,9]
#OutSampleDays=[10]
#Ticks=['AA.N', 'F.N','IBM.N' , 'JPM.N', 'KO.N', 'XRX.N']
#CleanData(DataFiles,Ticks,InSampleDays,OutSampleDays)
#RowNames=['Raw quotes and trades','Trades', 'Non missing price and volume', 'Trades between 9:35 and 15:55', 'Aggrageted trades', 'Without outliers','Without opening trades'   ] 
#Table=np.array(DataFrame.from_csv('table200810.csv', header=None).reset_index())
#print Table
#DataTable(Table, RowNames, Ticks,DataFiles[0][5:11])
#
##''' 2010 April ''' 
#DataFiles=['NYSE_201004_data_11.csv', 
#           'NYSE_201004_data_12.csv', 
#           'NYSE_201004_data_21.csv', 
#           'NYSE_201004_data_22.csv', 
#           'NYSE_201004_data_23.csv', 
#           'NYSE_201004_data_24.csv', 
#           'NYSE_201004_data_31.csv', 
#           'NYSE_201004_data_32.csv', 
#           'NYSE_201004_data_33.csv', 
#           'NYSE_201004_data_34.csv', 
#           'NYSE_201004_data_35.csv', 
#           'NYSE_201004_data_36.csv', 
#           'NYSE_201004_data_41.csv', 
#           'NYSE_201004_data_51.csv', 
#           'NYSE_201004_data_52.csv'] 
#InSampleDays=[23, 26,27,28,29]
#OutSampleDays=[30]                    
#Ticks=['AA.N', 'F.N','IBM.N' , 'JPM.N', 'KO.N', 'XRX.N']
##Ticks=['AA.N']
#CleanData(DataFiles,Ticks,InSampleDays,OutSampleDays)
#RowNames=['Raw quotes and trades','Trades', 'Non missing price and volume', 'Trades between 9:30 and 16:00', 'Aggrageted trades', 'Without outliers', 'Without opening trades'  ] 
#Table=np.array(DataFrame.from_csv('table201004.csv', header=None).reset_index())
#print Table
#DataTable(Table, RowNames, Ticks,DataFiles[0][5:11])
#
#
#Bounds=[0.002, 0.005, 0.002,0.002, 0.002, 0.002]
#DataFiles=['200810_AA.csv','200810_F.csv','200810_IBM.csv','200810_JPM.csv','200810_KO.csv','200810_XRX.csv']
#EmpiricalTickReturnHist(DataFiles)
#EmpiricalLogReturnHist(DataFiles,Bounds)
#
#Bounds=[0.0015, 0.002, 0.0006,0.0006, 0.00125, 0.00175]
#DataFiles=['201004_AA.csv','201004_F.csv','201004_IBM.csv','201004_JPM.csv','201004_Ko.csv', '201004_XRX.csv' ]
#EmpiricalTickReturnHist(DataFiles)
#EmpiricalLogReturnHist(DataFiles,Bounds)
#
#
#
#InFiles=['DataFull_2008103_10_AA.csv','DataFull_2008103_10_F.csv','DataFull_2008103_10_IBM.csv'
#        ,'DataFull_2008103_10_JPM.csv','DataFull_2008103_10_KO.csv','DataFull_2008103_10_XRX.csv']
#InSampleDays=[3,6,7,8,9]
#OutSampleDays=[10]         
#DescriptiveTable(InFiles,InSampleDays,OutSampleDays)
#
#InFiles=['DataFull_20100423_30_AA.csv','DataFull_20100423_30_F.csv','DataFull_20100423_30_IBM.csv'
#        ,'DataFull_20100423_30_JPM.csv','DataFull_20100423_30_KO.csv','DataFull_20100423_30_XRX.csv']
#InSampleDays=[23, 26,27,28,29]
#OutSampleDays=[30]            
#DescriptiveTable(InFiles,InSampleDays,OutSampleDays)
SkellamVSGaussian(8, 1)



#EmpiricalLogVsTick('200810_JPM.csv',0.002)
''' Extra stuff.... '''
#CheckData(DataFiles)

#EmpiricalTickReturnSeries(DataFiles, 10)
#VolAndDuration(DataFiles, 10)
#VolAndVolume(DataFiles, 10)
#EmpiricalDurationHist(DataFiles) 
#EmpiricalVolumeHist(DataFiles)  
 

#VolPersistence(DataFiles, 10)  
 
##Calculating durations
#trades['Duration']=trades['Time'].map(time_in_sec).diff()
##Calculating tick returns
#trades['TickReturn']=trades['Price'].diff()  /0.01   
#trades=trades[trades['Duration'] >=0 ]    
#  
#print trades[:20]  
#dtdf=trades['DateTime']
#dt = np.array(dtdf[dtdf.map(get_day)==10].map(dates.date2num))
#p=trades['Price']
#p =  np.array(p[dtdf.map(get_day)==10])
#v=trades['Volume']
#v =  np.array(v[dtdf.map(get_day)==10])
#
#
##fs = dates.date2num(dt) # converted
#
## matplotlib date format object
#hfmt = dates.DateFormatter('%m/%d %H:%M')
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#lns1=ax1.plot(dt,p, label='Price', linestyle='-', color='blue')
#ax1.set_ylabel('price', fontsize=8 )
#ax1.tick_params(axis='both', which='major', labelsize=8)
#ax2 = ax1.twinx()  
#lns2=ax2.plot(dt,v, label='Volume', linestyle='-', color='green')
#ax2.set_ylabel('price', fontsize=8 )
#ax2.xaxis.set_major_locator(dates.HourLocator(interval=1))
#ax2.xaxis.set_major_formatter(hfmt)
#lns = lns1+lns2
#labs = [l.get_label() for l in lns]
#ax2.legend(lns, labs, loc=2, prop={'size':8})
#ax2.tick_params(axis='both', which='major', labelsize=8)
#plt.title('Price and volume')
#plt.xlabel('Time')
#plt.xlim((min(dt), max(dt)))
#fig.savefig('pandv.png', bbox_inches=0)
#fig.savefig('pandv.pdf', bbox_inches=0) 
   

#date='20081002'
#time='01:01:02.960678'
#off=-1
#time2='03:01:02.960678'
#
#t1=datetime_offset(date,time, off)
#
#print  time_from_datetime(t1)