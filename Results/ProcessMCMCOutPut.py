    
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
        
        