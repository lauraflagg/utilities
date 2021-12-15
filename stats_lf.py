import numpy as np
import matplotlib.pyplot as plt

import csv

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

def FWHM_to_sigma(fwhm):
    return fwhm/np.sqrt(8*np.log(2))

def sigma_to_FWHM(sigma):
    return sigma*np.sqrt(8*np.log(2))

def parabola(x, bas, qua, cen):
    return bas+qua*(x-cen)**2

def chisq(arr1,arr2):
    x=np.sum((arr1-arr2)**2)
    return x
#chi sq func    

def chisqwithshift(f,g,lb,ub):
    l=ub-lb
    arrl=len(f)
    #a1=np.zeros(arrl+l+1)
    #a2=a1
    shifts=np.arange(lb, ub+1)
    shifts=shifts.astype(np.int)
    
    corrs=[]
    for shift in shifts:
        if shift>0:
            a1=g[0:arrl-shift]
            a2=f[shift:]
        elif shift==0:
            a1=g
            a2=f
        else:
            a1=g[-shift:]
            a2=f[0:arrl+shift]
            
        corrs.append(chisq(a1,a2))
        
    return corrs        

def iseven(n):
    if (n % 2)==0:
        res=True,1
    else:
        res=False,0
    return res

def is_ele_even(arr):
    x=[]
    for item in arr:
        x.append(iseven(item))
    return x

def everyother(arr,maxl=10):
    temp=arr
    while len(temp)>maxl:
        c=0
        temp2=[]
        while c<len(temp):
            temp2.append(temp[c])
            c=c+2
        temp=temp2    
    return temp    

def func2d(func, arr1, arr2, plot=False, println=False, printmins=True,xlab='xlab',ylab='ylab',gtit='gtitle',write='noname'):
    import seaborn as sns
    l1=len(arr1)
    l2=len(arr2)
    finals=np.zeros((l1,l2))
    i=0
    while i<l1:
        j=0
        while j<l2:
            finals[i][j]=func(arr1[i],arr2[j])
            if println==True:  
                print(arr1[i],arr2[j], finals[i][j])
            if write!='noname':
                with open(write, 'wb') as fp: 
                    writer = csv.writer(fp, delimiter=',') # 
                   # writer.writerow(["date", "power_of_telluric_calculatedvia_lmfit"]) # write header 
                    writer.writerows(finals)	        
                                
            j=j+1
        i=i+1
    if plot==True:                
        sns.set()
        lx=arr2
        sns.heatmap(finals, xticklabels=lx, yticklabels=arr1)
        if gtit != 'gtitle':
            plt.title(gtit)
        if xlab != 'xlab':
            plt.xlabel(xlab)
        if ylab != 'ylab':
            plt.ylabel(ylab)            
        plt.show()
        plt.close()
    if printmins==True:
        minloc=np.where(finals==np.min(finals))
        min2=arr2[minloc[1]]
        min1=arr1[minloc[0]]
        print ('minimized at:', min1, min2 )   

    return finals    

def xcor_fast(a,b,lb,ub, dispersion=0):
    '''a, b are the two arrays; 
    lb and ub are the lower and upper bounds of the pixel shift
    dispersion is dispersion in velocity
    This is less accurate than the other versions, but much faster,
    because it only calculates the stdev once per array'''
    a = (a - np.mean(a)) 
    b = (b - np.mean(b)) 
    lags=np.arange(lb,ub+1,1)
    cc= signal.correlate(a, b, 'full')
    ll=len(cc)
    lw=ub-lb+1
    #bord=int((ll-lw)/2)
    cent=int((ll-1)/2)
    corrs= cc[(cent+lb):(cent+ub+1)]/((len(arr1)-np.abs(lags))*(np.std(a)*np.std(b) ))
    if dispersion!=0:
        rvs=lags*dispersion
        corrs=corrs,rvs
    return corrs 

def xcor(f, g, lb, ub, dispersion=0):
    #print(version)
    ''''ub and lb are upper and lower bounds of pixel shift
    dispersion in velocity' '''
    #print(len(f))
    if np.abs(lb)>len(f) or ub>len(f):
        message='Bounds at '+str(lb)+', '+str(ub)+' too wide for array of length '+str(len(f))
        raise ValueError(message)

    l=ub-lb
    arrl=len(f)
    #a1=np.zeros(arrl+l+1)
    #a2=a1
    shifts=np.arange(lb, ub+1)
    shifts=shifts.astype(np.int)
    
    
    corrs=[]
    for shift in shifts:
        if shift>0:
            a1=g[0:arrl-shift]
            a2=f[shift:]
        elif shift==0:
            a1=g
            a2=f
        else:
            a1=g[-shift:]
            a2=f[0:arrl+shift]
        
        a1=a1-np.mean(a1)
        a2=a2-np.mean(a2)
        

        
        s1=np.sqrt(np.sum(a1**2)/len(a1))
        s2=np.sqrt(np.sum(a2**2)/len(a2))
        
        R_fg=np.sum((a1*a2))/len(a1)

        R=R_fg/(s1*s2)

        corrs.append(R)
    corrs=np.array(corrs)
    
    if dispersion!=0:
        rvs=np.arange(lb,ub+1,1)*dispersion
        corrs=corrs,rvs
    return corrs    

def xcor(f, g, lb, ub, dispersion=0):
    #ub and lb are upper and lower bounds of pixel shift
    #print(len(f))
    if np.abs(lb)>len(f) or ub>len(f):
        message='Bounds at '+str(lb)+', '+str(ub)+' too wide for array of length '+str(len(f))
        raise ValueError(message)

    l=ub-lb
    arrl=len(f)
    #a1=np.zeros(arrl+l+1)
    #a2=a1
    shifts=np.arange(lb, ub+1)
    shifts=shifts.astype(np.int)
    
    corrs=[]
    for shift in shifts:
        if shift>0:
            a1=g[0:arrl-shift]
            a2=f[shift:]
        elif shift==0:
            a1=g
            a2=f
        else:
            a1=g[-shift:]
            a2=f[0:arrl+shift]
            
        corrs.append(np.corrcoef(a1,a2)[0,1])
    
    #edit on 20190712
    corrs=np.array(corrs)
    
    if dispersion!=0:
        rvs=np.arange(lb,ub+1,1)*dispersion
        corrs=corrs,rvs
    return corrs    

def binarr(arr,wid=4):
    ll=len(arr)
    i=0
    binned=[]
    while i<ll-wid:
        binned.append(np.mean(arr[i:i+wid]))
        i=i+wid
    arr_out=np.array(binned)
    return arr_out

def coefftofit(x,coeff):
    l=len(coeff)
    #coeff is array output of polyfit
    i=0
    y=x*0
    while i<l:
        y=y+coeff[i]*x**(l-i-1)
        i=i+1
    return y

def BIC(daty,modely,k): #NEED TO TEST
    #k=number of variables
    n=len(daty)
    sse=np.sum((daty-modely)**2)    
    return n*np.ln(sse/n) + k*np.ln(n)


def bootstrap(data, num_samples, statistic, alpha): #need to test
    #","B"Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic.","B" 
    #from http://people.duke.edu/~ccc14/pcfb/analysis.html
    #example of statistic is np.mean or np.median
    #example usage: low, high = bootstrap(x, 100000, np.mean, 0.05)
    n = len(data) 
    idx = np.random.randint(0, n, (num_samples, n)) 
    samples = data[idx] 
    stat = np.sort(statistic(samples, 1)) 
    return (stat[int((alpha/2.0)*num_samples)], stat[int((1-alpha/2.0)*num_samples)])


def gaussian(x, amp, cen, wid, bas):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2)) + bas

def gaussian2D(x, y, amp, xcen, xwid, bas, ycen, ywid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    xs,ys=np.meshgrid(x,y)
    return amp * np.exp(-((xs-xcen)**2/(2*xwid**2))-((ys-ycen)**2 /(2*ywid**2))) + bas



def log_likelihood(f, g, lb, ub, dispersion=0,version='bl',maskvalue=1e6):
    #print(version)
    ''''ub and lb are upper and lower bounds of pixel shift
    version='bl' (default) or 'z' or 'xc' or 'xcor'
    maxvalue will mask out all values of a certain value'''
    #print(len(f))
    if np.abs(lb)>len(f) or ub>len(f):
        message='Bounds at '+str(lb)+', '+str(ub)+' too wide for array of length '+str(len(f))
        raise ValueError(message)

    l=ub-lb
    arrl=len(f)
    #a1=np.zeros(arrl+l+1)
    #a2=a1
    shifts=np.arange(lb, ub+1)
    shifts=shifts.astype(np.int)
    
    
    corrs=[]
    for shift in shifts:
        if shift>0:
            a1=g[0:arrl-shift]
            a2=f[shift:]
        elif shift==0:
            a1=g
            a2=f
        else:
            a1=g[-shift:]
            a2=f[0:arrl+shift]
        
        a1=a1-np.mean(a1)
        a2=a2-np.mean(a2)
        
        okaylocs=np.where((a1!=maskvalue) & (a2!=maskvalue))
        
        s1=np.sqrt(np.sum(a1[okaylocs]**2)/len(a1[okaylocs]))
        s2=np.sqrt(np.sum(a2[okaylocs]**2)/len(a2[okaylocs]))
        
        R_fg=np.sum((a1[okaylocs]*a2[okaylocs]))/len(a1[okaylocs])
        logging.debug('LOG FUNC')
        logging.debug(R_fg)
        logging.debug(len(a1)*np.log(2*np.pi))
        

        if version=='bl':          
            
            R=-len(a1[okaylocs])/2.*np.log(s1**2+s2**2-2*R_fg)-len(a1[okaylocs])*np.log(2*np.pi)
            #R=-len(a1)/2.*np.log(np.log(s1*s2)+np.log(s1/s2+s2/s1-2*(R_fg)/(s1*s2)))-len(a1)*np.log(2*np.pi)
        elif version=='xcor' or version=='xc':
            R=R_fg/(s1*s2)
        else:
            logging.debug((R_fg)**2)
            logging.debug((s1*s2)**2)
            logging.debug(1-(R_fg)**2/(s1*s2)**2)
            logging.debug(np.log(1-(R_fg)**2/(s1*s2)**2))
            R=-len(a1[okaylocs])/2.*np.log(1-(R_fg)**2/(s1*s2)**2)-len(a1[okaylocs])*np.log(2*np.pi)
        logging.debug(R_fg)
        corrs.append(R)
    
    #edit on 20190712
    corrs=np.array(corrs)
    
    if dispersion!=0:
        rvs=np.arange(lb,ub+1,1)*dispersion
        corrs=corrs,rvs
    return corrs  

def log_likelihood_zucker(f, g, lb, ub, dispersion=0,version='z',maskvalue=0):
    #print(version)
    ''''ub and lb are upper and lower bounds of pixel shift
    version='bl' (default) or 'z' or 'xc' or 'xcor'
    maxvalue will mask out all values of a certain value'''
    #print(len(f))
    if np.abs(lb)>len(f) or ub>len(f):
        message='Bounds at '+str(lb)+', '+str(ub)+' too wide for array of length '+str(len(f))
        raise ValueError(message)

    l=ub-lb
    arrl=len(f)
    #a1=np.zeros(arrl+l+1)
    #a2=a1
    shifts=np.arange(lb, ub+1)
    shifts=shifts.astype(np.int)
    
    
    corrs=[]
    for shift in shifts:
        if shift>0:
            a1=g[0:arrl-shift]
            a2=f[shift:]
        elif shift==0:
            a1=g
            a2=f
        else:
            a1=g[-shift:]
            a2=f[0:arrl+shift]
        
        a1=a1-np.mean(a1)
        a2=a2-np.mean(a2)
        
        okaylocs=np.where((a1!=maskvalue) & (a2!=maskvalue))
        
        s1=np.sqrt(np.sum(a1[okaylocs]**2)/len(a1[okaylocs]))
        s2=np.sqrt(np.sum(a2[okaylocs]**2)/len(a2[okaylocs]))
        
        R_fg=np.sum((a1[okaylocs]*a2[okaylocs]))/len(a1[okaylocs])
        logging.debug('LOG FUNC')
        logging.debug(R_fg)
        logging.debug(len(a1)*np.log(2*np.pi))
        

        if version=='bl':          
            
            R=-len(a1[okaylocs])/2.*np.log(s1**2+s2**2-2*R_fg)-len(a1[okaylocs])*np.log(2*np.pi)
            #R=-len(a1)/2.*np.log(np.log(s1*s2)+np.log(s1/s2+s2/s1-2*(R_fg)/(s1*s2)))-len(a1)*np.log(2*np.pi)
        elif version=='xcor' or version=='xc':
            R=R_fg/(s1*s2)
        else:
            logging.debug((R_fg)**2)
            logging.debug((s1*s2)**2)
            logging.debug(1-(R_fg)**2/(s1*s2)**2)
            logging.debug(np.log(1-(R_fg)**2/(s1*s2)**2))
            R=-len(a1[okaylocs])/2.*np.log(1-(R_fg)**2/(s1*s2)**2)-len(a1[okaylocs])*np.log(2*np.pi)
        logging.debug(R_fg)
        corrs.append(R)
    
    #edit on 20190712
    corrs=np.array(corrs)
    
    if dispersion!=0:
        rvs=np.arange(lb,ub+1,1)*dispersion
        corrs=corrs,rvs
    return corrs  