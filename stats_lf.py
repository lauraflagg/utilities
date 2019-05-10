import numpy as np
import matplotlib.pyplot as plt

import csv

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
        res=True
    else:
        res=False
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
                print( arr1[i],arr2[j], finals[i][j])
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
        print( 'minimized at:', min1, min2 )   

    return finals    

def xcor(f, g, lb, ub, dispersion=0):
    #ub and lb are upper and lower bounds of pixel shift
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
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic.""" 
    #from http://people.duke.edu/~ccc14/pcfb/analysis.html
    #example of statistic is np.mean or np.median
    #example usage: low, high = bootstrap(x, 100000, np.mean, 0.05)
    n = len(data) 
    idx = np.random.randint(0, n, (num_samples, n)) 
    samples = data[idx] 
    stat = np.sort(statistic(samples, 1)) 
    return (stat[int((alpha/2.0)*num_samples)], stat[int((1-alpha/2.0)*num_samples)])
    