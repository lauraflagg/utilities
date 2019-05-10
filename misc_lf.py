import numpy as np

def appendarrs(first, *args):
    if type(first) is 'list':
        for item in args:
            first=first+item
    else:        
        for item in args:
            first=np.concatenate((first, item))
    return first


def removebad(arr_all,arr_bad):
    #input will be numpy arrays
    goodlocs=[]
    i=0
    while i<len(arr_all):
        if (arr_all[i] not in arr_bad)==True:
            goodlocs.append(i)
        i=i+1
    arr_good=arr_all[goodlocs]
    return arr_good, goodlocs