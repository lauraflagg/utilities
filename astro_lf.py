import numpy as np
#import spectres
from scipy.interpolate import interp1d

def wn2wl(arr,units='microns'):
    #assume wave number is in 1/cm
    #arr should be numpy arr
    return 10.0**4/arr


def readirafwl(head):
    #input like:
    #	a0v=fits.open(files[i])
    #	heads=a0v[0].header
    for item in head:
        x=2

def deltawltovelo(delta,wl):
    vel=3.*10.**5*delta/wl
    return vel

def veltodeltawl(vel,wl):
    delta=vel*wl/(3*10**5)
    return delta

def rotbroad(wl, flux, vsini, eps=0.6,nr=10,ntheta=100, dif=0):

    # based on CMJ's program rotint.pro edited May 11 1994
    # ;  This routine reads in a spectrum, s, on a wavelength scale, w, and a vsini
    # ;  with which to rotationally broaden the spectrum.  The rotationally broadened
    # ;  spectrum is returned in ns.  Parameters that can be set are the coefficient
    # ;  of the limb darkening law, eps (0.6 default), the number of radial steps on
    # ;  the disk, nr (default = 10), and the maximum number of steps in angle around
    # ;  the disk, ntheta (default = 100).  Final optional parameter dif allows for 
    # ;  differential rotation according to the law Omeg(th)/Omeg(eq) = (1. - dif/2
    # ;  - (dif/2) cos(2 th)).  Dif = .675 nicely reproduces the law proposed by
    # ;  Smith, 1994, A&A, in press. to unify WTTS and CTTS.  Dif = .23 is similar to
    # ;  observed solar differential rotation.  Note: the th in the above expression
    # ;  is the stellar co-latitude, not the same as the integration variable used
    # ;  below.  This is a disk integration routine.

    #note from LF: using interpolate instead of spline



    final_spec=np.zeros(len(flux))

    ns=np.zeros(len(flux))
    tarea=0.0

    dr=1./nr
    j=0
    while j<nr:
        r=dr/2.+j*dr
        area=((r+dr/2.)**2-(r-dr/2.)**2)/int(ntheta*r)*(1.-eps+eps*np.cos(np.arcsin(r)))
        k=0
        while k<int(ntheta*r)-1:
            th=np.pi/int(ntheta*r)+k*2.*np.pi/int(ntheta*r)
            if dif!=0:
                vl=vsini*r*np.sin(th)*(1.-dif/2.-dif/2.*np.cos(2.*np.arccos(r*np.cos(th))))
                interped = interp1d(wl+wl*vl/(3*10**5), flux)
                ns=ns+area*interped(wl) #interpol(flux,wl+wl*vl/3.e5,wl) 
                tarea=tarea+area
            else:
                vl=r*vsini*np.sin(th)
                interped = interp1d(wl+wl*vl/(3*10**5), flux,fill_value='extrapolate')
                ns=ns+area*interped(wl)#     ns=ns+area*interpol(flux,wl+wl*vl/3.e5,wl) 
                tarea=tarea+area
            k=k+1
        j=j+1



    ns=ns/tarea
    return ns


def findbests2n(spec,width=20,edge=10,p=100):
    #p is percentile from 0 to 100
    #for best use 100 or 99 or somethinglike that
    #width is in pixels
    #edge in pixels; for a non-blaze corrected, the best s2n will never be at ends
    l=len(spec)

    s2nstemp=[]
    means=[]
    stds=[]
    i=edge

    while i<(l-width-edge):
        rangeend=i+width
        arr=spec[i:rangeend]
        avg=np.mean(arr)
        std=np.std(arr,ddof=1)
        s2n=np.mean(arr)/np.std(arr,ddof=1)
        s2nstemp.append(s2n)
        means.append(avg)
        stds.append(std)
        i=i+1  

    s2nstemp=np.nan_to_num(s2nstemp)
    return np.percentile(s2nstemp,p)

def equilibtemp(T_s,dist,albedo,R_star,e=1,slowrotator=True):
    #e is emmisivity
    #dist in AU
    #R_star in solar radii

    r=R_star*6.957*10**5
    d=dist*149597871


    if slowrotator==True:
        factor=2
    else:
        factor=4
    T=T_s*((1-albedo)/factor*e)**0.25*(r/d)**0.5

    return T