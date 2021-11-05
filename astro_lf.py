import numpy as np
import spectres
from scipy.interpolate import interp1d
from astropy import constants, units

c_kms=constants.c.value/10**3

def coadd1d(master_wl,wl_arr,fl_arr):
    '''program resamplespectra to common wavelength scale and then coadds them
    master_wl is a 1d arr, wl_arr is a 2darr, fl_arr is a 2d arr with each row a spec to be coadded'''
    fluxes=[]
    for w,f in zip(wl_arr,fl_arr):
        new_f=spectres.spectres(master_wl,w,f)
        fluxes.append(new_f)
    fluxes=np.array(fluxes)
    
    return np.sum(fluxes,0)

def vel2wl(vel_arr,central_wl):
    '''returns wl 
        inputs are vel_arr and central_wl'''
    wl_arr=(vel_arr*central_wl)/c_kms+central_wl
    return wl_arr
    
    


def wl2vel(wl_arr,central_wl=None):
    '''returns vel in km/s
    
        inputs are wl_arr and central_wl'''
    if central_wl==None:
        central_wl=np.average(wl_arr)

    #print(c_kms)
    vel_arr=(wl_arr-central_wl)*c_kms/central_wl
    return vel_arr

def wn2wl(arr,units='microns'):
    '''assume wave number is in 1/cm'''
    return 10.0**4/np.array(arr)


def deltawltovelo(delta,wl):
    '''input is delta, wl
    returns velocity in km/s'''
    vel=c_kms*delta/wl
    return vel

def veltodeltawl(vel,wl):
    delta=vel*wl/(c_kms)
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

    for j in range(nr):
        r=dr/2.+j*dr
        area=((r+dr/2.)**2-(r-dr/2.)**2)/int(ntheta*r)*(1.-eps+eps*np.cos(np.arcsin(r)))
        for k in range(int(ntheta*r)-1):
            th=np.pi/int(ntheta*r)+k*2.*np.pi/int(ntheta*r)
            if dif!=0:
                vl=vsini*r*np.sin(th)*(1.-dif/2.-dif/2.*np.cos(2.*np.arccos(r*np.cos(th))))
                interped = interp1d(wl+wl*vl/(3*10**5), flux,fill_value='extrapolate')
                ns=ns+area*interped(wl) #interpol(flux,wl+wl*vl/3.e5,wl) 
                tarea=tarea+area
            else:
                vl=r*vsini*np.sin(th)
                interped = interp1d(wl+wl*vl/(3*10**5), flux,fill_value='extrapolate')
                ns=ns+area*interped(wl)#     ns=ns+area*interpol(flux,wl+wl*vl/3.e5,wl) 
                tarea=tarea+area

    ns=ns/tarea
    return ns


def findbests2n(spec,width=20,edge=10,p=100):
    '''inputs: spec,width, edge,p
    p is percentile from 0 to 100
    for best use 100 or 99 or somethinglike that
    width is in pixels
    edge in pixels; for a non-blaze corrected, the best s2n will never be at ends'''
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

def findcontlevel(spec,width=20,edge=10,p=100):
    '''inputs: spec,width, edge,p
    p is percentile from 0 to 100
    for best use 100 or 99 or somethinglike that
    width is in pixels
    edge in pixels; for a non-blaze corrected, the best s2n will never be at ends'''
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
        #s2n=np.mean(arr)/np.std(arr,ddof=1)
        #s2nstemp.append(s2n)
        means.append(avg)
        stds.append(std)
        i=i+1  
    
    stds=np.array(stds)
    stds[np.where(stds==0)]=1e5
    lowest_sd_loc0=np.where(stds==np.min(stds))
    #print(lowest_sd_loc0)
    lowest_sd_loc00=lowest_sd_loc0[0]
    lowest_sd_loc=lowest_sd_loc00[0]
    #print(lowest_sd_loc0,lowest_sd_loc00,lowest_sd_loc)

    cont=np.average(spec[lowest_sd_loc+edge:lowest_sd_loc+width+edge])
    return cont

def equilibtemp(T_s,dist,albedo,R_star,e=1,slowrotator=True):
    '''inputs: T_s, dist, albedo,R_star,e,slowrotator
    e is emmisivity
    dist in AU
    R_star in solar radii'''

    r=R_star*6.957*10**5
    d=dist*149597871


    if slowrotator==True:
        factor=2
    else:
        factor=4
    T=T_s*((1-albedo)/factor*e)**0.25*(r/d)**0.5

    return T



def createwlscale(delt,low=1000,up=1450):
    '''creates a wavelengths cale that's evenly spaced in velocity space --not wavelength -- for easy cross-correlating
    low and up in wavelength
    delt in km/s'''
    wlc=low
    wl=[wlc]
    while wlc<=up:
        #wlc is wl current
        delta=(delt/(c_kms))*wlc
        wl.append(wlc+delta)
        wlc=wlc+delta
    wl=np.array(wl)

    return wl

def getorbitpars(period=None,m_s=None,sep=None,v_p=None,m_p=None,e=0.,i=90):
    '''m_s in solar masses
    m_p in jupiter masses
    period in days
    sep in au
    v_p in km/s'''
    sec_in_day=60*60*24
    
    if m_p==None:
        m_p=m_s*1e-5
    
    if v_p==None:
        if period==None:
            period=np.sqrt(4*np.pi**2*(sep*constants.au.value)**3/(constants.G.value*(m_s*constants.M_sun.value+m_p*constants.M_jup.value)))/sec_in_day
        if sep==None:
            sep=((period*sec_in_day)**2*(constants.G.value*(m_s*constants.M_sun.value+m_p*constants.M_jup.value))/(np.pi**2*4))**.333333333333/constants.au.value
        v_p=np.pi*2*sep*constants.au.value/(period*sec_in_day)*1e-3
        k_s=203*(period)**(-.33333333333)*m_p*np.sin(2*np.pi/360*i)/(m_s+9.548e-4*m_p)**.666666666*1./np.sqrt(1-e**2)
        k_p=v_p*np.sin(2*np.pi*i/360)
        
    return {'v_p':v_p*units.km/units.s,'sep':sep*units.AU,'period':period*units.day,'k_s':k_s*units.m/units.s,'k_p':k_p*units.km/units.s}

def p_to_s_massratio(m_p,m_s):
    '''    m_s, in solar masses
    m_p, in jupiter masses'''
    return m_s*constants.M_sun/(m_p*constants.M_jup)

def calc_asini(m_p,m_s,a_p):
    '''calculates the asini of a star from the orbital parameters of a planet
    m_s, in solar masses
    m_p, in jupiter masses
    a_p in au
    return a_s*sini in Gm'''
    s2p_massratio=p_to_s_massratio(m_p=m_p,m_s=m_s)
    a_km=constants.au.value*1e-3*a_p
    a_gm=a_km*1e-6
    
    return a_gm/s2p_massratio


class wlunit():
    def __init__(self,names,conversion):
        '''names is a list of strings that can be used to identify the wl unit
        conversion is the factor compared to m; i.e. cm would be 1e-2'''
        self.names=names
        self.conversion=conversion
        self.identifier=names[0]
        
wl_unit_choices={'aa':wlunit(['aa','AA','angstroms','Angstroms'],1e-10),'microns':wlunit(['microns','Microns','mum'],1e-6),'nm':wlunit(['nm'],1e-9)}
             
    

def guess_unit(wl):
    if wl<100:
        return 'microns'
    elif wl<900:
        return 'nm'
    else:
        return 'aa'

def find_unit(wl_unit):
    for item in wl_unit_choices:
        temp=wl_unit_choices[item]
        if wl_unit in temp.names:
            return item
    
    

def v_curve(g,pd,ai,e,w):
    #pro vcurve,g,pd,ai,e,w,fout
    
    # This procedure creates a 2 x 360 array containing the 
    # radial velocity data as a function of time for a spectroscopic
    # binary.  Input g (center of mass velocity of system), pd
    # (period of orbit in days), ai (asin(i) of orbit in giga-meters),
    # e (eccentricity of orbit), w (longitude of periastron), and
    # fout (name of output file).
    # Formulae from Heintz (1978), Danby (1990), and math tables.
    # LAP 9.1.96
    
    # print statement to prompt input:
    
    #print,'type g(km/s),pd(days),ai(Gm),e(ecc), w(deg) & output (filename)'
    
    #port from IDL by Laura Flagg    
    fout=np.zeros((2,360))
    out=np.zeros((2,360))
   
   # define pi and convert to useful units
   
    pi=np.pi
   # convert to seconds:
    p=pd*8.64e4
   # convert to km:
    asi=ai*1e6
   # convert to radians:
    wr=(w/360.0)*2*pi
   
   #stop
   
    t1=g
    t2=2*pi/p #units of s^-1, angular frequency
    t3=asi #distaance
    ed=(1-e**2)
    t4=np.sqrt(1/ed)
    t5=e*np.cos(wr) #large if wr is close to 0, dimensionless
    erat=(np.sqrt((1-e)/(1+e)))
   
   # print,sqrt((1-e)/(1+e))
   #stop
    for ji in range(0,360):
   
        j=ji+540
        jr=(j/360.0)*2*pi #the angle in radians
   
        t6=np.cos(jr+wr) #combine the positiona angle with the angle of periastron
   
        v=t1+(t2*t3*t4*(t5+t6))
   
        prt1=(e*np.sin(jr)*np.sqrt(ed)/(1+e*np.cos(jr)))
        prt2=2*np.arctan(erat*np.tan(jr/2.0))
        if ji == 0:
            prt2=-pi      
         #because of weird idl, they get -pi, while python gets pi, so correcting that
        jj=j-540
        t =(1/t2)*(prt2-prt1)
    
        out[0,jj]=t/p
        out[1,jj]=v
    #  if ji == 180:
     #    print prt1, prt2, t, t2, t6,p, out[0,jj], jj, jr, np.tan(jr/2.0)

         
   
   #   print,j,jr,t6,v,out(0,jj),out(1,jj)
   
      
    if out[0,0] < 0.0:
        out[0,]=out[0,] - out[0,0]
      
   
    fout[0,0:180]=out[0,180:]-out[0,180]
    fout[0,180:]=out[0,0:180]+out[0,180]
    fout[1,0:180]=out[1,180:]
    fout[1,180:]=out[1,0:180]
   
    return fout


def do_vcurve(phase,par):#
    #function func_citau,phase,par
    #
    #  This function computes the radial velocity amplitude for a single
    #  line spectroscopic binary, using the code vcurve.pro supplied by
    #  Lisa Prato.  
    #  INPUTS:
    #    phase - The orbital phase for the desired points
    #    par(0) - g: center of mass velocity of system in m/s
    #    par(1) - pd: period of orbit in days
    #    par(2) - ai: asin(i) of orbit in giga-meters
    #    par(3) - e: eccentricity of orbit
    #    par(4) - w: longitude of periastron
    #    par(5) - ph0: Phase offset    
#  OUTPUTS:
#    function returns the velocity of the star in m/s
#
#  HISTORY:
#    10-Apr-2014 CMJ - Written, based on func420.pro for XO-3b
#    25-Feb-2007 CMJ - Written
#    13-Mar-2007 CMJ - Added Phase offset term
#

    # Set up variables
    g = par[0]/1000.
    pd = par[1]
    par[2] = abs(par[2])
    ai = par[2]
    par[3] = abs(par[3])
    e = par[3]
    par[4] = par[4] % 360.
    w = par[4]
    par[5] = (par[5]+20.) % 1.
    ph0 = par[5]
    #if ph0 lt 0. then ph0 = 0.
    #if ph0 gt 1. then ph0 = 1.

    fout=v_curve(g,pd,ai,e,w)

    # Interpolate onto phases and return
    #
    #vel = 1.d3*interpol(reform(fout(1,*)),reform(fout(0,*)),((phase+ph0) mod 1.)) 
    a=fout[1]
    b=fout[0]
    c=((phase+ph0) % 1.)
    vel=np.interp(c,b,a)*1000.
    #in m/s
    #vel(10:20) = vel(10:20) + par(6)             # adjust HET velocities



    return vel

def coadd_echelle(w_matrix,f_matrix,e_matrix,wls,ledgecut=5,redgecut=5):
    flux_list=[]
    unc_list=[]
    spec_count=np.zeros_like(wls)
    flux=np.zeros_like(wls)
    unc=np.zeros_like(wls)

    for w,f,e in zip(w_matrix,f_matrix,e_matrix):

        interped_fl,interped_unc=spectres.spectres(wls,w[ledgecut:-(redgecut+1)],f[ledgecut:-(redgecut+1)],spec_errs=e[ledgecut:-(redgecut+1)],fill=0,verbose=False)

        interped_region=np.where(interped_fl!=0)
        spec_count[interped_region]+=1
        flux_list.append(interped_fl)
        unc_to_use=np.nan_to_num(interped_unc)
        unc_list.append(unc_to_use)
    for i,item in enumerate(wls):


        fs=np.array([cnt[i] for cnt in flux_list])
        us=np.array([cnt[i] for cnt in unc_list])
        good=np.where((us!=0) & (fs!=0))
        ws=1./us[good]**2



        if len(ws)>0:
            w_sum=np.sum(ws)
            f_avg=np.average(fs[good],weights=ws**2)
            flux[i]=f_avg

            #if i==29000 or i==28400:
                #   print('i,fs,us,good,ws,f_avg',i,fs,us,good,ws,f_avg)

            if len(ws)==1:
                unc[i]=us[good]
            else:
                #self.unc[i]=np.sqrt(np.sum(us[good]**2*ws/np.sum(ws))/(len(ws)-1))
                unc[i]=np.sqrt(1./np.sum(1./us[good]**2))

    return flux,unc