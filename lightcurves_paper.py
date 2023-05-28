import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

#############################################################################
KP=3.08568025*pow(10.,19); #// in meter.
G=6.67384*pow(10.,-11.0);#// in [m^3/s^2*kg].
velocity= 3.0*pow(10.0,8.0);#//velosity of light
Msun=1.98892*pow(10.,30);# //in [kg].
Rsun= 6.957*pow(10.0, 8.0); #///solar radius [meter]
AU=1.4960*pow(10.0,11.0);
year=365.2421875;
###############################################################################
plt.clf()


fif=open("./residualsA.txt","w")
fif.close()


nr=int(49)
f0=open("./files/IBH_ROMAN.txt","r")
nn=sum(1 for line in f0)
par=np.zeros(( nn , nr ))
par=np.loadtxt("./files/IBH_ROMAN.txt") 

nam=["icon","b", "l", #3
   "l.struc", "l.mass", "l.Dl", "l.vl", #7
    "s.struc", "s.cl", "s.mass", "s.Ds", "s.Tstar", "s.Rstar", "s.logl", "s.type", "s.col", "s.vs", "s.Mab[1]", 
    "s.Mab[4]", "s.Map[1]", "s.Map[4]", #21
    "s.magb[1]", "s.magb[4]", "s.blend[1]", "s.blend[4]","s.Nblend[1]", "s.Nblend[4]", "Ext[1]", "Ext[4]", #29
    "l.tE", "l.RE(AU)", "l.t0", "l.mul", "l.Vt", "l.u0", "l.ksi", "s.opt(e6)", "rostar", "tetE" , #39
    "flagf","flagdet", "dchi", "ndw", "direction", "mus1", "mus2", "s.tetv","mul1", "mul2", "piE"] #46
'''
for i in range(nr):
    plt.figure()
    #if(i==60 or i==50): plt.hist(np.log10(par[:,i]),35,color="g", histtype="bar")
    plt.hist(par[:,i],35,color="g", histtype="bar")
    plt.title(str(nam[i]))
    plt.grid(True)
    plt.savefig("./histo/Histo{0:d}.jpg".format(i), dpi=200)
    print "********** HISTO IS PLOTTED *************",  i
'''
#################################################################
for i in range(nn):
    
    icon, lat, lon=       int(par[i,0]), par[i,1], par[i,2]
    strucl, Ml, Dl, vl=    par[i,3], par[i,4], par[i,5], par[i,6]
    strucs, cl, mass, Ds,Tstar, Rstar, logl= par[i,7], par[i,8], par[i,9], par[i,10], par[i,11], par[i,12], par[i,13]
    types, col, vs, MI, MW149, mI, mW149=    par[i,14], par[i,15],par[i,16], par[i,17], par[i,18], par[i,19], par[i,20]
    magbI, mbs, blendI, fb, NbI, Nbw, ExI, ExW= par[i,21], par[i,22],par[i,23], par[i,24], par[i,25], par[i,26], par[i,27], par[i,28]
    tE, RE, t0, mul, Vt, u0, opd,ros,tetE=par[i,29],par[i,30],par[i,31], par[i,32], par[i,33], par[i,34], par[i,35],par[i,36], par[i,37]
    flagf,flagdet, dchi, ndw,li, mus1, mus2=int(par[i,38]), par[i,39],par[i,40], par[i,41], par[i,42], par[i,43], par[i,44]
    xi, mul1, mul2, piE=  par[i,45], par[i,46],par[i,47], par[i,48]
    pirel= float(1.0/Dl-1.0/Ds)## marcs
    print("Parameters: i,  icon, tE, tetE, piE :   ", i,  icon,  tE,    tetE,    piE)
    if(flagf>0): 
        f1=open("./files/dat_{0:d}.dat".format(icon),"r")
        nd= sum(1 for line in f1)  
        dat=np.zeros((nd,7)); 
        dat= np.loadtxt("./files/dat_{0:d}.dat".format(icon))
        
        
        f2=open("./files/mag_{0:d}.dat".format(icon),"r")
        nm=sum(1 for line in f2)  
        mag1=np.zeros((nm,22));  mag=np.zeros((nm,22))  
        mag1= np.loadtxt("./files/mag_{0:d}.dat".format(icon))
        nz=0; 
        for k in range(nm): 
            if(mag1[k,21]==0): nz+=1
        k=0; 
        for j in range(nm):
            if(mag1[j,21]==0): 
                mag[k,:]=mag1[nz-j-1,:]
            else:   
                mag[k,:]=mag1[j,:]
            k+=1    
        
        dat1=np.zeros((nd,7));   d1=0         
        dat2=np.zeros((nd,7));   d2=0
        for j in range(nd): 
            if(dat[j,6]>0):
                dat1[d1,:]= dat[j,:]
                d1+=1 
            if(dat[j,6]<0):
                dat2[d2,:]= dat[j,:]
                d2+=1     
        errm= np.mean(dat[:,2])
        erra= np.mean(dat[:,5])        
        ###############################################################33    
        '''
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        #plt.errorbar(dat1[:d1,0]/year,dat1[:d1,1],yerr=dat1[:d1,2], fmt=".", markersize=1.4,  color='red', ecolor='gray', elinewidth=0.3, capsize=0)
        #plt.errorbar(dat2[:d2,0]/year,dat2[:d2,1],yerr=dat2[:d2,2], fmt=".", markersize=1.8,  color='red', ecolor='gray', elinewidth=0.3, capsize=0)
        plt.plot(mag[:,0]/year, mag[:,2], "k:", lw=1.2, label=r"$\rm{Magnification}$")
        plt.plot(mag[:,0]/year, mag[:,1], "k--",lw=1.2, label=r"$\rm{Magnification}+\rm{parallax}$")
        plt.xlabel(r"$\rm{time(yrs)}$", fontsize=18)
        plt.ylabel(r"$\rm{W149}-\rm{magnitude}$", fontsize=18)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.title(r"$t_{\rm E}\rm{(days)}=$"+str(round(tE,1))+ r"$,~~\theta_{\rm E}\rm{(mas)}=$"+str(round(tetE,2))+r"$,~~\pi_{\rm E}=$"+str(round(piE,3)),fontsize=18, color="k")
        plt.gca().invert_yaxis()
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        #plt.legend()
        #ax1.legend(prop={"size":16.5})
        ax1.legend(prop={"size":14.5})
        fig=plt.gcf()
        fig.savefig("./lightc/li_{0:d}.jpg".format(icon),dpi=200)
        ###############################################################33    
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        #plt.errorbar(dat[:,3],dat[:,4],yerr=dat[:,5],xerr=dat[:,5], fmt=".", markersize=1.2, color='red', ecolor='gray', elinewidth=0.3, capsize=0)
        #plt.errorbar(dat1[:d1,3],dat1[:d1,4],yerr=dat1[:d1,5],xerr=dat1[:d1,5], fmt=".", markersize=1.2, color='red', ecolor='gray', elinewidth=0.3, capsize=0)
        #plt.errorbar(dat2[:d2,3],dat2[:d2,4],yerr=dat2[:d2,5],xerr=dat2[:d2,5], fmt=".", markersize=1.5, color='red', ecolor='gray', elinewidth=0.3, capsize=0)
        plt.plot(mag[:,9], mag[:,10], "b:", lw=1.5)
        plt.plot(mag[:,5], mag[:,6], "b--", lw=1.5, label=r"$\rm{source(undeflected)}+\rm{parallax}$")
        plt.plot(mag[:,5]+mag[:,13], mag[:,6]+mag[:,14], "b-.", lw=1.5, label=r"$\rm{source(deflected)}+\rm{parallax}$")
        plt.plot(mag[:,5]-mag[:,7], mag[:,6]-mag[:,8], "--", color='darkred', lw=1.5, label=r"$\rm{lens-source(undeflected)}+\rm{parallax}$")
        plt.plot(mag[:,11], mag[:,12], "m:",lw=1.5)
        plt.plot(mag[:,7],  mag[:,8],  "m--",lw=1.5, label=r"$\rm{Lens} + \rm{parallax}$")
        plt.plot(0.0,0.0,'ko',markersize=5.0)
        plt.plot(mag[:,15]+mag[:,17], mag[:,16]+mag[:,18],"-.", color='k',lw=1.8, label=r"$\rm{Deflection}$")
        plt.xlim([-1.5*tetE, 1.5*tetE])
        plt.ylim([-1.5*tetE, 1.5*tetE])
        plt.xlabel(r"$\rm{x~position(mas)}$", fontsize=18)
        plt.ylabel(r"$\rm{y~position(mas)}$", fontsize=18)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.title(r"$u_{0}=$"+str(round(u0,2))+ r"$,~~m_{\rm{base}}\rm{(mag)}=$"+str(round(mbs,2))+r"$,~~t_{0}(\rm{years})=$"+str(round(t0/year,1)),fontsize=18, color="k")
        plt.gca().invert_yaxis()
        #plt.legend()
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        ax1.legend(prop={"size":12.5})
        fig=plt.gcf()
        fig.savefig("./lightc/astro_{0:d}.jpg".format(icon),dpi=200)
        print( "**** Lightcurve was plotted, No:  ",  icon  )
        '''
        ###############################################################33   
        ''' 
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        
        #plt.errorbar(dat[:,3],dat[:,4],yerr=dat[:,5],xerr=dat[:,5], fmt=".", markersize=0.5, color='red', ecolor='gray', elinewidth=0.1, capsize=0)
        plt.plot(mag[:,13], mag[:,14],"-", color='r',lw=1.5, label=r"$\rm{Astrometric}~\rm{Deflection}$")
        plt.plot(mag[:,15], mag[:,16],"--", color='g',lw=1.5, label=r"$\rm{Deflection}+\rm{Parallax}$")
        plt.plot(mag[:,17], mag[:,18],":", color='b',lw=1.5)##, label=r"$\rm{Deflection}+\rm{Parallax}$")
        #plt.xlim([ -0.8*tetE, 0.8*tetE ])
        #plt.ylim([ -0.8*tetE, 0.8*tetE ])
        plt.xlabel(r"$\rm{x~position(mas)}$", fontsize=17)
        plt.ylabel(r"$\rm{y~position(mas)}$", fontsize=17)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.title(r"$M_{\rm{l}}(M_{\odot})=$"+str(round(Ml,1))+ r"$,~~D_{\rm{l}}\rm{(kpc)}=$"+str(round(Dl,1))+r"$,~~\pi_{\rm{rel}}=$"+str(round(pirel,3))+r"$,~~\theta_{\rm E}=$"+str(round(tetE,3))+r"$,~~\pi_{\rm E}=$"+str(round(piE,3)),fontsize=14, color="k")
        plt.gca().invert_yaxis()
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        ax1.legend(prop={"size":12.5})
        fig=plt.gcf()
        fig.savefig("./lightc/Shift_{0:d}.jpg".format(icon),dpi=200)
        print( "**** Lightcurve was plotted, No:  ",  icon  )
        '''
        ###############################################################
        '''
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.plot(mag[:,13], mag[:,14],"-", color='b',lw=1.5, label=r"$\rm{Astrometric}~\rm{Deflection}$")
        plt.plot(mag[:,15], mag[:,16],"--", color='r',lw=1.5, label=r"$\rm{Deflection}+\rm{Parallax}$")
        #plt.xlim([ -0.8*tetE, 0.8*tetE ])
        #plt.ylim([ -0.8*tetE, 0.8*tetE ])
        plt.xlabel(r"$\rm{x~position(mas)}$", fontsize=17)
        plt.ylabel(r"$\rm{y~position(mas)}$", fontsize=17)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.title(r"$M_{\rm{l}}(M_{\odot})=$"+str(round(Ml,1))+ r"$,~~D_{\rm{l}}\rm{(kpc)}=$"+str(round(Dl,1))+r"$,~~\pi_{\rm{rel}}=$"+str(round(pirel,3))+r"$,~~\pi_{\rm E}=$"+str(round(piE,3)),fontsize=16, color="k")
        plt.gca().invert_yaxis()
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        ax1.legend(prop={"size":12.5})
        fig=plt.gcf()
        fig.savefig("./lightc/BShift_{0:d}.jpg".format(icon),dpi=200)
        print( "**** Lightcurve was plotted, No:  ",  icon  )
        '''
        
        ###############################################################33    
    
        #res0= np.mean(abs((mag[:,1]-mag[:,2])/mag[:,1]) ) 
        #res1= np.mean( np.sqrt((mag[:,13]-mag[:,15])**2 + (mag[:,14]-mag[:,16])**2)/np.sqrt(mag[:,13]**2+ mag[:,14]**2) ) 
        amp1= np.max(abs(mag[:,1]-mag[:,2]))## mag 
        amp2= np.max(np.sqrt((mag[:,13]-mag[:,15])**2+(mag[:,14]-mag[:,16])**2) ) ##mas 

        result=np.array([Dl, Ds, Ml, pirel, tetE, piE, tE, errm, erra, amp1 , amp2 ])
        fif=open("./residualsA.txt","a+")
        np.savetxt(fif, result.reshape((1,11)),fmt="%.6f   %.6f   %.6f   %.8f   %.8f  %.8f   %.5f %.8f  %.8f  %.8f  %.8f") 
        #fif.write("\n***************************************\n")
        fif.close()
        







