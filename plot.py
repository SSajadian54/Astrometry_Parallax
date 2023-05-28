import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
####################################################################
##Dl, Ds, Ml, pirel, tetE, piE, tE, errm, erra, amp1 , amp2 
##0    1   2   3       4    5    6   7     8     9     10
f0=open("./residuals.txt","r")
nn=sum(1 for line in f0)
par=np.zeros(( nn , 11))
par=np.loadtxt("./residuals.txt") 
par[:,0]=par[:,0]/par[:,1]
par[:,2]= np.log10(par[:,2])
par[:,3]= np.log10(par[:,3])
par[:,4]= np.log10(par[:,4])
par[:,5]= np.log10(par[:,5])
par[:,6]= np.log10(par[:,6])

nam=[r"$x_{\rm{ls}}$", r"$D_{\rm s}(\rm{kpc})$", r"$\log_{10}[M_{\rm l}(M_{\odot})]$", r"$\log_{10}[\pi_{\rm{rel}}(\rm{mas})]$", r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$", r"$\log_{10}[\pi_{\rm E}]$", r"$\log_{10}[t_{\rm E}(\rm{days})]$"]

nam2=[r"$\log_{10}[\Delta_{\rm{max}} m/\sigma_{\rm m}]$",r"$\log_{10}[\Delta_{\rm{max}} \delta \theta_{\rm c}/\sigma_{\rm{a}}]$"]

x1=[0.01,2.0,0.5, -3.3,-0.5,-3.5, 1.3]
x2=[0.98,16.0,2.71, 0.5,2.2, -1.0, 3.5]


par[:,8]=50.0*0.001

###################################################################
for i in range(7): 
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    plt.plot(par[:,i], np.log10(par[:,9]/par[:,7]),"bo", markersize=4.0, label=str(nam2[0]))
    plt.plot(par[:,i], np.log10(par[:,10]/par[:,8]),"ro", markersize=4.0, label=str(nam2[1]))
    plt.xlabel(str(nam[i]), fontsize=18)
    plt.ylabel(r"$\rm{Relative}~\rm{Residual}$", fontsize=18)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.legend()
    plt.legend(prop={"size":16.5})
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./plot_{0:d}.jpg".format(i),dpi=200)
    
j =int(0)    
dat=np.zeros((nn,11))    
for i in range(nn): 
    if(float(par[i,10]/par[i,8]) > float(par[i,9]/par[i,7])): 
        dat[j,:]=par[i,:]
        j+=1

print ("fraction:  ",  j*100.0/nn  )        
        
for i in range(7): 
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    plt.plot(par[:,i], par[:,10]/par[:,8]/(par[:,9]/par[:,7]) , "bo", markersize=4.5)##, label=str(nam2[0]))
    plt.plot(dat[:j,i],dat[:j,10]/dat[:j,8]/(dat[:j,9]/dat[:j,7]) , "ro", markersize=4.5)##, label=str(nam2[0]))
    
    plt.xlabel(str(nam[i]), fontsize=18)
    plt.ylabel(r"$\pi_{\rm{a}, n}/\pi_{\rm{m}, n}$", fontsize=18)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.yscale('log')
    plt.xlim(float(x1[i]), float(x2[i]) )
    plt.legend()
    plt.legend(prop={"size":16.5})
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./plotC_{0:d}.jpg".format(i),dpi=200)
    





