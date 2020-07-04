import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import matplotlib
import random

def read_data(name):
    
    f=open(name,'r')
    lines=f.readlines()[0].split('\r')
    conc=[]
    Act=[]
    Stdev=[]
    raw_act=[[],[],[]]
    for i in lines:
        line= i.split('\t')
        acti=[float(j) for j in line[1:]]
        raw_act[0].append(float(line[1]))
        raw_act[1].append(float(line[2]))
        raw_act[2].append(float(line[3]))

        conc.append(float(line[0]))
        Act.append(np.mean(acti))
        Stdev.append(np.std(acti)/np.sqrt(len(acti)-1))
        
    return conc, Act, Stdev,raw_act

  
def EC50(params, conc):
   Maxi,EC50,hill=params
   
   Y=1+(Maxi-1)/(1+(conc/EC50)**(-1*hill))
   return Y
   
def EC50_error(params,conc,Act):  
   
   Error=[]
   for c,i in enumerate(conc):
       Error.append(Act[c]-EC50(params,i))
   return Error
   
      
def plot_data(name,conc, Act, Stdev,raw_act,color, line, facecolor):

        plt.errorbar(conc, np.array(Act), yerr=np.array(Stdev), marker='o',linestyle='None',color=color,markeredgecolor=color,markerfacecolor=facecolor)
        plt.xscale('log')
        plt.show()
        ini_guess=[12,0.0000001,1]
        #plt.xlim(0.0000000001, 0.00002)
        fit,success=optimize.leastsq(EC50_error, ini_guess, args=(conc,Act), maxfev=10000000)
        x=np.arange(min(conc),0.00005,min(conc) )
        simul=[]

        print np.log10(fit[1]), fit
        for i in x:
            
            simul.append(EC50(fit, i))
        plt.plot(x,np.array(simul),color=color, linestyle=line)
        plt.xlim(min(conc),0.00005)
        plt.ylim(1,17.5)
       # plt.axvline(fit[1],color='black')
        

if __name__=='__main__':
    plt.figure()
    conc, Act, Stdev,raw_act=read_data('SR1_1ng.txt')
    plot_data('SR1_1ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
    conc, Act, Stdev,raw_act=read_data('SR1bulk_1ng.txt')
    plot_data('SR1bulk_1ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'dashed', 'white')
    conc, Act, Stdev,raw_act=read_data('SR1_L184_1ng.txt')
    plot_data('SR1_L184_1ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')

    plt.savefig('SR1_1ng.EPS')
    
    plt.figure()
    conc, Act, Stdev,raw_act=read_data('SR1_0.25ng.txt')
    plot_data('SR1_0.25ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4' )
    conc, Act, Stdev,raw_act=read_data('SR1bulk_0.25ng.txt')
    plot_data('SR1bulk_0.25ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'dashed', 'white')
    conc, Act, Stdev,raw_act=read_data('SR1_L184_0.25ng.txt')
    plot_data('SR1_4ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')

    plt.savefig('SR1_0.25ng.EPS')
    
    plt.figure()
    conc, Act, Stdev,raw_act=read_data('SR1_4ng.txt')
    plot_data('SR1_4ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
    conc, Act, Stdev,raw_act=read_data('SR1bulk_4ng.txt')
    plot_data('SR1bulk_4ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'dashed', 'white')        
    conc, Act, Stdev,raw_act=read_data('SR1_L184E_4ng.txt')
    plot_data('SR1_4ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'dashed','black')
    plt.savefig('SR1_4ng.EPS')
    plt.close('all')
#    
    
#    conc, Act, Stdev,raw_act=read_data('SR2_4ng.txt')
#    plot_data('SR2_4ng.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
#    
#    conc, Act, Stdev,raw_act=read_data('SR2R161D.txt')
#    plot_data('SR2R161D.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
#    
#    conc, Act, Stdev,raw_act=read_data('SR2D194R.txt')
#    plot_data('SR2D194R.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
#
#    conc, Act, Stdev,raw_act=read_data('SR2F246R.txt')
#    plot_data('SR2F246R.txt',conc, Act, Stdev,raw_act,'#0261B7', 'solid','#0261B7')
#    plt.savefig('SR2_mutants.EPS')
#    
#    
    

    
#    conc, Act, Stdev,raw_act=read_data('SR2R161D.txt')
#    plot_data('SR2R161D.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
#    
#    conc, Act, Stdev,raw_act=read_data('SR2D194R.txt')
#    plot_data('SR2D194R.txt',conc, Act, Stdev,raw_act,'#B15FA4', 'solid','#B15FA4')
#
#    conc, Act, Stdev,raw_act=read_data('SR2F246R.txt')
#    plot_data('SR2F246R.txt',conc, Act, Stdev,raw_act,'#0261B7', 'solid','#0261B7')
#    plt.savefig('SR2_mutants.EPS')

