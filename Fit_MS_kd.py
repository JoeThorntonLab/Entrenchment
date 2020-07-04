import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy import optimize


def load_data(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    mz=[]
    inten=[]
    monomer=0
    dimer=0
    for i in lines:
        mz.append(float(i.split()[0]))
        inten.append(float(i.split()[1]))
        if mz[-1]>2000 and mz[-1]<3200:
            monomer+=inten[-1]
        if mz[-1]>3200 and mz[-1]<4200:
            dimer+=inten[-1]
    print dimer*2/(monomer+dimer*2), filename
    return mz,inten
    
def load_fit(filename):
    print filename
    f=open(filename,'r')
    lines=f.readlines()
    if len(lines)==2:
        mono=float(lines[0].split()[1])
        dimer=float(lines[1].split()[1])
        return dimer*2/(mono+dimer*2)
    else:
        return 0
    
def plot_spectrum(mz,inten,filename):
    inten=np.array(inten)/max(inten)#normalize
    plt.figure()
    plt.plot(mz,inten)
    plt.xlim(2000,5000)
    plt.savefig(filename.split('.')[0]+'.eps')
    
def KD_dimer ( M0, Kd):
    '''Function to obtain a KD for a monomer/dmer equilibrium if the amount of dimer can be measured and the absolute concentration of monomer is known'''
    D= (4.*M0+Kd-np.sqrt(8.*M0*Kd+Kd**2))/8.
    return D

def KD_dimer_err (Kd, M0, D ):
    '''Error function fo KD function above'''
    err= KD_dimer(M0, Kd)/M0- D/M0
    return err**2

def fit_kds(files):

    SR1_dimer=[]
    SR1_mol=[]
    bulk_dimer=[]
    bulk_mol=[]
    SR1184E_dimer=[]
    SR1184E_mol=[]
    AF_mol=[]
    AF_dimer=[]
    for i in files:
      if i.split('/')[0]=='Georg':
        if i.split('/')[1].split('_')[1] not in ['184E','BULK']:
            SR1_dimer.append(load_fit(i+'/'+i.split('/')[1][:-11]+'peaks.dat'))
            conc=i.split('/')[1].split('_')[1]
            if conc[-2:]=='NM':
                SR1_mol.append(float(conc[:-2])/1000)
            else:
                SR1_mol.append(float(conc[:-2]))
        elif i.split('/')[1].split('_')[1]=='184E':
            SR1184E_dimer.append(load_fit(i+'/'+i.split('/')[1][:-11]+'peaks.dat'))
            conc=i.split('/')[1].split('_')[2]
            if conc[-2:]=='NM':
                SR1184E_mol.append(float(conc[:-2])/1000)
            else:
                SR1184E_mol.append(float(conc[:-2]))
        elif i.split('/')[1].split('_')[1]=='BULK':
            bulk_dimer.append(load_fit(i+'/'+i.split('/')[1][:-11]+'peaks.dat'))
            conc=i.split('/')[1].split('_')[2]
            if conc[-2:]=='NM':
                bulk_mol.append(float(conc[:-2])/1000)
            else:
                bulk_mol.append(float(conc[:-2]))
      else:
            print i
            AF_dimer.append(load_fit(i+'/'+i.split('/')[1][:-11]+'peaks.dat'))
            conc=i.split('/')[1].split('_')[3].split('.')[0]
            print conc
            if conc[-2:]=='NM':
                AF_mol.append(float(conc[:-2])/1000)
            else:
                AF_mol.append(float(conc[:-2]))
    Kd=1
    Kd_SR1, success= optimize.leastsq(KD_dimer_err, Kd, args=(np.array(SR1_mol),np.multiply(np.array(SR1_dimer),np.array(SR1_mol))/2), maxfev=10000000) 
    Kd_bulk, success= optimize.leastsq(KD_dimer_err, Kd, args=(np.array(bulk_mol),np.multiply(np.array(bulk_dimer),np.array(bulk_mol))/2), maxfev=10000000) 
    Kd_184E, success= optimize.leastsq(KD_dimer_err, Kd, args=(np.array(SR1184E_mol),np.multiply(np.array(SR1184E_dimer),np.array(SR1184E_mol))/2), maxfev=10000000) 
    #kd_AF, sucess=     optimize.leastsq(KD_dimer_err, Kd, args=(np.array(AF_mol),     np.multiply(np.array(AF_dimer),np.array(AF_mol))/2), maxfev=10000000) 
    print Kd_SR1, Kd_bulk,Kd_184E
    #print kd_AF, AF_dimer, AF_mol
    x_sim=np.arange(0.1,25,0.1)
    fit_SR1=[]
    fit_bulk=[]
    fit_184E=[]
    fit_AF=[]
    for i in x_sim:
        fit_SR1.append(2*KD_dimer(i, Kd_SR1)/i)
        fit_bulk.append(2*KD_dimer(i, Kd_bulk)/i)
        fit_184E.append(2*KD_dimer(i, Kd_184E)/i)
        #fit_AF.append(2*KD_dimer(i, kd_AF)/i)
    plt.plot(x_sim,fit_SR1)    
    plt.plot(x_sim,fit_bulk)
    plt.plot(x_sim,fit_184E)
    #plt.plot(x_sim,fit_AF)
    
    #plt.plot(AF_mol,AF_dimer,'o')
    plt.plot(SR1_mol,SR1_dimer,'o')
    plt.plot(SR1184E_mol,SR1184E_dimer,'o')
    plt.plot(bulk_mol,bulk_dimer,'o')
    plt.xscale('log')
    plt.savefig('LBD_kds.eps')
    plt.show()
    



if __name__=='__main__':
    plt.close('all')
    files=glob.glob('Georg/*unidecfiles')
    
    #files=glob.glob('SR1-AF_filutaion/*unidecfiles')
    fit_kds(files)
    #files=glob.glob('190401_SR1_raw/*.txt')
    #for i in files:
    #    mz,inten=load_data(i)
    #    plot_spectrum(mz,inten,i)
    