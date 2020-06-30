import glob
import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn
import scipy
import pickle
import collections
import random
from shutil import copyfile

def find_dimers():
    f=open('Ahnert.txt', 'r')
    lines=f.readlines()[0].split('\r')
    f.close()
    dimers=[]
    for i in lines:
        if i.split('\t')[1]=='C2' and i.split('\t')[6]=='C' and i.split('\t')[4]=='1' and i.split('\t')[3]=='1':
            dimers.append(i.split('\t')[0])

    f=open('dimers.txt','w')
    string=''
    for i in dimers:
        string+=i+','
    f.writelines(string[:-1])
    f.close()
    
def find_monomers():
    f=open('Ahnert.txt', 'r')
    lines=f.readlines()[0].split('\r')
    f.close()
    monomers=[]
    for i in lines:
        if i.split('\t')[1]=='M' and i.split('\t')[6]=='C' and i.split('\t')[4]=='1':
            monomers.append(i.split('\t')[0])

    f=open('monomers.txt','w')
    string=''
    for i in monomers:
        string+=i+','
    f.writelines(string[:-1])
    f.close()    
    

def flatten_pdb(pdb):
    '''This function maes dimers out of biological assembly pdbs and strips out ligands'''
    f=open(pdb,'r')
    lines=f.readlines()
    f.close()
    f=open(pdb.split('.')[0]+'.flat.pdb','w')
    for i in lines:
            if i.split()[0]=='ATOM':
                chain_curr= list(i)[21]
                break
    chain='A'
    Atom=0
    for i in lines:
        line=i
        if i.split()[0]=='ATOM' or i.split()[0]=='TER':# or i.split()[0]=='HETATM':
            Atom+=1
        if i.split()[0]=='MODEL' and i.split()[1]=='2':
            chain='B'
            chain_curr='B'
        elif i.split()[0]=='ATOM' and list(i)[21]!=chain_curr:
            chain='B'
            chain_curr='B'
        if i.split()[0]=='ATOM' or i.split()[0]=='TER':# or i.split()[0]=='HETATM':
            line=list(i)
            line[21]=chain
            line=''.join(line)
            
            if Atom<10:
                line=line[:10]+str(Atom)+line[11:]
            elif Atom<100 and Atom>9:
                line=line[:9]+str(Atom)+line[11:]
            elif Atom<1000 and Atom>99:
                line=line[:8]+str(Atom)+line[11:]
            elif Atom<10000 and Atom>999:
                line=line[:7]+str(Atom)+line[11:]
            elif Atom<100000 and Atom>9999:
                line=line[:6]+str(Atom)+line[11:]
            
            f.writelines(line)

        elif i.split()[0]!='MODEL' and i.split()[0]!='ENDMDL' and i.split()[0]!='HETATM':
            
            f.writelines(line)
def strip_monomers(pdb):   
    '''This strips monomers of ligands'''  
    print pdb
    f=open(pdb,'r')
    print 'yes'
    lines=f.readlines()
    f.close()
    f=open(pdb.split('.')[0]+'.strip.pdb', 'w')
    for i in lines:
        if  i.split()[0] !='HETATM':
            f.writelines(i)
    f.close()
    

def make_monomers(pdb):
    '''this makes monomers out of dimers and strips them of all heteroatoms'''
    f=open(pdb,'r')
    lines=f.readlines()
    f.close()
    fA=open(pdb.split('.')[0]+'_A.pdb', 'w')
    fB=open(pdb.split('.')[0]+'_B.pdb', 'w')
    A=0
    B=0
    for i in lines:
        if i.split()[0]!='ATOM':# and i.split()[0] !='HETATM':
            fA.writelines(i)
            fB.writelines(i)
        else :#if i.split()[0] !='HETATM':
            if i.split()[4]=='A':
                fA.writelines(i)
                if i.split()[0] !='HETATM':
                    A+=1
            if i.split()[4]=='B':
                fB.writelines(i)  
                if i.split()[0] !='HETATM':
                    B+=1   
            
    fA.close()
    fB.close()  
    if abs(A-B)>200:
        print 'chains unequal,ignoring', A,B
        os.remove(pdb.split('.')[0]+'_A.pdb')
        os.remove(pdb.split('.')[0]+'_B.pdb')
        return False
    else:
        return True
def write_script(pdb):
    f=open(pdb+'.sh','w')
    f.writelines('source /Applications/ccp4-7.0/bin/ccp4.setup-sh\n') #make sure shell understands ccp4
    f.writelines('areaimol XYZIN %s XYZIN2 %s XYZOUT %s <<eof-area\n' % (pdb, pdb.split('.')[0]+'_A.pdb',pdb.split('.')[0]+'_A.out'))
    f.writelines('DIFFMODE COMPARE\n')
    f.writelines('OUTPUT\n')
    f.writelines('END\n')
    f.writelines('eof-area')
    f.close()
    
def write_monomer_script(pdb):
    f=open(pdb+'.sh','w')
    f.writelines('source /Applications/ccp4-7.0/bin/ccp4.setup-sh\n') #make sure shell understands ccp4
    f.writelines('areaimol XYZIN %s XYZOUT %s <<eof-area\n' % (pdb, pdb[:-4]+'.out'))
    f.writelines('MODE NOHOH\n')
    f.writelines('OUTPUT\n')
    f.writelines('END\n')
    f.writelines('eof-area')
    f.close()    
      
             
def make_monomer_dimer_dic(blastfile):
    '''This finds monomer dimer pairs for Figure 4'''
    f=open(blastfile,'r')
    lines=f.readlines()
    f.close()
    dic={}
    to_delete=set()#store keys of monomers to remove because they are identical to a dimer
    for i in lines:
        monomer=i.split()[0].split('|')[0].split(':')[0]
        
        dimer=i.split()[1].split('|')[0].split(':')[0]
        start1=int(i.split()[6])
        end1=int(i.split()[7])
        start2=int(i.split()[8])
        end2=int(i.split()[9])
        if float(i.split()[2])==100:
            # If they are identical, ignore this pair
            to_delete.add(monomer)
        if monomer not in dic.keys() and float(i.split()[2])<70 and int(i.split()[5])<10 :
            dic[monomer]=[[dimer,[start1,end1,start2,end2]]]
        elif monomer in dic.keys():
            if float(i.split()[2])<70 and dimer not in dic[monomer]:
                dic[monomer].append([dimer,[start1,end1,start2,end2]])

    for i in list(to_delete):
        try:
            del dic[i]
        except KeyError:
            pass
    return dic
    
    
def make_monomer_monomer_dic(blastfile):
    '''Get monomers that are similar to each other?'''
    f=open(blastfile,'r')
    lines=f.readlines()
    f.close()
    dic={}
    for i in lines:
        monomer1=i.split()[0].split('|')[0].split(':')[0]
        monomer2=i.split()[1].split('|')[0].split(':')[0]
        start1=int(i.split()[6])
        end1=int(i.split()[7])
        start2=int(i.split()[8])
        end2=int(i.split()[9])
        
        if monomer1 not in dic.keys() and float(i.split()[2])<98:
            dic[monomer1]=[[monomer2,[start1,end1,start2,end2]]]
        elif float(i.split()[2])<98:
            dic[monomer1].append([monomer2,[start1,end1,start2,end2]])

    return dic    

def calc_areadiff(pdb):
    '''This executes the areaimol scripts'''
    write_monomer_script(pdb)
    subprocess.check_output('bash %s.sh'%(pdb), shell=True)
    os.remove(pdb+'.sh')#cleanup
    
   

def analyze_minimum_distance(picklefile, start, end):
    '''Analzes pickeled output from DynamXL. Uses start and end from alignment as input'''
    f=open(picklefile)
    min_dist=[]
    A=pickle.load(f)
    for d,i in enumerate(A[1][0]):
        temp=[]
        for c,j in enumerate(i):
            if j>0:
                if A[0][c]>=int(start)  or A[0][d]>=int(start):
                        
                        if A[0][c]<=int(end)  or A[0][d]<=int(end):
                                temp.append(j)
        if len(temp)>0:
              min_dist.append(min(temp))
    return np.average(min_dist)


def compare_homologous_monomers(mono_dim_dic,monomer_dic, recalc=False):
    '''This compares the surfaces of monomers to their dimeric homologs'''
    ASA_Dic={'ALA':129.0,'ARG':274.0,'ASN':195.0,'ASP':193.0,'CYS':167.0,'GLU':223.0,'GLN':225.0,'GLY':104.0,'HIS':224.0,'ILE':197.0,'LEU':201.0,'LYS':236.0,'MET':224.0,'PHE':240.0,'PRO':159.0,'SER':155.0,'THR':172.0,'TRP':285.0,'TYR':263.0,'VAL':174.0}
    Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP']
    # what I defined as hydrophobic
    SASA=[]
    used_monomers=set(['1B6C'])
    used_dimers=set(['1B6C'])
    distances_diff=[]

    resi_diff=[]
    f=open('monomer_dimer_output.txt','w')
    f.writelines('Monomer_pdb_code'+' '+'Dimer_pdb_code'+' '+'HydrophobicSASA_monomer-HydrophobicSASA_dissociated_dimer'+' '+'Exposed_hydrophobic_resiudes_monomer-Exposed_hydrophobic_residues_dissociated_dimer'+' '+'Mean_minimum_distance_exposed_hydrophobic_resiudes_monomer-Mean_minimum_distance_exposed_hydrophobic_resiudes_dissociated_dimer'+'\n')

    for i in mono_dim_dic.keys():
        if i not in used_monomers and mono_dim_dic[i][0][0] not in used_dimers:
         try:
            filename='monomers/'+i.lower()+'.pdb'
            filename_dimer='nonredundant/'+mono_dim_dic[i][0][0]+'_A.pdb'    

            start_mono=mono_dim_dic[i][0][1][0]
            end_mono=mono_dim_dic[i][0][1][1]
            start_dim=mono_dim_dic[i][0][1][2]
            end_dim=mono_dim_dic[i][0][1][3]
            tot_recalc=False
            if tot_recalc==True:
                strip_monomers(filename)
                filename='monomers/'+i.lower()+'.strip.pdb' #changing the filename to use the stripped version of the monomers
                calc_areadiff(filename)
                flatten_pdb('nonredundant/'+mono_dim_dic[i][0][0]+'.pdb1')
                make_monomers('nonredundant/'+mono_dim_dic[i][0][0]+'.flat.pdb')
                calc_areadiff('nonredundant/'+mono_dim_dic[i][0][0]+'.flat.pdb')
            filename='monomers/'+i.lower()+'.strip.pdb' #changing the filename to use the stripped version of the monomers

            mono_SASA,mono_resi,mono_start_atom,mono_end_atom =analyze_surfaces(filename,start_mono,end_mono)
            dim_SASA,dim_resi,dimer_start_atom,dimer_end_atom=analyze_surfaces(filename_dimer,start_dim,end_dim)    
            

            monohydro_temp=0
            dimhydro_temp=0
            monohydroSASA_temp=0
            dimhydroSASA_temp=0
            len_mono=end_mono-start_mono
            len_dim=end_dim-start_dim
            mono_hydro_resis=set([]) #This is for DynamXL's calculation
            dimer_hydro_resis=set([]) #This is for DynamXL's calculation

            if abs(len_mono-len_dim)<10: #Don't accept a length difference of more than 10 residues

                for c,j in enumerate(mono_SASA):
                    try:
                        if   float(j)/ASA_Dic[mono_resi[c]]>0.2 and mono_resi[c] in Hydro:
                                        monohydro_temp+=1
                                        monohydroSASA_temp+=j
                                        mono_hydro_resis.add(mono_resi[c]+'.CA')
                    except KeyError:
                        pass #ignore residues with weird names
                for c,j in enumerate(dim_SASA):
                    try:
                        if   float(j)/ASA_Dic[dim_resi[c]]>0.2 and dim_resi[c] in Hydro:
                                        dimhydro_temp+=1        
                                        dimhydroSASA_temp+=j
                                        dimer_hydro_resis.add(dim_resi[c]+'.CA')
                    except KeyError:
                        pass #ignore residues with weird names
                
                recalc=False # if this is set to True, you recalculate the DynamXL distances. This takes a while

                if recalc==True:
                    filename='monomers/'+i.lower()+'.strip.pdb' #changing the filename to use the stripped version of the monomers

                    subprocess.check_output('python ../../../dynamxl_src/dynamxl.py -pdb %s -sel %s -m 2 -o %s.pkl '%(filename, ','.join(mono_hydro_resis),filename), shell=True)
                    subprocess.check_output('python ../../../dynamxl_src/dynamxl.py -pdb %s -sel %s -m 2 -o %s.pkl '%(filename_dimer, ','.join(dimer_hydro_resis),filename_dimer), shell=True)
                    min_dis_mono=analyze_minimum_distance(filename+'.pkl',mono_start_atom,mono_end_atom) 
                    min_dis_dimer=analyze_minimum_distance(filename_dimer+'.pkl',dimer_start_atom, dimer_end_atom)

                else:

                    min_dis_mono=analyze_minimum_distance(filename+'.pkl',mono_start_atom,mono_end_atom) 
                    min_dis_dimer=analyze_minimum_distance(filename_dimer+'.pkl',dimer_start_atom, dimer_end_atom)
                
                
                distances_diff.append(min_dis_dimer-min_dis_mono)
                SASA.append(monohydroSASA_temp-dimhydroSASA_temp)
                resi_diff.append(monohydro_temp-dimhydro_temp)
                used_dimers.add(mono_dim_dic[i][0][0]) #only ever use each dimer once
                f.writelines(i+' '+mono_dim_dic[i][0][0]+' '+str(SASA[-1])+' '+str(resi_diff[-1])+' '+str(distances_diff[-1])+'\n')
                # Write all the pdb files youo used into a separate folder
                copyfile('nonredundant/'+mono_dim_dic[i][0][0]+'.flat.pdb', 'used_files_monomer_dimer/'+mono_dim_dic[i][0][0]+'.flat.pdb')
                copyfile('nonredundant/'+mono_dim_dic[i][0][0]+'.flat.pdb.out', 'used_files_monomer_dimer/'+mono_dim_dic[i][0][0]+'.flat.pdb.out')
                copyfile('nonredundant/'+mono_dim_dic[i][0][0]+'_A.pdb', 'used_files_monomer_dimer/'+mono_dim_dic[i][0][0]+'_A.pdb')
                copyfile('nonredundant/'+mono_dim_dic[i][0][0]+'_A.out', 'used_files_monomer_dimer/'+mono_dim_dic[i][0][0]+'A.out')
                
                copyfile('monomers/'+i.lower()+'.strip.pdb', 'used_files_monomer_dimer/'+i+'.strip.pdb')
                copyfile('monomers/'+i+'.strip.out', 'used_files_monomer_dimer/'+i+'strip.out')
                copyfile('monomers/'+i.lower()+'.strip.pdb.pkl', 'used_files_monomer_dimer/'+i+'.strip.pdb.pkl')

                copyfile(filename+'.pkl', 'used_files_monomer_dimer/'+mono_dim_dic[i][0][0]+'.pkl')

                try:
                    for j in monomer_dic[i]: #don't use any monomer in the homology group again
                        used_monomers.add(j[0])
                except KeyError:
                    pass

            
         except IOError:
            pass
         except ZeroDivisionError:
            pass
         except subprocess.CalledProcessError:
            pass
         
    f.close()
    fig, ax1 = plt.subplots()
    seaborn.set_style("whitegrid", {'axes.grid' : False})
    seaborn.set_style("ticks")
    plt.axvline(np.median(SASA),color='black')
    cum_SASA=[len(sorted(SASA)[:c+1])/float(len(sorted(SASA))) for c,i in enumerate(sorted(SASA))]
    print scipy.stats.wilcoxon(SASA)

    ax1.hist(np.array(SASA))
    plt.xlabel('difference in exposed hydrophobic surface area')
    ax2=plt.twinx()
    ax2.plot(sorted(SASA),cum_SASA, 'k')
    plt.savefig('Figure4D.eps')
    

    fig, ax1 = plt.subplots()
    x=np.array(distances_diff)[~np.isnan(distances_diff)]
    print scipy.stats.wilcoxon(x)
    cum_dist=[len(sorted(x)[:c+1])/float(len(sorted(x))) for c,i in enumerate(sorted(x))]
    ax1.hist(x)
    plt.xlabel('difference in minimum distance betweeb hydrophobics')
    ax2=plt.twinx()
    ax2.plot(sorted(x),cum_dist, 'k')
    plt.axvline(np.median(x),color='black')
    plt.savefig('Figure4F.eps')



    
    fig, ax1 = plt.subplots()
    cum_resi=[len(sorted(resi_diff)[:c+1])/float(len(sorted(resi_diff))) for c,i in enumerate(sorted(resi_diff))]
    ax1.hist(resi_diff, bins=10)
    plt.xlabel('difference in exposed hydrophobic sites')

    ax2=plt.twinx()
    ax2.plot(sorted(resi_diff),cum_resi, 'k')    
    plt.axvline(np.median(resi_diff),color='black')
    plt.savefig('Figure4E.eps')
    


    
    print scipy.stats.wilcoxon(resi_diff),np.std(resi_diff)
    plt.show()
   
   
def count_hydrophobic(hydro_list, resi_list):
    '''Counts ghydrophobic sites in a sequence'''
    cnt=collections.Counter(resi_list)
    hydros=0
    for i in hydro_list:
        hydros+=cnt[i]
    return hydros
 
def get_buried(mono_resi,dim_resi,mono_SASA,dim_SASA):       
    SASA_temp=0
    SASA_tot_temp=0
    no_buried_temp=0    
    ASA_Dic={'ALA':129.0,'ARG':274.0,'ASN':195.0,'ASP':193.0,'CYS':167.0,'GLU':223.0,'GLN':225.0,'GLY':104.0,'HIS':224.0,'ILE':197.0,'LEU':201.0,'LYS':236.0,'MET':224.0,'PHE':240.0,'PRO':159.0,'SER':155.0,'THR':172.0,'TRP':285.0,'TYR':263.0,'VAL':174.0}
    Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP']
    exposed=[]
    buried=[]

    for c,j in enumerate(mono_SASA):
        try:
                if  (j-dim_SASA[c])/ASA_Dic[mono_resi[c]]>0.1: # residue is entirely buried
                    buried.append(mono_resi[c])
                    no_buried_temp+=1
                    if mono_resi[c] in Hydro:
                        SASA_temp+=j-dim_SASA[c]
                    SASA_tot_temp+=j-dim_SASA[c]

                if   float(dim_SASA[c])/ASA_Dic[mono_resi[c]]>0.2 and dim_SASA[c]==j:
                    exposed.append(mono_resi[c])
        except KeyError:
            pass
    return no_buried_temp,buried,SASA_temp,SASA_tot_temp,exposed
    
def get_hydro_fracs(buried,exposed):
    buried_cnt=collections.Counter(buried)
    exposed_cnt=collections.Counter(exposed)
    Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP']

    Exposed_hyro=0
    Buried_hydro=0
    for key in buried_cnt:
            buried_cnt[key]/=float(len(buried)) 
            if key in Hydro:
                Buried_hydro+=buried_cnt[key]
    for key in exposed_cnt:
            exposed_cnt[key]/=float(len(exposed))         
            if key in Hydro:
                Exposed_hyro+=exposed_cnt[key]
    return  Buried_hydro,Exposed_hyro

def tally_surfaces():
    files=glob.glob('nonredundant/*.pdb1')
    os.mkdir('surface_tally_used_files')
    all_exposed=[]
    all_buried=[]
    Frac_interface=[]
    Frac_exposed=[]
    PDB_GC_dic,species_dic=get_species('dimer_uniprot.fasta')
    GCs=[]
    GC_hdyro_frac=[]
    buried_list=[]
    SASA_list=[]
    SASA_tot_list=[]
    buried_tot=[]
    f=open('surface_comparison_output.txt','w')
    f.writelines('PDB_code  Fraction_hydrophobic_sites_at_interface    Fraction_hydrophobic_sites_on_surface'+'\n')
    f1=open('used_dimers_interfaces.txt','w')
    f1.writelines('PDB_code  SASA_buried_at_interface_by_hydrophobic_sites    Number_of_hydrophobic_sites_buried_at_interface'+'\n')

    for i in files:
        SASA_temp=0
        SASA_tot_temp=0
        try:
            mono_SASA,mono_resi,mono_start,mono_end=analyze_surfaces(i[:-5]+'_A.pdb')
            dim_SASA,dim_resi,dim_start,dim_end=analyze_surfaces(i[:-5]+'.flat.pdb')
            if sum(mono_SASA)<0 or sum(dim_SASA)<0:
                print i,'old file, ignoring'
            else:
                no_buried_temp=0    
        
                Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP']
                exposed=[]
                buried=[]
                if len(mono_resi)-len(dim_resi)==0: #if there are the same number of sites in the extracted monomer and one chain of the dimer

                        no_buried_temp,buried,SASA_temp,SASA_tot_temp,exposed=get_buried(mono_resi,dim_resi,mono_SASA,dim_SASA)
                        
                        buried_tot.append(no_buried_temp)
                        buried_list.append(count_hydrophobic(Hydro, buried))
                        SASA_list.append(SASA_temp)
                        SASA_tot_list.append(SASA_tot_temp)
                        all_buried=all_buried+buried
                        all_exposed=all_exposed+exposed

                        
                        
                        A=[]
                        B=[]
                        try:
                            GCs.append(PDB_GC_dic[i.split('\\')[1][:4]])
                        except KeyError:
                            print 'dont have species in GC dic', i.split('\\')[1][:4]
                        Buried_hydro,Exposed_hyro=get_hydro_fracs(buried,exposed)
        
                        Frac_interface.append(Buried_hydro)
                        Frac_exposed.append(Exposed_hyro)
                        ######### These ugly hard coded slashes will fail on a Mac####
                        f1.writelines(i.split('\\')[1].split('.')[0]+' '+str(SASA_list[-1])+' '+str(buried_list[-1])+'\n')
                        f.writelines(i.split('\\')[1].split('.')[0]+' '+str(Frac_interface[-1])+' '+str(Frac_exposed[-1])+'\n')
                        copyfile(i[:-5]+'_A.pdb', 'surface_tally_used_files/'+i.split('\\')[1][:-5]+'_A.pdb')
                        copyfile(i[:-5]+'_A.out', 'surface_tally_used_files/'+i.split('\\')[1][:-5]+'_A.out')
                        copyfile(i[:-5]+'.flat.pdb', 'surface_tally_used_files/'+i.split('\\')[1][:-5]+'.flat.pdb')
                        copyfile(i[:-5]+'.flat.pdb.out', 'surface_tally_used_files/'+i.split('\\')[1][:-5]+'.flat.pdb.out')
                        print i

            
        except IOError:
            print i, 'cant find output'
    aas=[]
    f.close()
    f1.close()
    print scipy.stats.wilcoxon(Frac_interface,Frac_exposed)
    print np.median(Frac_interface), np.median(Frac_exposed)
    print np.average(GCs), len(GCs), len(Frac_interface)

    
    plt.figure()
    plt.plot(buried_tot,Frac_interface,'o', alpha=0.5) 
    plt.xlabel('Number of buried sites at interface')
    plt.ylabel('Fraction of LFWVICM at interface')
    plt.savefig('Ext_data5B.eps')
    
    plt.figure()
    plt.hist(GCs, color='grey')
    plt.ylabel('GC content of source organism')
    plt.ylabel('Number of structures')
    plt.savefig('Ext_data5D.eps')
    
    buried_cnt=collections.Counter(all_buried)
    exposed_cnt=collections.Counter(all_exposed)
    
    for key in exposed_cnt:
        buried_cnt[key]/=float(len(all_buried)) 
        exposed_cnt[key]/=float(len(all_exposed))
        A.append(exposed_cnt[key])
        B.append( buried_cnt[key])
        aas.append(key)
    seaborn.set_style("ticks")
    fig, ax1 = plt.subplots()


    ax2=plt.twinx()
    GC_content_grid(0.2, 0.7, 0.05, buried_tot,ax2)
    equi_frac1=[]
    equi_frac2=[]
    equi_frac3=[]
    equi_frac4=[]
    for i in range (200):
        frac1,frac2,frac3,frac4=run_Equi_model(100,buried_tot)    
        equi_frac1.append(frac1)
        equi_frac2.append(frac2)
        equi_frac3.append(frac3)
        equi_frac4.append(frac4)
    plt.errorbar([41.6,38,50.8,67],[np.median(equi_frac1),np.median(equi_frac2),np.median(equi_frac3),np.median(equi_frac4)],yerr=[np.std(equi_frac1),np.std(equi_frac2),np.std(equi_frac3),np.std(equi_frac4)], fmt='o',color='red')
    #ax3=plt.twiny()
    plt.savefig('Ext_data5C.eps')
    
    print scipy.stats.wilcoxon(Frac_interface,Frac_exposed)
    print np.median(Frac_interface), np.median(Frac_exposed)
    print np.average(GCs), len(GCs), len(Frac_interface) 
    
    
    
    plt.figure()
    plt.hist(Frac_interface,alpha=0.5,bins=20)#,orientation="horizontal")
    plt.hist(Frac_exposed,alpha=0.5,bins=20)#,orientation="horizontal")
    plt.errorbar([np.median(equi_frac1),np.median(equi_frac2),np.median(equi_frac3),np.median(equi_frac4)],[41.6,38,50.8,67],xerr=[np.std(equi_frac1),np.std(equi_frac2),np.std(equi_frac3),np.std(equi_frac4)], fmt='o',color='red')

    plt.xlim(0,1)
    plt.savefig('Figure4C.eps') 
    
    plt.figure()
    plt.hist(np.array(Frac_interface)-np.array(Frac_exposed))
    plt.xlabel('Interface hydrophobicity -surface hydrophobicity')
    plt.ylabel('No of dimers')
    plt.savefig('Ext_data5A.eps')
    

    
     #get values for SR1
    SR1_mono_SASA,SR1_mono_resi,SR1_mono_start,SR1_mono_end=analyze_surfaces('nonredundant/SR1_A.pdb')
    SR1_dim_SASA,SR1_dim_resi,SR1_dim_start,SR1_dim_end=analyze_surfaces('nonredundant/SR1.flat.pdb')
    SR1_no_buried_temp,SR1_buried,SR1_SASA_temp,SR1_SASA_tot_temp,SR1_exposed=get_buried(SR1_mono_resi,SR1_dim_resi,SR1_mono_SASA,SR1_dim_SASA)
    print get_hydro_fracs(SR1_buried,SR1_exposed)
    #
    
    # get values for SR2
    SR2_mono_SASA,SR2_mono_resi,SR2_mono_start,SR2_mono_end=analyze_surfaces('nonredundant/SR2_A.pdb')
    SR2_dim_SASA,SR2_dim_resi,SR2_dim_start,SR2_dim_end=analyze_surfaces('nonredundant/SR2.flat.pdb')
    SR2_no_buried_temp,SR2_buried,SR2_SASA_temp,SR2_SASA_tot_temp,SR2_exposed=get_buried(SR2_mono_resi,SR2_dim_resi,SR2_mono_SASA,SR2_dim_SASA)
    print get_hydro_fracs(SR2_buried,SR2_exposed)

    plt.figure()
    plt.hist(buried_list)
    plt.axvline(x= count_hydrophobic(Hydro, SR1_buried))
    plt.axvline(x= count_hydrophobic(Hydro, SR2_buried))

    plt.xlabel('Number of buried hydrophobic sites at interface')
    plt.ylabel('Number of dimers')
    plt.savefig('Figure4B.eps')
    print 'Fracation more hydrophobic sites',len([i for i in buried_list if i>=count_hydrophobic(Hydro, SR1_buried)])/float(len(buried_list))
    print 'Fracation more hydrophobic sites',len([i for i in buried_list if i>=count_hydrophobic(Hydro, SR2_buried)])/float(len(buried_list))
    print 'total number of dimers', len(SASA_list)
    

    ##### FIgure 4A######
    plt.figure()
    plt.hist(SASA_list) 
    plt.axvline(x= SR1_SASA_temp)
    plt.axvline(x= SR2_SASA_temp)
    print 'Fracation more hydrophobic area',(float(len([i for i in SASA_list if i>=SR1_SASA_temp])))/len(SASA_list)
    print 'Fracation more hydrophobic area',(float(len([i for i in SASA_list if i>=SR2_SASA_temp])))/len(SASA_list)

    plt.xlabel('Buried hydrophobic surface area at interface')
    
    plt.ylabel('Number of dimers')
    plt.savefig('Figure4A.eps')

    
 


    plt.show()
    return A,B,aas
    
def analyze_surfaces(pdb,start=-1000,end=10000):
    #This is just to calculate amino acid propensities at surface or not default for start 
    f=open(pdb[:-4]+'.out','r')
    lines=f.readlines()
    resin=-20
    temp=0
    SASA=[]
    resis=[]
    resicount=1
    start_atom=-10
    end_atom=100000
    for c,i in enumerate(lines):
        if i.split()[0]=='ATOM' :
            if resicount<start: # if you're below the start resiude, do nothing but counting residues. Default for start is so high that you're always counting when it's not specified
                if resin !=float(i[22:26]):
                    resicount+=1
                resin=float(i[22:26])
                
                
            elif resicount>=start and resicount<=end: #if you're above the start resiude, start recording SASA and aa type
                if resicount==start or resicount==1:
                    start_atom=i[5:11]
                if resicount<=end:
                    end_atom=i[5:11]
                if float(i[22:26])<resin and pdb.split('.')[1]=='flat': #if your;re a dimer stop after the first chain
                    return SASA[1:],resis[1:],start_atom,end_atom
                if float(i[22:26])!=resin: #when you encounter a new resiude number, record total for old number
                    resis.append(lines[c-1][17:20])
                    resicount+=1

                    SASA.append(temp)

                    resin=float(i[22:26])
                    temp=float(i[60:66].split()[0])
                else:
                    temp+=float(i[60:66].split()[0])
    return  SASA[1:], resis[1:],start_atom,end_atom
    
def codon_simulator(number,Freqs):
    gencode = {
    'ATA':'ILE', 'ATC':'ILE', 'ATT':'ILE', 'ATG':'MET',
    'ACA':'THR', 'ACC':'THR', 'ACG':'THR', 'ACT':'THR',
    'AAC':'ASN', 'AAT':'ASN', 'AAA':'LYS', 'AAG':'LYS',
    'AGC':'SER', 'AGT':'SER', 'AGA':'ARG', 'AGG':'ARG',
    'CTA':'LEU', 'CTC':'LEU', 'CTG':'LEU', 'CTT':'LEU',
    'CCA':'PRO', 'CCC':'PRO', 'CCG':'PRO', 'CCT':'PRO',
    'CAC':'HIS', 'CAT':'HIS', 'CAA':'GLN', 'CAG':'GLN',
    'CGA':'ARG', 'CGC':'ARG', 'CGG':'ARG', 'CGT':'ARG',
    'GTA':'VAL', 'GTC':'VAL', 'GTG':'VAL', 'GTT':'VAL',
    'GCA':'ALA', 'GCC':'ALA', 'GCG':'ALA', 'GCT':'ALA',
    'GAC':'ASP', 'GAT':'ASP', 'GAA':'GLU', 'GAG':'GLU',
    'GGA':'GLY', 'GGC':'GLY', 'GGG':'GLY', 'GGT':'GLY',
    'TCA':'SER', 'TCC':'SER', 'TCG':'SER', 'TCT':'SER',
    'TTC':'PHE', 'TTT':'PHE', 'TTA':'LEU', 'TTG':'LEU',
    'TAC':'TYR', 'TAT':'TYR', 'TAA':'_', 'TAG':'_',
    'TGC':'CYS', 'TGT':'CYS', 'TGA':'_', 'TGG':'TRP'}

    nucleotides=['A','T','C','G']

    Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP','TYR']
    Seq=''
    usage,total=codon_usage()
    Hydro_aas=0
    total_aas=0
    aas=[]
    seq=''
    for i in range (number):
        codon=''
        codon_accept=False
        while codon_accept==False:
            for j in range (3):
                codon=codon+np.random.choice(nucleotides,p=Freqs)
            seq=seq+codon
            prob=usage[codon]/total[gencode[codon]] # This probability would reflect codon usage, and could be used to implement codon bias in teh GC model. Not currently used
            draw=random.uniform(0, 1)
            if draw<=1:# if you wanted to use codon bias, substituet the 1 for prob here. Currently this accepts every codon
                codon_accept=True
                if gencode[codon]!='_':
                    aas.append(gencode[codon])
                if gencode[codon] in Hydro:
                    Hydro_aas+=1.0
                if gencode[codon]!='_':
                    total_aas+=1
                        
            else:
                codon=''
        Seq=Seq+codon
    
    return Hydro_aas/ total_aas ,aas, seq  
        
def Equi_model(branch_length, number):
    #mouse
    if number==1:
        Q1=np.array([[-0.1435687 ,  0.04495586,  0.07105926,  0.02755359],
       [ 0.04495586, -0.1435687 ,  0.02755359,  0.07105926],
       [ 0.24033653,  0.07128626, -0.35643129,  0.04480851],
       [ 0.07128626,  0.24033653,  0.04480851, -0.35643129]])  
        P=scipy.linalg.expm(Q1*branch_length)

    #yeast Zhu... Petrov 2014 PNAS   
    elif number ==2:
        Q2=np.array([[-0.15834166,  0.03146853,  0.07192807,  0.05494505],
       [ 0.03146853, -0.15834166,  0.05494505,  0.07192807],
       [ 0.17482517,  0.09090909, -0.34165834,  0.07592408],
       [ 0.09090909,  0.17482517,  0.07592408, -0.34165834]])  
        P=scipy.linalg.expm(Q2*branch_length)
 
    
    # E.coli    Lee.... Foster 2012 PNAS
    elif number ==3:
        Q3=np.array([[-0.22222222,  0.03535354,  0.10606061,  0.08080808],
       [ 0.03535354, -0.22222222,  0.08080808,  0.10606061],
       [ 0.17676768,  0.06565657, -0.27777778,  0.03535354],
       [ 0.06565657,  0.17676768,  0.03535354, -0.27777778]])        
        P=scipy.linalg.expm(Q3*branch_length)
    
    #Pseudomonas aeruginosa  Dettman ... Kassen   BMC Genomics 2016
    elif number== 4:
        Q4= np.array([[-0.13613614,  0.02252252,  0.06806807,  0.04554555],
       [ 0.02252252, -0.13613614,  0.04554555,  0.06806807],
       [ 0.26176176,  0.07957958, -0.36386386,  0.02252252],
       [ 0.07957958,  0.26176176,  0.02252252, -0.36386386]])
        P=scipy.linalg.expm(Q4*branch_length)

    return P   

def run_Equi_model(branch_length,buried_list):
    
    size=int(random.choice(buried_list))*3
    #mouse
    Freqs1=[0.291899747,0.29220859,0.20793226,0.207959403]
    #yeast
    Freqs2=[0.31,0.31,0.19,0.19]
    # E.coli K12
    Freqs3=[0.246,0.246, 0.254,0.254]
    #Pseudomonas aeruginosa
    Freqs4=[0.165,0.165,0.335,0.335]
    freq,aas, input_seq=codon_simulator(size, Freqs1)   
    P=Equi_model(branch_length,1)
    hydro_aa1,total_aa1,aa_seq1=evolve_sequence(input_seq,P)
    
    freq,aas, input_seq=codon_simulator(size, Freqs2)   
    P=Equi_model(branch_length,2)
    hydro_aa2,total_aa2,aa_seq2=evolve_sequence(input_seq,P)
    
    freq,aas, input_seq=codon_simulator(size, Freqs3)   
    P=Equi_model(branch_length,3)
    hydro_aa3,total_aa3,aa_seq3=evolve_sequence(input_seq,P)
    
    freq,aas, input_seq=codon_simulator(size, Freqs4)   
    P=Equi_model(branch_length,4)
    hydro_aa4,total_aa4,aa_seq4=evolve_sequence(input_seq,P)
    
    return hydro_aa1/total_aa1,hydro_aa2/total_aa2,hydro_aa3/total_aa3,hydro_aa4/total_aa4,

def evolve_sequence(input_seq,P):
    gencode = {
    'ATA':'ILE', 'ATC':'ILE', 'ATT':'ILE', 'ATG':'MET',
    'ACA':'THR', 'ACC':'THR', 'ACG':'THR', 'ACT':'THR',
    'AAC':'ASN', 'AAT':'ASN', 'AAA':'LYS', 'AAG':'LYS',
    'AGC':'SER', 'AGT':'SER', 'AGA':'ARG', 'AGG':'ARG',
    'CTA':'LEU', 'CTC':'LEU', 'CTG':'LEU', 'CTT':'LEU',
    'CCA':'PRO', 'CCC':'PRO', 'CCG':'PRO', 'CCT':'PRO',
    'CAC':'HIS', 'CAT':'HIS', 'CAA':'GLN', 'CAG':'GLN',
    'CGA':'ARG', 'CGC':'ARG', 'CGG':'ARG', 'CGT':'ARG',
    'GTA':'VAL', 'GTC':'VAL', 'GTG':'VAL', 'GTT':'VAL',
    'GCA':'ALA', 'GCC':'ALA', 'GCG':'ALA', 'GCT':'ALA',
    'GAC':'ASP', 'GAT':'ASP', 'GAA':'GLU', 'GAG':'GLU',
    'GGA':'GLY', 'GGC':'GLY', 'GGG':'GLY', 'GGT':'GLY',
    'TCA':'SER', 'TCC':'SER', 'TCG':'SER', 'TCT':'SER',
    'TTC':'PHE', 'TTT':'PHE', 'TTA':'LEU', 'TTG':'LEU',
    'TAC':'TYR', 'TAT':'TYR', 'TAA':'_', 'TAG':'_',
    'TGC':'CYS', 'TGT':'CYS', 'TGA':'_', 'TGG':'TRP'}
    nucleotide_dic={'A':0,'T':1,'C':2,'G':3}
    nucleotides=['A','T','C','G']

    new_seq=[]
    Hydro=['LEU', 'PHE', 'ILE',  'MET', 'VAL', 'CYS','TRP']
    hydro_aa=0
    total_aa=0
    aa_seq=[]
    for i in input_seq:
        probs=P[nucleotide_dic[i]]
        change=np.random.choice(nucleotides,p=np.array(probs)/sum(probs))
        new_seq.append(change)
    new_seq=''.join(new_seq)
    for i in range(0,len(new_seq)-3,3):
        codon=gencode[new_seq[i]+new_seq[i+1]+new_seq[i+2]]
        if codon in Hydro:
            hydro_aa+=1.0
        if codon!='_':
            total_aa+=1
            aa_seq.append(codon)
    return hydro_aa,total_aa,aa_seq
def codon_usage():
    '''This is a codon bias model that's actually not used in the current implementation'''
    usage={'TTT':17.2,'TCT':16.2,'TAT':12.2,'TGT':11.4,
    'TTC':21.8,'TCC':18.1,'TAC':16.1,'TGC':12.3,'TTA':6.7,
    'TCA':11.8,'TAA':1.0,'TGA':1.6,'TTG':13.4,'TCG':4.2,
    'TAG':0.8,'TGG':12.5,'CTT':13.4,'CCT':18.4,'CAT':10.6,
    'CGT':4.7,'CTC':20.2,'CCC':18.2,'CAC':15.3,'CGC':9.4,
    'CTA':8.1,'CCA':17.3,'CAA':12.0,'CGA':6.6,'CTG':39.5,
    'CCG':6.2,'CAG':34.1,'CGG':10.2,'ATT':15.4,'ACT':13.7,
    'AAT':15.6,'AGT':12.7,'ATC':22.5,'ACC':19.0,'AAC':20.3,
    'AGC':19.7,'ATA':7.4,'ACA':16.0,'AAA':21.9,'AGA':12.1,
    'ATG':22.8,'ACG':5.6,'AAG':33.6,'AGG':12.2,'GTT':10.7,
    'GCT':20.0,'GAT':21.0,'GGT':11.4,'GTC':15.4,'GCC':26.0,
    'GAC':26.0,'GGC':21.2,'GTA':7.4,'GCA':15.8,'GAA':27.0,
    'GGA':16.8,'GTG':28.4,'GCG':6.4,'GAG':39.4,'GGG':15.2}    
    
    gencode = {
    'ATA':'ILE', 'ATC':'ILE', 'ATT':'ILE', 'ATG':'MET',
    'ACA':'THR', 'ACC':'THR', 'ACG':'THR', 'ACT':'THR',
    'AAC':'ASN', 'AAT':'ASN', 'AAA':'LYS', 'AAG':'LYS',
    'AGC':'SER', 'AGT':'SER', 'AGA':'ARG', 'AGG':'ARG',
    'CTA':'LEU', 'CTC':'LEU', 'CTG':'LEU', 'CTT':'LEU',
    'CCA':'PRO', 'CCC':'PRO', 'CCG':'PRO', 'CCT':'PRO',
    'CAC':'HIS', 'CAT':'HIS', 'CAA':'GLN', 'CAG':'GLN',
    'CGA':'ARG', 'CGC':'ARG', 'CGG':'ARG', 'CGT':'ARG',
    'GTA':'VAL', 'GTC':'VAL', 'GTG':'VAL', 'GTT':'VAL',
    'GCA':'ALA', 'GCC':'ALA', 'GCG':'ALA', 'GCT':'ALA',
    'GAC':'ASP', 'GAT':'ASP', 'GAA':'GLU', 'GAG':'GLU',
    'GGA':'GLY', 'GGC':'GLY', 'GGG':'GLY', 'GGT':'GLY',
    'TCA':'SER', 'TCC':'SER', 'TCG':'SER', 'TCT':'SER',
    'TTC':'PHE', 'TTT':'PHE', 'TTA':'LEU', 'TTG':'LEU',
    'TAC':'TYR', 'TAT':'TYR', 'TAA':'_', 'TAG':'_',
    'TGC':'CYS', 'TGT':'CYS', 'TGA':'_', 'TGG':'TRP'}
    AA={}
    for i in gencode.keys():
        try:
            AA[gencode[i]]+=usage[i]
        except KeyError:
            AA[gencode[i]]=usage[i]
    return usage, AA
def stdev(x):
    
    return np.std(x)


    
def run_compare_surfaces() :
    A,B,keys=tally_surfaces()
    plt.figure()
    number=[]
    aas=[]
    
    
    for i in range (100):
        result=run_Equi_model(100)
        aas=aas+result[1]
        number.append(result[0])
        
    random_cnt=collections.Counter(aas)
    
    C=[]
    for key in keys:
        random_cnt[key]/=float(len(aas))
        C.append(random_cnt[key])

    print random_cnt
    plt.hist(number)
    print np.median(number)
    plt.figure()
    plt.plot(B,C,'o')
    print scipy.stats.spearmanr(np.array(A), np.array(C))
    print scipy.stats.spearmanr(np.array(B), np.array(C))

    for i,key in enumerate(keys):
        plt.annotate(key,(B[i],C[i]))
    
    plt.show()   
    
def GC_content_grid(mini, maxi, step,buried_list,ax2):
    GC_content=np.arange(mini,maxi+step,step)
    ratios=[]
    ratios_std=[]
    for i in GC_content: 
        Freqs=[(1-i)/2,(1-i)/2,i/2,i/2]
        ratios_temp=[]
        for j in range(200):
            size=int(random.choice(buried_list))*3
            ratio,aas,seq=codon_simulator(size,Freqs)
            
            ratios_temp.append(ratio)
        ratios.append(np.average(ratios_temp))
        ratios_std.append(np.std(ratios_temp))
    
    print ratios,GC_content
    ax2.errorbar(np.array(GC_content)*100,ratios,yerr=ratios_std, color='black')
    #plt.xlabel('GC content')
    plt.ylabel('Simulated fraction of hydrophobic amino acids')
    #plt.savefig('GC_simulation.eps')
    #plt.show()
        
    
def get_species(fasta_file):
    # gets the species names out of a uniprot fasta file
    f=open(fasta_file,'r')
    lines=f.readlines()
    f.close()
    species=[]
    PDB_species_dic={}
    GC_dic=make_GC_dic('GC_content.txt')
    dimer_map=make_dimer_map('dimer_mapping.tab')
    species_dic={}
    for i in lines:
        if i[0]=='>':
            uniprot_id=i.split('|')[1]
            curr_spec=i.split('OS=')[1].split()[0:2]
            try:
                PDB_species_dic[dimer_map[uniprot_id]]=GC_dic[' '.join(curr_spec)]
                species_dic[dimer_map[uniprot_id]]=' '.join(curr_spec)
            except KeyError:
                pass
            species.append(' '.join(curr_spec))
    print len(species)
    species_cnt=collections.Counter(species)
    print species_cnt
    return PDB_species_dic,species_dic

    

def make_dimer_map(mapping_file):
    f=open(mapping_file,'r')
    lines=f.readlines()
    f.close()
    map_dic={}
    for i in lines[1:]:
        map_dic[i.split()[1]]=i.split()[0].upper()
    return map_dic

def make_GC_dic(GC_file):
    f=open(GC_file,'r')
    lines=f.readlines()
    f.close()
    GC_dic={}
    for i in lines:
        GC_dic[' '.join(i.split()[0:2])]=float(i.split()[2])
    return GC_dic
    

if __name__=='__main__':
    plt.close('all')
    monomer_dic=make_monomer_monomer_dic('monomer_monomer_blast.txt')
    dic=make_monomer_dimer_dic('monomer_dimer_blast.txt')
    compare_homologous_monomers(dic,monomer_dic)
    tally_surfaces()
