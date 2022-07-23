# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 21:23:36 2020

@author: ASUS
"""

#from Bio import AlignIO
#count = AlignIO.convert("opuntia1.aln", "clustal", "opuntia.sth","stockholm")
#print("Converted %i alignments" % count)

import os
from Bio.Align.Applications import ClustalwCommandline
clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="obn1.fasta")

assert(os.path.isfile(clustalw_exe), "Clustal W executable missing")
stdout, stderr = clustalw_cline()

from Bio import AlignIO
align = AlignIO.read("obn1.aln", "clustal")
print(align)

from Bio import Phylo
tree = Phylo.read("obn1.dnd", "newick")
Phylo.draw_ascii(tree)

#from Bio import AlignIO
#count = AlignIO.convert("obn1.aln", "clustal", "obn1.sth","stockholm")
#print("Converted %i alignments" % count)

from Bio import AlignIO
alignment = AlignIO.read("obn1.sth", "stockholm")
print("Number of rows: %i" % len(alignment)) 
#Number of rows: 
for record in alignment:
    print("%s - %s" % (record.seq, record.id))


print("alignment")
print(alignment)
#
import numpy as np
from Bio import AlignIO
alignment = AlignIO.read("obn1.sth", "stockholm")
align_array = np.array([list(rec) for rec in alignment], np.character)
print("Array shape %i by %i" % align_array.shape)

align_array = np.array([list(rec) for rec in alignment], np.character, order="F")

#matriks hasil dari clustalw
print('array')
print (align_array)

#from Bio import AlignIO
#count = AlignIO.convert("obn0.fasta", "fasta", "obn0.sth","stockholm")
#print("Converted %i alignments" % count)

import numpy as np
from Bio import AlignIO
alignment1 = AlignIO.read("obn0.sth", "stockholm")
align_array1 = np.array([list(rec) for rec in alignment1], np.character)
print("Array shape %i by %i" % align_array.shape)

align_array1 = np.array([list(rec) for rec in alignment1], np.character, order="F")
#Bukan hasil dari Clustal atau merupakan populais awal
print('array 2')
print (align_array1)

add=len(align_array[0])-len(align_array1[0])

for i in range(add):
    align_array1 = np.append(align_array1,[['-'],
                                           ['-'],
                                           ['-'],
                                           ['-']],axis=1)

#Membentuk Individu Ke-i (Matriks) berdasarkan 2 matriks sebelumnya

#Generate r dist uniform
#q = int(input("nilai q"))
#
#R = np.zeros ((q,1))
#for i in range(q):
#    for 

import random


#inisialisasi populasi
seedprob = 0.5


mpop = []

#Mengisi matriks ke-i
for i in range(10):
    j=0
    while j<=3:
        R= random.uniform(0,1)
        if R<=seedprob:
            mpop.append(align_array[j])
        else:
            mpop.append(align_array1[j])
        j+=1
#np.savetxt('mpop.txt',mpop,fmt='%s')
        
def ingen(matrix,i,j):
    matrix=[]
    for a in range(i,j):
#        print(a)
        matrix.append(mpop[a])
    return matrix

gen=[]
for i in range(10):
    geno=[]
    gen.append(ingen(geno,i*4,i*4+4))
#    np.savetxt('gen%s.txt'%i,gen[i],fmt='%s')
gennn=np.array(gen)
print('gennya')
#print(gennn)

import numpy as np

iterasii=int(input("banyak iterasi = "))
for x in range (iterasii):
    print('populasi awal')
    print(gennn)
    f=[]
    sub=[]
    for i in range(len(gen)):
        for j in range(len(gen[i])):
            for k in range(j+1,len(gen[i])):
                for l in range(len(gen[i][j])):
                    if gen[i][j][l]!=gen[i][k][l]:
                        e = 1
                    else:
                        e = 4
                    sub.append(e)
                f.append(sub)
                sub=[]
    
    print()
    #print('f')
    #print(f)
    sume=[]
    for i in range(len(f)):
        subsum=[]
        subsum.append(np.sum(f[i]))
        sume.append(subsum)
        
    def infit(matr,i,j):
        matr=[]
        for t in range(i,j):
    #        print(a)
            matr.append(sume[t])
        return matr    
    
    fit=[]
    for i in range(10):
        fitt=[]
        fit.append(infit(fitt,i*6,i*6+6))
        
    sumf=[]
    for i in range(len(fit)):
        subsumf=[]
        subsumf.append(np.sum(fit[i]))
        sumf.append(subsumf)
    #    subsume=[]
    #    subsume.append(subsum(f[i]))
    
    print()
    print('fitness')
    print(sumf)
    #print('sume')
    #print(sume)
    
    #print()
    #import numpy as np
    #print(np.sum(sume))
    
    #srt = np.array(sumf)
    #print(srt)
    
    urut= np.argsort(sorted(range(len(sumf)),key=sumf.__getitem__))
    #def argsort(1):
    #    return sorted(range(len(1)), key=1.__getitem__)
    #for i in range(len(srt)):
    #    c=np.argsort(srt)
    #    urut.append(c)
    #np.argsort(srt)
    #print('urut')
    #print(urut)
    
    urut1=[]
    for i in range(len(urut)):
        urut1.append(urut[i]+1)
    print('urutan')
    print(urut1)
    
    #sorted_list = sorted(sumf)
    ##sorted_list = sorted(sumf, reverse = True)
    #
    #print('rank')
    ##print(sorted_list)
    #srl = np.sort(sorted_list)
    #print(srl)
    ##print(srt)
    
    
    
    #fitrank =np.array(range(10,100,10))
    #print(fitrank)
    #fitrank1 =np.array(range(1,11))
    #print(fitrank1)
    #
    #fr =[]
    #for i in range (len(sumf)):
    ##    frr = []
    #    if sumf [i]<= sumf[i+0]:
    ##        frr.append (fitrank1[i])
    #         fr.append(fitrank1[i])
    #
    #
    #print(fr)       
        
    prob=np.zeros((10,1))
    for i in range(10):
        prob[i] = urut1[i]/len(urut1)
    print('peluang adalah')
    print('==========')
    print(prob)
    print( )  
    
    kum = 0
    kumulatif = np.zeros((10,1))
    for i in range (10):
        kum = kum + prob[i,0]
        kumulatif[i,0]=kum
    print('kumulatif adalah')
    print('==========')
    print(kumulatif)
    print( )  
       
    kumrw = 0
    kumulatifrw = np.zeros((10,1))
    for i in range (10):
        kumrw = (kumulatif[i,0]/kumulatif[9])
        kumulatifrw[i,0]=kumrw
    print('kumulatif rw adalah')
    print('==========')
    print(kumulatifrw)
    print( )   
    
    RR = np.zeros((10,1))
    for i in range (10):
        rr = random.uniform(0,1)
        RR[i,0]=rr
    #print('nilai r')
    #print('==========')
    print(RR)
    
    #RR = np.zeros((10,1))
    #for i in range (10):
    #    r = random.uniform(0,1)
    #    RR[i,0]=r
    #print('nilai r')
    #print('==========')
    #print(RR)
    #print( )
    
    
    
    #print( )
    
    mpopn=[]
    
    looprr=0
    loopkum=0
    
    while looprr<=9:
        while loopkum<=9:
    #        print('looprr ',looprr,'loopkum ',loopkum )
            if RR[looprr]<= kumulatifrw[loopkum]:
                mpopn.append(gen[loopkum])
                looprr+=1
                loopkum=0
            else:
                loopkum+=1
            if looprr==10:
                break
                                         
    ##mpopn = []
    ##for i in range (10):
    ##    j=i
    ##    while j<=9:
    ###        RR = random.uniform(0,1)
    ##        if RR[i]<=kumulatifrw[j]:
    ###            mpopn[i]=gen[j]
    ##            mpopn.append(gen[j])
    ##            break
    #
    #
    
    #print(np.array(gen))        
    print("POPULASI BARU")
    print('==========')
    #print(len(mpopn))
    #print()
    #print(np.array(mpopn))
    mpopna=np.array(mpopn)
    #print(mpopna)
    
                
    ##ic = 2
    ##cross = []
    ##for i in range (ic):
    ##    for j in range (4-ic):
    ##        cr =cross.append(mpopn[0][i+2])
    ##        cr1 =cross.append(mpopn [i][i+2])
    ###        cross.append(cr)
    ##        cross.append(cr1)
    ##
    ##print('coba cros', cross)
    #        
    ##cr=[]
    ##for j in range(len(mpopn)):
    ##    cr.append(mpopn[j])
    ##
    ##for i in range(len(cr)):
    ##    if np.mod(i,2)==0:
    ##        for k in range(2):
    ##            temp=[]
    ##            temp=cr[i][k]
    ##            cr[i][k]=cr[i+1][k]
    ##            cr[i+1][k]=temp
    #
    cr=[]
    for j in range(len(mpopna)):
        cr.append(mpopna[j])
    
    for i in range(len(cr)):
        if np.mod(i,2)==0:
            for k in range(2):
                temp=[]
                temp=cr[i][k]
                cr[i][k]=cr[i+1][k]
                cr[i+1][k]=temp
    print('udah di crossover')
    #print(np.array(cr))
    #
    #
    from random import uniform
    from random import randint
    probmut = 0.7
    rmut0=np.zeros((10,1))
    for i in range (10):
        rmutasi= uniform(0,1)
        rmut0[i,0]=rmutasi
    
    print('rmutasi')
    print(rmut0)
    print()
    
    semp=np.zeros((len(cr),len(cr[0]),len(cr[0][0])),dtype=np.chararray)
    for i in range(len(cr)):
        if rmut0[i,0]<= probmut:
            for j in range(len(cr[i])):
                for m in range(len(cr[i][j])):
                    if cr[i][j][m]!='-':
                        semp[i][j][m]=cr[i][j][m]
                    else:
                        for k in range(m,len(cr[i][j])-1):
                                semp[i][j][k]=cr[i][j][k+1]
                        break
        else:
            semp[i]= cr[i]
            
    print(semp)
    
    sump=np.zeros((len(cr),len(cr[0]),len(cr[0][0])),dtype=np.chararray)
    for i in range(len(cr)):
        t = randint(10,30)
        a=1
        print(t)
        for m in range (len(cr[i])):
            sump[i][m][t]='-'
            if rmut0[i,0]<= probmut:
                for j in range(len(cr[i][m])):
                    if t!=0:
                        for k in range(t):
                            sump[i][m][k]=semp[i][m][k]
                            for l in range(t+1,len(cr[0][0])):
                                sump[i][m][l]=semp[i][m][l-1]
                    elif t==len(cr[i][m])-1:
                        for k in range(t):
                            sump[i][m][k]=semp[i][m][k]
                    else:
                        while a<(len(cr[0][0])):
                            sump[i][m][a]=semp[i][m][a-1]
                            a+=1
                            break
            else:
                sump[i][m]=cr[i][m]
    print(sump)
    
    
    popgab=np.concatenate((gen,sump))
    print('pop gabungan')
    print(popgab)
    
    fb=[]
    subb=[]
    for i in range(len(popgab)):
        for j in range(len(popgab[i])):
            for k in range(j+1,len(popgab[i])):
                for l in range(len(popgab[i][j])):
                    if popgab[i][j][l]!=popgab[i][k][l]:
                        e = 1
                    else:
                        e = 4
                    subb.append(e)
                fb.append(subb)
                subb=[]
    
    print()
    #print('f')
    #print(f)
    sumeb=[]
    for i in range(len(fb)):
        subsumb=[]
        subsumb.append(np.sum(fb[i]))
        sumeb.append(subsumb)
        
    def infitb(matrb,i,j):
        matrb=[]
        for t in range(i,j):
    #        print(a)
            matrb.append(sumeb[t])
        return matrb   
    
    fitb=[]
    for i in range(20):
        fittb=[]
        fitb.append(infitb(fittb,i*6,i*6+6))
        
    sumfb=[]
    for i in range(len(fitb)):
        subsumfb=[]
        subsumfb.append(np.sum(fitb[i]))
        sumfb.append(subsumfb)
    #    subsume=[]
    #    subsume.append(subsum(f[i]))
    
    print()
    print('fitness')
    print(sumfb)
    
    urutb= np.argsort(sorted(range(len(sumfb)),key=sumfb.__getitem__))
    urut1b=[]
    for i in range(len(urutb)):
        urut1b.append(urutb[i]+1)
    print('urutan')
    print(urut1b)
    
    probb=np.zeros((20,1))
    for i in range(20):
        probb[i] = urut1b[i]/len(urut1b)
    print('peluang adalah')
    print('==========')
    print(probb)
    print( )  
    
    kumb = 0
    kumulatifb = np.zeros((20,1))
    for i in range (20):
        kumb = kumb + probb[i,0]
        kumulatifb[i,0]=kumb
    print('kumulatif adalah')
    print('==========')
    print(kumulatifb)
    print( )  
       
    kumrwb = 0
    kumulatifrwb = np.zeros((20,1))
    for i in range (20):
        kumrwb = (kumulatifb[i,0]/kumulatifb[19])
        kumulatifrwb[i,0]=kumrwb
    print('kumulatif rw adalah')
    print('==========')
    print(kumulatifrwb)
    print( )   
    
    RRb = np.zeros((10,1))
    for i in range (10):
        rrb = random.uniform(0,1)
        RRb[i,0]=rrb
    #print('nilai r')
    #print('==========')
    print(RRb)
    
    
    mpopnb=[]
    looprrb=0
    loopkumb=0
    while looprrb<=9:
        while loopkumb<=19:
    #        print('looprr ',looprr,'loopkum ',loopkum )
            if RRb[looprrb]<= kumulatifrwb[loopkumb]:
                mpopnb.append(popgab[loopkumb])
                looprrb+=1
                loopkumb=0
            else:
                loopkumb+=1
            if looprrb==10:
                break
    print('populasi baru final')
    print(mpopnb)
    
    fb2=[]
    subb2=[]
    for i in range(len(mpopnb)):
        for j in range(len(mpopnb[i])):
            for k in range(j+1,len(mpopnb[i])):
                for l in range(len(mpopnb[i][j])):
                    if mpopnb[i][j][l]!=mpopnb[i][k][l]:
                        e = 1
                    else:
                        e = 4
                    subb2.append(e)
                fb2.append(subb2)
                subb2=[]
    
    print()
    #print('f')
    #print(f)
    sumeb2=[]
    for i in range(len(fb2)):
        subsumb2=[]
        subsumb2.append(np.sum(fb2[i]))
        sumeb2.append(subsumb2)
        
    def infitb2(matrb2,i,j):
        matrb2=[]
        for t2 in range(i,j):
    #        print(a)
            matrb2.append(sumeb2[t2])
        return matrb2   
    
    fitb2=[]
    for i in range(10):
        fittb2=[]
        fitb2.append(infitb2(fittb2,i*6,i*6+6))
        
    sumfb2=[]
    for i in range(len(fitb2)):
        subsumfb2=[]
        subsumfb2.append(np.sum(fitb2[i]))
        sumfb2.append(subsumfb2)
    #    subsume=[]
    #    subsume.append(subsum(f[i]))
    
    print()
    print('fitness')
    print(sumfb2)
    
    urutb2= np.argsort(sorted(range(len(sumfb2)),key=sumfb2.__getitem__))
    print(urutb2)
    
    gennn = mpopnb
############################################################################

#from Bio import AlignIO
#count = AlignIO.convert("mpop11.fasta", "fasta", "mpop11.sth","stockholm")
#print("Converted %i alignments" % count)

#alignment3 = AlignIO.read("mpop11.sth", "stockholm")
#align_array = np.array([list(rec) for rec in alignment3], np.character)
#print("Array shape %i by %i" % align_array.shape)
#
#align_array3 = np.array([list(rec) for rec in alignment3], np.character, order="F")
#
##matriks hasil dari clustalw
#print('array')
#print (align_array3)

#number_columns = 23
#sample_string = '-----------------------'
#
#l = [list(sample_string[i:i+number_columns]) for i in range(0, len(sample_string), number_columns)]
#matrix = [s if len(s) == number_columns else s+[None]*(number_columns-len(s)) for s in l]
