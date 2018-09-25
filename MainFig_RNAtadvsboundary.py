from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas
from collections import Counter
from matplotlib.patches import Rectangle
import seaborn as sns

def computeMatrix(bwfile,boundarylist,chrn,winsize,res,chrs_l):
    mm = np.zeros((len(boundarylist),2*winsize+1))
    blist = boundarylist[(boundarylist>winsize)&(boundarylist<(chrs_l-winsize*res)/res)]
    for i in range(0,len(blist)):
        for j in range(0,winsize+1):
            mm[i,winsize-j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]-j)*res),int((blist[i]-j+1)*res))))/res)
            mm[i,winsize+j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]+j)*res),int((blist[i]+j+1)*res))))/res)
    return mm

def domainlevel(domain,ll):
        domain=domain.astype(int)
        for o in range(0,np.shape(domain)[0]):
                for p in range(domain[o,0],(domain[o,1]+1)):
                        ll[p]+=1
        return ll

def getlevel(tads):
    ftads = tads[(tads[:,1]-tads[:,0]).argsort()[::-1],:]
    rtads = tads[(tads[:,1]-tads[:,0]).argsort(),:]
    flevel = np.ones(len(tads))
    rlevel = np.ones(len(tads))
    for i in range(0,len(tads)):
        rn = []
        fn = []
        for j in range(0,i):
            if rtads[i,0]<=rtads[j,0] and rtads[i,1]>=rtads[j,1]:
               rn.append(rlevel[j])
            if ftads[i,0]>=ftads[j,0] and ftads[i,1]<=ftads[j,1]:
               fn.append(flevel[j])
        if len(rn)>=1:
           rlevel[i] = max(rn)+1
        if len(fn)>=1:
           flevel[i] = max(fn)+1
    return (np.column_stack((ftads,flevel,rlevel[::-1])))

#####Calculate Background Percentage####
celltype = 'Gm12878'
resn='10kb'
res = 10000

chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]

DPrawout1b=open('/storage/home/lua137/work/TADcalling/tadlevelboundary1.txt','wa')
DPrawout2b=open('/storage/home/lua137/work/TADcalling/tadlevelboundary2.txt','wa')
DPrawout3b=open('/storage/home/lua137/work/TADcalling/tadlevelboundary3.txt','wa')
DPrawout4b=open('/storage/home/lua137/work/TADcalling/tadlevelboundary4.txt','wa')
DPrawout5b=open('/storage/home/lua137/work/TADcalling/tadlevelboundary5.txt','wa')
DPrawoutsolob=open('/storage/home/lua137/work/TADcalling/tadlevelboundarysolo.txt','wa')

DPrawout1TAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDE1.txt','wa')
DPrawout2TAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDE2.txt','wa')
DPrawout3TAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDE3.txt','wa')
DPrawout4TAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDE4.txt','wa')
DPrawout5TAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDE5.txt','wa')
DPrawoutsoloTAD=open('/storage/home/lua137/work/TADcalling/tadlevelINSIDEsolo.txt','wa')

DPrawoutnob=open('/storage/home/lua137/work/TADcalling/tadlevelgap.txt','wa')

for chrnum in range(1,23):
        dpTAD_raw = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/'+celltype+'/10kb/'+'dp_raw_pen0.1_newest.chr'+str(chrnum),sep='\t',header=None)
        dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
	dpTAD_rawfr = getlevel(dpTAD_rawa)
	tadarea = np.asarray(domainlevel(dpTAD_rawa,[0]*(chrs_length[chrnum-1]/res+1)))
	dpTAD_rawsolo = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,2]==1),0:2]
	soloarea = np.asarray(domainlevel(dpTAD_rawsolo,[0]*(chrs_length[chrnum-1]/res+1)))
	dpTAD_rawsoloa = np.intersect1d(np.where(tadarea==1)[0],np.where(soloarea==1)[0])
        dpTAD_raw1 = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,3]>1),0:2]
	l1area = np.asarray(domainlevel(dpTAD_raw1,[0]*(chrs_length[chrnum-1]/res+1)))
        dpTAD_raw1a = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,3]>1),0:2]
        dpTAD_raw2a = dpTAD_rawfr[dpTAD_rawfr[:,2]==2,0:2]
        dpTAD_raw3a = dpTAD_rawfr[dpTAD_rawfr[:,2]==3,0:2]
        dpTAD_raw4a = dpTAD_rawfr[dpTAD_rawfr[:,2]==4,0:2]
        dpTAD_raw5a = dpTAD_rawfr[dpTAD_rawfr[:,2]>=5,0:2]
        dpTAD_rawb = np.unique(dpTAD_rawa.flatten())
        dpTAD_raw1b = np.unique(dpTAD_raw1a.flatten())
        dpTAD_raw2b = np.unique(dpTAD_raw2a.flatten())
        dpTAD_raw3b = np.unique(dpTAD_raw3a.flatten())
        dpTAD_raw4b = np.unique(dpTAD_raw4a.flatten())
        dpTAD_raw5b = np.unique(dpTAD_raw5a.flatten())
	dpTAD_rawsolob = np.unique(dpTAD_rawsolo.flatten())
	dpTAD_raw1all = np.intersect1d(np.where(tadarea==1)[0],np.where(l1area==1)[0])
        dpTAD_raw2all = np.where(tadarea==2)[0]
        dpTAD_raw3all = np.where(tadarea==3)[0]
        dpTAD_raw4all = np.where(tadarea==4)[0]
        dpTAD_raw5all = np.where(tadarea>=5)[0]
	
	dpTAD_raw1tad = np.setdiff1d(dpTAD_raw1all,dpTAD_rawb)	
	dpTAD_raw2tad = np.setdiff1d(dpTAD_raw2all,dpTAD_rawb)
	dpTAD_raw3tad = np.setdiff1d(dpTAD_raw3all,dpTAD_rawb)
	dpTAD_raw4tad = np.setdiff1d(dpTAD_raw4all,dpTAD_rawb)
	dpTAD_raw5tad = np.setdiff1d(dpTAD_raw5all,dpTAD_rawb)
	dpTAD_rawsolotad = np.setdiff1d(dpTAD_rawsoloa,dpTAD_rawb)
	dpTAD_rawnob = np.where(tadarea==0)[0]
	#print tadarea
	print chrs_length[chrnum-1]/res+1,len(dpTAD_rawnob)

        for i in dpTAD_raw1b:
                DPrawout1b.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw2b:
                DPrawout2b.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw3b:
                DPrawout3b.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw4b:
                DPrawout4b.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw5b:
                DPrawout5b.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_rawsolob:
                DPrawoutsolob.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')

        for i in dpTAD_rawnob:
                DPrawoutnob.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\n')

        for i in dpTAD_raw1tad:
                DPrawout1TAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw2tad:
                DPrawout2TAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw3tad:
                DPrawout3TAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw4tad:
                DPrawout4TAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw5tad:
                DPrawout5TAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_rawsolotad:
                DPrawoutsoloTAD.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')

