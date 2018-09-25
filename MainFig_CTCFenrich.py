import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math
import pandas
from StringIO import StringIO
from collections import Counter
from matplotlib.patches import Rectangle
import seaborn as sns
import pyBigWig
def computeMatrix(bwfile,boundarylist,chrn,winsize,res,chrs_l):
    mm = np.zeros((len(boundarylist),2*winsize+1))
    blist = boundarylist[(boundarylist>winsize)&(boundarylist<(chrs_l-winsize*res)/res)]
    for i in range(0,len(blist)):
        for j in range(0,winsize+1):
            mm[i,winsize-j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]-j)*res),int((blist[i]-j+1)*res))))/res)
            mm[i,winsize+j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]+j)*res),int((blist[i]+j+1)*res))))/res)
    return mm
chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
res = 10000
bw = pyBigWig.open("/storage/home/lua137/work/TADcalling/data/ENCODE-rad21.signal.bigWig")
#bw = pyBigWig.open("/storage/home/lua137/work/TADcalling/data/E116-Ctcf.fc.signal.bigwig")
#bw = pyBigWig.open("/storage/home/lua137/work/TADcalling/data/ENCODE-smc3.fc.signal.bigWig")

cDIall = np.empty((0,21))
cArrowall = np.empty((0,21))
crGMAPball = np.empty((0,21))
cdpTAD_rawball = np.empty((0,21))
cdpTAD_rawbl1 = np.empty((0,21))
cdpTAD_rawbl2 = np.empty((0,21))
cdpTAD_corball = np.empty((0,21))
cdpTAD_corbl1 = np.empty((0,21))
cdpTAD_corbl2 = np.empty((0,21))
cTADtreeTADball = np.empty((0,21))
cICTADball = np.empty((0,21))
cdpTAD_rawsoloball = np.empty((0,21))
cdpTAD_rawhierball = np.empty((0,21))
cdpTAD_rawoverball = np.empty((0,21))
cdpTAD_hollowball = np.empty((0,21))
cdpTAD_pknormball = np.empty((0,21))

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

nhier =0
nsolo = 0
for chrnum in range(2,23):
   if chrnum in [9]: 
	continue
   else:
	DItad = pandas.read_table('/storage/home/lua137/work/TADcalling/DI_TAD/hg19/GM12878/10kb/GM12878_10kb_chr'+str(chrnum)+'.add.DI.out.7col.final',sep='\t',header=None)
	DI=DItad.loc[:,1:2].values/res
	DIb=np.unique(DI.flatten())
	DIb = DIb[~np.isnan(DIb)]
	print len(DIb),len(DI)
	Arrowhead = pandas.read_table('/storage/home/lua137/work/TADcalling/juicer/Arrowhead.Gm12878.10kb.KR.chr'+str(chrnum),sep='\t',header=None)
	Arrow=Arrowhead.loc[:,1:2].values/res
	Arrowb=np.unique(Arrow.flatten())
	print len(Arrowb),len(Arrow)
	import glob, os
	flist=[]
	os.chdir('/storage/home/lua137/work/TADcalling/TADtree/final_alg/10kb/Gm12878/chr'+str(chrnum))
	for file in glob.glob("N*.txt"):
	    flist.append(int(file.split('.')[0].split('N')[1]))
	TADtree = pandas.read_table('/storage/home/lua137/work/TADcalling/TADtree/final_alg/10kb/Gm12878/chr'+str(chrnum)+'/N'+str(max(flist))+'.txt',sep='\t',header=0)
	TADtreeTAD = TADtree[['start','end']].values-1
	TADtreeTADb = np.unique(TADtreeTAD.flatten())
	print len(TADtreeTADb),len(TADtreeTAD)
	rG = pandas.read_table('/storage/home/lua137/work/TADcalling/rGMAP/GM12878_combined_10000_chr'+str(chrnum)+'.rGMAPTAD',sep='\t',header=None)
	rGMAP=rG.loc[:,0:1].values/res
	rGMAPb=np.unique(rGMAP.flatten())
	print len(rGMAPb), len(rGMAP)
#	ICFinder = pandas.read_table('/storage/home/lua137/work/TADcalling/IC-Finder/IC-Finder/Gm12878/chr'+str(chrnum)+'.domain',sep='\t',header=None)
#	ICTAD=ICFinder.values-1
#	ICTADb=np.unique(ICTAD.flatten())
#	len(ICTADb),len(ICTAD)
	dpTAD_raw = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/'+'dp_raw_pen0.1_newest.chr'+str(chrnum),sep='\t',header=None)
        dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
        dpTAD_rawfr = getlevel(dpTAD_rawa)
        dpTAD_rawsoloa = dpTAD_rawfr[(dpTAD_rawfr[:,3]==1)&(dpTAD_rawfr[:,2]==1),0:2]
	dpTAD_rawhiera = dpTAD_rawfr[np.logical_or((dpTAD_rawfr[:,3]!=1),(dpTAD_rawfr[:,2]!=1)),0:2]
	dpTAD_rawb = np.unique(dpTAD_rawa.flatten())
        dpTAD_rawsolob1 = np.unique(dpTAD_rawsoloa.flatten())
        dpTAD_rawhierb1 = np.unique(dpTAD_rawhiera.flatten())
	dpTAD_rawsolob = np.setdiff1d(dpTAD_rawsolob1,dpTAD_rawhierb1)
	dpTAD_rawhierb = np.setdiff1d(dpTAD_rawhierb1,dpTAD_rawsolob1)
	dpTAD_rawoverb = np.intersect1d(dpTAD_rawsolob1,dpTAD_rawhierb1)
	print len(np.unique(dpTAD_rawb)),len(dpTAD_rawa)
	nhier += np.shape(dpTAD_rawhiera)[0]
	nsolo += np.shape(dpTAD_rawsoloa)[0]
	
        dpTAD_hollow = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/'+'dp_raw_pen0.1_hollow.chr'+str(chrnum),sep='\t',header=None)
        dpTAD_hollowa1 = dpTAD_hollow.loc[dpTAD_hollow[5]>0,:].values[:,0:3]
        dpTAD_hollowa = dpTAD_hollowa1[dpTAD_hollowa1[:,2]>0,0:2]-1
        dpTAD_hollowb = np.unique(dpTAD_hollowa.flatten())
        print len(np.unique(dpTAD_hollowb)),len(dpTAD_hollowa)

        dpTAD_pknorm = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/'+'dp_raw_pen0.1_pknorm.chr'+str(chrnum),sep='\t',header=None)
        dpTAD_pknorma1 = dpTAD_pknorm.loc[dpTAD_pknorm[5]>0,:].values[:,0:3]
        dpTAD_pknorma = dpTAD_pknorma1[dpTAD_pknorma1[:,2]>0,0:2]-1
        dpTAD_pknormb = np.unique(dpTAD_pknorma.flatten())
        print len(np.unique(dpTAD_pknormb)),len(dpTAD_pknorma)

	cDIall = np.append(cDIall,computeMatrix(bw,DIb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	cArrowall = np.append(cArrowall,computeMatrix(bw,Arrowb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	crGMAPball = np.append(crGMAPball,computeMatrix(bw,rGMAPb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	cdpTAD_rawball = np.append(cdpTAD_rawball,computeMatrix(bw,dpTAD_rawb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawsoloball = np.append(cdpTAD_rawsoloball,computeMatrix(bw,dpTAD_rawsolob,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
        #cdpTAD_rawhierball = np.append(cdpTAD_rawhierball,computeMatrix(bw,dpTAD_rawhierb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawoverball = np.append(cdpTAD_rawoverball,computeMatrix(bw,dpTAD_rawoverb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_corball = np.append(cdpTAD_corball,computeMatrix(bw,dpTAD_corb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_pknormball = np.append(cdpTAD_pknormball,computeMatrix(bw,dpTAD_pknormb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_hollowball = np.append(cdpTAD_hollowball,computeMatrix(bw,dpTAD_hollowb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawbl1 = np.append(cdpTAD_rawbl1,computeMatrix(bw,np.unique(dpTAD_rawa1[dpTAD_rawa1[:,2]==1,0:2].flatten())-1,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawbl2 = np.append(cdpTAD_rawbl2,computeMatrix(bw,np.unique(dpTAD_rawa1[dpTAD_rawa1[:,2]==2,0:2].flatten())-1,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawbl3 = np.append(cdpTAD_rawbl2,computeMatrix(bw,np.unique(dpTAD_rawa1[dpTAD_rawa1[:,2]==3,0:2].flatten())-1,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	#cdpTAD_rawbl4 = np.append(cdpTAD_rawbl2,computeMatrix(bw,np.unique(dpTAD_rawa1[dpTAD_rawa1[:,2]>3,0:2].flatten())-1,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
	cTADtreeTADball = np.append(cTADtreeTADball,computeMatrix(bw,TADtreeTADb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)
#	cICTADball = np.append(cICTADball,computeMatrix(bw,ICTADb,'chr'+str(chrnum),10,10000,chrs_length[chrnum-1]), axis=0)

print nhier,nsolo
out = np.array([np.mean(cdpTAD_rawball, axis=0),np.mean(cDIall, axis=0),np.mean(cArrowall, axis=0),np.mean(crGMAPball, axis=0),np.mean(cTADtreeTADball, axis=0)])
np.savetxt('/storage/home/lua137/work/TADcalling/Rad21enrichment.wholegenome.compare.txt',out)
'''
plt.figure(6)
fig,ax = plt.subplots(1)
ax.plot(np.mean(cdpTAD_rawball, axis=0),c='m',label='Our method (n = {1:d})'
             ''.format(0, np.shape(cdpTAD_rawball)[0]))
#ax.plot(np.mean(cdpTAD_corball, axis=0),c='m',label='Our method (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_corball)[0]))
#ax.plot(np.mean(cdpTAD_pknormball, axis=0),c='g',label='DP_pknorm (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_pknormball)[0]))
ax.plot(np.mean(cDIall, axis=0),label='DIcaller (n = {1:d})'
             ''.format(0, np.shape(cDIall)[0]))
ax.plot(np.mean(cArrowall, axis=0),c='r',label='Arrowhead (n = {1:d})'
             ''.format(0, np.shape(cArrowall)[0]))
#ax.plot(np.mean(cdpTAD_hollowball, axis=0),c='g',label='DP_rawhollowall (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_hollowball)[0]))
#ax.plot(np.mean(cdpTAD_rawhierball, axis=0),c='y',label='DP_rawhierall (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawhierball)[0]))
#ax.plot(np.mean(cdpTAD_rawsoloball, axis=0),c='g',label='DP_rawsoloall (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawsoloball)[0]))
#ax.plot(np.mean(cdpTAD_rawoverball, axis=0),c='b',label='DP_rawoverall (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawoverball)[0]))
#ax.plot(np.mean(cdpTAD_corbl1, axis=0),c='k',linestyle='dashed',label='DP_cor > 1 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_corbl1)[0]))
#ax.plot(np.mean(cdpTAD_rawbl1, axis=0),c='g',linestyle='dashed',label='DP_raw > 1 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawbl1)[0]))
#ax.plot(np.mean(cdpTAD_corbl2, axis=0),c='k',linestyle=':',label='DP_cor > 2 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_corbl2)[0]))
#ax.plot(np.mean(cdpTAD_rawbl2, axis=0),c='b',label='DP_raw = 2 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawbl2)[0]))
#ax.plot(np.mean(cdpTAD_rawbl3, axis=0),c='k',label='DP_raw = 3 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawbl3)[0]))
#ax.plot(np.mean(cdpTAD_rawbl4, axis=0),c='y',label='DP_raw >= 4 (n = {1:d})'
#             ''.format(0, np.shape(cdpTAD_rawbl4)[0]))
ax.plot(np.mean(crGMAPball, axis=0),c='c',label='rGMAP (n = {1:d})'
             ''.format(0, np.shape(crGMAPball)[0]))
ax.plot(np.mean(cTADtreeTADball, axis=0),c='g',label='TADtree (n = {1:d})'
             ''.format(0, np.shape(cTADtreeTADball)[0]))
#ax.plot(np.mean(cICTADball, axis=0),c='y',label='ICFinder (n = {1:d})'
#             ''.format(0, np.shape(cICTADball)[0]))
ax.set_ylabel('CTCF signal')
ax.legend(loc="upper right",prop={'size': 8})
ax.axes.get_xaxis().set_visible(False)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 26}

matplotlib.rc('font', **font)
plt.savefig('/storage/home/lua137/work/TADcalling/CTCFenrichment.wholegenome.compare.png',dpi=400)
plt.close()
'''
