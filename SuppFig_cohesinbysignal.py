import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math
import pandas
from scipy.stats import gaussian_kde
from StringIO import StringIO
from collections import Counter
from matplotlib.patches import Rectangle
import pyBigWig

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

def computeMatrix(bf1,tadlist1,chrn,res,hm):
    hm = np.copy(hm)
    bwfile1 = np.copy(bf1)
    mm = np.zeros((np.shape(tadlist1)[0],3))
    tadlist=tadlist1[(tadlist1[:,1]-tadlist1[:,0]).argsort(),:]
    tadlist=tadlist.astype(int)
    for i in range(0,np.shape(tadlist)[0]):
	if abs(tadlist[i,1]-tadlist[i,0])>1:
            #mm[i,0]=max(0,np.nansum(np.array(bwfile1.values(chrn, int(tadlist[i,0]*res),int(tadlist[i,1]*res))))/res/(tadlist[i,1]-tadlist[i,0]))
	    mm[i,0]=max(0,np.nanmean(bwfile1[(tadlist[i,0]+1):tadlist[i,1]]))
	    bwfile1[(tadlist[i,0]+1):tadlist[i,1]]=np.nan
            mm[i,1]=np.nanmean(hm[(tadlist[i,0]+1):tadlist[i,1],(tadlist[i,0]+1):tadlist[i,1]])
	    hm[(tadlist[i,0]+1):tadlist[i,1],(tadlist[i,0]+1):tadlist[i,1]] = np.nan
	    mm[i,2]=tadlist[i,2]
	    #if np.log2(mm[i,0]+1e-10) > 2.5:
	    #	print tadlist[i,:]
    return mm

#mm[domain[o,0]:(domain[o,1]+1),domain[o,0]:(domain[o,1]+1)] = np.nan
   

chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
res = 10000
resn = '10kb'
#bw1 = pyBigWig.open("/storage/home/lua137/work/TADcalling/data/ENCODE-rad21.signal.bigWig")

bedfile = pandas.read_table("/storage/home/lua137/work/TADcalling/data/ENCODE-rad21.signal.bigWig.10kb.tab",sep='\t',header=None)
winfile = pandas.read_table("/storage/home/lua137/work/TADcalling/data/ENCODE-rad21.signal.bigWig.10kb.bed",sep='\t',header=None)

#bedfile = pandas.read_table("/storage/home/lua137/work/TADcalling/data/E116-Ctcf.fc.signal.bigwig.10kb.tab",sep='\t',header=None)
#winfile = pandas.read_table("/storage/home/lua137/work/TADcalling/data/ENCODE-rad21.signal.bigWig.10kb.bed",sep='\t',header=None)

a = winfile.loc[:,0].values

cdpTAD_rawball = np.empty((0,3))
levelall = np.array([])
for chrnum in range(1,23):
	mat = np.load('/storage/home/lua137/scratch/Gm12878/10kb/chr'+str(chrnum)+'_dpmean.npy')
        where_are_NaNs = np.isnan(mat)
        mat[where_are_NaNs] = 0
	bw1 = bedfile.loc[np.where(a==('chr'+str(chrnum)))[0],4].values
	print len(bw1)
	dpTAD_raw = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/'+resn+'/dp_raw_pen0.1_newest.chr'+str(chrnum),sep='\t',header=None)
	dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
        dpTAD_rawfr = getlevel(dpTAD_rawa)
	#print dpTAD_rawfr.astype(int)
	#print dpTAD_rawa1.astype(int)
	#dpTAD_raw1a = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,3]>1),0:2]
        #dpTAD_rawn1a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=1),0:2]
        #dpTAD_raw2a = dpTAD_rawfr[dpTAD_rawfr[:,2]==2,0:2]
        #dpTAD_rawn2a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=2),0:2]
        #dpTAD_raw3a = dpTAD_rawfr[dpTAD_rawfr[:,2]==3,0:2]
        #dpTAD_rawn3a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=3),0:2]
        #dpTAD_raw4a = dpTAD_rawfr[dpTAD_rawfr[:,2]>=4,0:2]
        #dpTAD_rawn4a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=4),0:2]
        #dpTAD_raw5a = dpTAD_rawfr[dpTAD_rawfr[:,2]>=5,0:2]
        #dpTAD_rawn5a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=5),0:2]
        #dpTAD_rawsoloa = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,2]==1),0:2]

	#print len(np.unique(dpTAD_rawb)),len(dpTAD_rawa)
	#cdpTAD_rawball = np.append(cdpTAD_rawball,computeMatrix(bw1,dpTAD_rawa,'chr'+str(chrnum),res,mat), axis=0)

        Arrowhead = pandas.read_table('/storage/home/lua137/work/TADcalling/juicer/Arrowhead.Gm12878.10kb.KR.chr'+str(chrnum),sep='\t',header=None)
        Arrow=Arrowhead.loc[:,1:2].values/res
        Arrowb=np.unique(Arrow.flatten())
        #print len(Arrowb),len(Arrow)
	Arrowsort = Arrow[(Arrow[:,1]-Arrow[:,0]).argsort()]
	cdpTAD_rawball = np.append(cdpTAD_rawball,computeMatrix(bw1,dpTAD_rawfr,'chr'+str(chrnum),res,mat), axis=0)
	#levelall = np.append(levelall,dpTAD_rawfr[:,2][::-1])	

#if len(levelall) != len(cdpTAD_rawball):
#	print 'different length'
#	exit()
cdpTAD_rawball = cdpTAD_rawball[cdpTAD_rawball[:,0]>0,:]
#cdpTAD_rawball = cdpTAD_rawball[cdpTAD_rawball[:,2]==3,:]
#levelall = levelall[cdpTAD_rawball[:,0]>0]
x = cdpTAD_rawball[:,1]
y = np.log2(cdpTAD_rawball[:,0]+1e-10)
from scipy.stats.stats import pearsonr
print pearsonr(x, y)


plt.figure(6)
fig,ax = plt.subplots(1)
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
colordict = {1:'r',2:'g',3:'b',4:'c',5:'m'}
clist=[]
#for i in levelall:
#    if i >=5:
#	clist.append(colordict[5])
#    else:
#	clist.append(colordict[i])

ax.scatter(x,y,c=z,label='n = {1:d}'
             ''.format(0, np.shape(cdpTAD_rawball)[0]))
#ax.scatter(0,0,label='pearsonr={1}'
#            ''.format(2, pearsonr(x, y),c='w'))
ax.set_ylabel('Rad21 signal',size=15)
ax.set_xlabel('TAD mean signal',size=15)
ax.legend(loc="upper right",prop={'size': 12})
plt.savefig('/storage/home/lua137/work/TADcalling/Cohesinbysignalall.png',dpi=400)#,format='pdf')
plt.close()
