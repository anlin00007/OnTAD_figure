import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math
import pandas
from StringIO import StringIO
from collections import Counter
from matplotlib.patches import Rectangle
import seaborn as sns
chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
res = 10000
cDIall = 0
cArrowall = 0
crGMAPball = 0
cdpTAD_rawball = 0
cdpTAD_corball = 0
cTADtreeTADball = 0
cICTADball = 0
def coverage(domain,ll):
	domain=domain.astype(int)
        for o in range(0,np.shape(domain)[0]):
                for p in xrange(domain[o,0],(domain[o,1]+1)):
			ll[p]=1
	return ll
nbins = 0 
for chrnum in range(2,23):
   if chrnum == 9: 
	continue
   else:
	DItad = pandas.read_table('/storage/home/lua137/work/TADcalling/DI_TAD/hg19/GM12878/10kb/GM12878_10kb_chr'+str(chrnum)+'.add.DI.out.7col.final',sep='\t',header=None)
	DI=DItad.loc[:,1:2].values/res
	DI=DI[~np.isnan(DI).any(axis=1)]
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
	ICFinder = pandas.read_table('/storage/home/lua137/work/TADcalling/IC-Finder/IC-Finder/Gm12878/chr'+str(chrnum)+'.domain',sep='\t',header=None)
	ICTAD=ICFinder.values-1
	ICTADb=np.unique(ICTAD.flatten())
	len(ICTADb),len(ICTAD)
	dpTAD_raw = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/'+'dp_raw_pen0.1_newest.chr'+str(chrnum),sep='\t',header=None)
	dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
	dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
	dpTAD_rawb = np.unique(dpTAD_rawa.flatten())
	print len(np.unique(dpTAD_rawb)),len(dpTAD_rawa)
	cDIall = cDIall+sum(coverage(DI,[0]*(chrs_length[chrnum-1]/res+1)))
	cArrowall = cArrowall+sum(coverage(Arrow,[0]*(chrs_length[chrnum-1]/res+1)))
	crGMAPball = crGMAPball+sum(coverage(rGMAP,[0]*(chrs_length[chrnum-1]/res+1)))
	cdpTAD_rawball = cdpTAD_rawball+sum(coverage(dpTAD_rawa,[0]*(chrs_length[chrnum-1]/res+1)))
	cTADtreeTADball = cTADtreeTADball+sum(coverage(TADtreeTAD,[0]*(chrs_length[chrnum-1]/res+1)))
	cICTADball = cICTADball+sum(coverage(ICTAD,[0]*(chrs_length[chrnum-1]/res+1)))
        nbins = nbins + chrs_length[chrnum-1]/res+1

print np.array([cDIall,cArrowall,crGMAPball,cdpTAD_rawball,cTADtreeTADball,cICTADball])/float(nbins)
plt.figure(6)
fig,ax = plt.subplots(1)
ax.bar(np.arange(5),np.array([cDIall,cArrowall,crGMAPball,cdpTAD_rawball,cTADtreeTADball])/float(nbins),0.35)
plt.xticks(np.arange(5), ('DI','Arrowhead','rGMAP','dp_raw','TADtree'),size=12)

# Left Y-axis labels
plt.ylabel('% covered by TAD', {'color': 'k', 'fontsize': 15})
plt.yticks((0, 0.2, 0.4, 0.6, 0.8), ('0', '20%', '40%', '60%', '80%'), color='k',size=12)#,fontname="Arial", size=12,fontweight='bold')
plt.savefig('/storage/home/lua137/work/TADcalling/coverage.png',dpi=400)
plt.close()
