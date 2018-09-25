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
def state_enrichment(ideas, boundary):
    pos1 = boundary
    pos1 = pos1.astype(int)
    ideas1 = ideas[pos1,:]
    return ideas1.sum(axis=0)

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
celltype = 'Huvec'
resn='10kb'
res = 10000

ideasfinal = pandas.read_table('/storage/home/lua137/work/TADcalling/data/'+celltype+'_IDEAS36',delimiter='\t')
name = ideasfinal['name'].values
nbin = (ideasfinal['chromEnd'].values - ideasfinal['chromStart'].values+1)/200
d = Counter(list(ideasfinal['name'].values))
for i in range(0,len(nbin)):
    d[name[i]]+=(nbin[i]-1)
background = {}
for n in d.keys():
    background[n] = d[n] / float(sum(nbin))
######Calculate countmat for each chr######
statedict = dict()
chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
for i in range(len(background.keys())):
    statedict[background.keys()[i]] = i
print background.keys()

DPrawexpall = np.array([0.0]*len(background.keys()))
DPrawexp1 = np.array([0.0]*len(background.keys()))
DPrawexp2 = np.array([0.0]*len(background.keys()))
DPrawexp3 = np.array([0.0]*len(background.keys()))
DPrawexp4 = np.array([0.0]*len(background.keys()))
DPrawexp5 = np.array([0.0]*len(background.keys()))
DPrawexpsolo = np.array([0.0]*len(background.keys()))
DPrawexpnob = np.array([0.0]*len(background.keys()))

for chrnum in range(1,23):
        for chrnb in range(chrnum,chrnum+1):
            chrlength = chrs_length[chrnb-1]
            chrnumfull='chr'+str(chrnb)
            chromfinal = ideasfinal.loc[ideasfinal['chrom']==chrnumfull,['chrom','chromStart','chromEnd','name']].values
            mat = np.zeros(shape=(chrlength/res+1,len(background.keys())))
            for ideasrow in range(0,np.shape(chromfinal)[0]):
                for pos in range(int(chromfinal[ideasrow,1]),int(chromfinal[ideasrow,2])+1,200):
                    rownum = pos / res
                    state = chromfinal[ideasrow,3]
                    mat[rownum,statedict[state]] += 1
            #np.set_printoptions(suppress=True)
            #np.savetxt('/Users/linan/Desktop/HiCtest/'+celltype+'_'+chrnum+'_'+resn+'_36states_IDEAS.countmat', mat,fmt='%.14f')
        
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
	dpTAD_raw1a = np.intersect1d(np.where(tadarea==1)[0],np.where(l1area==1)[0])
        dpTAD_raw2a = np.where(tadarea==2)[0]
        dpTAD_raw3a = np.where(tadarea==3)[0]
        dpTAD_raw4a = np.where(tadarea==4)[0]
        dpTAD_raw5a = np.where(tadarea>=5)[0]
	
	dpTAD_rawnob = np.where(tadarea==0)[0]
	#print tadarea
	print chrs_length[chrnum-1]/res+1,len(dpTAD_rawnob)
	dpTAD_all = np.where(tadarea>0)[0]
	DPrawexpall += state_enrichment(mat, dpTAD_all)
        DPrawexp1 += state_enrichment(mat, dpTAD_raw1a)
        DPrawexp2 += state_enrichment(mat, dpTAD_raw2a)
        DPrawexp3 += state_enrichment(mat, dpTAD_raw3a)
        DPrawexp4 += state_enrichment(mat, dpTAD_raw4a)
        DPrawexp5 += state_enrichment(mat, dpTAD_raw5a)
        DPrawexpsolo += state_enrichment(mat, dpTAD_rawsoloa)
	DPrawexpnob += state_enrichment(mat, dpTAD_rawnob)

DPrawobsall = np.array([background[k]*np.sum(DPrawexpall) for k in background.keys()])
DPrawobs1 = np.array([background[k]*np.sum(DPrawexp1) for k in background.keys()])
DPrawobs2 = np.array([background[k]*np.sum(DPrawexp2) for k in background.keys()])
DPrawobs3 = np.array([background[k]*np.sum(DPrawexp3) for k in background.keys()])
DPrawobs4 = np.array([background[k]*np.sum(DPrawexp4) for k in background.keys()])
DPrawobs5 = np.array([background[k]*np.sum(DPrawexp5) for k in background.keys()])
DPrawobssolo = np.array([background[k]*np.sum(DPrawexpsolo) for k in background.keys()])
DPrawobsnob = np.array([background[k]*np.sum(DPrawexpnob) for k in background.keys()])
#for k in range(0,len(background.keys())):
#	print DPrawexpsolo[k],background.keys()[k]
enrichmat = np.stack([(DPrawexpnob+1)/(DPrawobsnob+1),(DPrawexpsolo+1)/(DPrawobssolo+1),(DPrawexp1+1)/(DPrawobs1+1),(DPrawexp2+1)/(DPrawobs2+1),(DPrawexp3+1)/(DPrawobs3+1),(DPrawexp4+1)/(DPrawobs4+1),(DPrawexp5+1)/(DPrawobs5+1),(DPrawexpall+1)/(DPrawobsall+1)], axis=0)
np.savetxt('/storage/home/lua137/work/TADcalling/script/1DenrichmentmatTAD'+celltype+'.txt',enrichmat,delimiter='\t')
names=['Gaps','Singletons','level1','level2','level3','level4','level>=5','AllTAD']
#enrichmattable = pandas.DataFrame(data=enrichmat,columns=background.keys(),index=names)
#g = sns.clustermap(enrichmattable,vmin=0,vmax=4,row_cluster=False)
#plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
#fi = g.fig
#fi.savefig('/storage/home/lua137/work/TADcalling/1DenrichmentheatmapbyTADlevel'+celltype+'.png',dpi=400)
#plt.close()

