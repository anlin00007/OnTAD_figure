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
celltype = 'Gm12878'
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

DPrawexp = np.array([0.0]*len(background.keys()))
DPrawexp1tad = np.array([0.0]*len(background.keys()))
DPrawexp2tad = np.array([0.0]*len(background.keys()))
DPrawexp3tad = np.array([0.0]*len(background.keys()))
DPrawexp4tad = np.array([0.0]*len(background.keys()))
DPrawexp5tad = np.array([0.0]*len(background.keys()))
DPrawexpsolotad = np.array([0.0]*len(background.keys()))

DPrawexp1b = np.array([0.0]*len(background.keys()))
DPrawexp2b = np.array([0.0]*len(background.keys()))
DPrawexp3b = np.array([0.0]*len(background.keys()))
DPrawexp4b = np.array([0.0]*len(background.keys()))
DPrawexp5b = np.array([0.0]*len(background.keys()))
DPrawexpsolob = np.array([0.0]*len(background.keys()))

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

        DPrawexp1tad += state_enrichment(mat, dpTAD_raw1tad)
        DPrawexp2tad += state_enrichment(mat, dpTAD_raw2tad)
        DPrawexp3tad += state_enrichment(mat, dpTAD_raw3tad)
        DPrawexp4tad += state_enrichment(mat, dpTAD_raw4tad)
        DPrawexp5tad += state_enrichment(mat, dpTAD_raw5tad)
        DPrawexp1b += state_enrichment(mat, dpTAD_raw1b)
        DPrawexp2b += state_enrichment(mat, dpTAD_raw2b)
        DPrawexp3b += state_enrichment(mat, dpTAD_raw3b)
        DPrawexp4b += state_enrichment(mat, dpTAD_raw4b)
        DPrawexp5b += state_enrichment(mat, dpTAD_raw5b)
        DPrawexpsolotad += state_enrichment(mat, dpTAD_rawsolotad)
        DPrawexpsolob += state_enrichment(mat, dpTAD_rawsolob)
	DPrawexpnob += state_enrichment(mat, dpTAD_rawnob)

DPrawobs1tad = np.array([background[k]*np.sum(DPrawexp1tad) for k in background.keys()])
DPrawobs2tad = np.array([background[k]*np.sum(DPrawexp2tad) for k in background.keys()])
DPrawobs3tad = np.array([background[k]*np.sum(DPrawexp3tad) for k in background.keys()])
DPrawobs4tad = np.array([background[k]*np.sum(DPrawexp4tad) for k in background.keys()])
DPrawobs5tad = np.array([background[k]*np.sum(DPrawexp5tad) for k in background.keys()])
DPrawobssolotad = np.array([background[k]*np.sum(DPrawexpsolotad) for k in background.keys()])
DPrawobs1b = np.array([background[k]*np.sum(DPrawexp1b) for k in background.keys()])
DPrawobs2b = np.array([background[k]*np.sum(DPrawexp2b) for k in background.keys()])
DPrawobs3b = np.array([background[k]*np.sum(DPrawexp3b) for k in background.keys()])
DPrawobs4b = np.array([background[k]*np.sum(DPrawexp4b) for k in background.keys()])
DPrawobs5b = np.array([background[k]*np.sum(DPrawexp5b) for k in background.keys()])
DPrawobssolob = np.array([background[k]*np.sum(DPrawexpsolob) for k in background.keys()])
DPrawobsnob = np.array([background[k]*np.sum(DPrawexpnob) for k in background.keys()])
#for k in range(0,len(background.keys())):
#	print DPrawexpsolo[k],background.keys()[k]
enrichmat = np.stack([(DPrawexpnob+1)/(DPrawobsnob+1),(DPrawexpsolob+1)/(DPrawobssolob+1),(DPrawexpsolotad+1)/(DPrawobssolotad+1),(DPrawexp1b+1)/(DPrawobs1b+1),(DPrawexp1tad+1)/(DPrawobs1tad+1),(DPrawexp2b+1)/(DPrawobs2b+1),(DPrawexp2tad+1)/(DPrawobs2tad+1),(DPrawexp3b+1)/(DPrawobs3b+1),(DPrawexp3tad+1)/(DPrawobs3tad+1),(DPrawexp4b+1)/(DPrawobs4b+1),(DPrawexp4tad+1)/(DPrawobs4tad+1),(DPrawexp5b+1)/(DPrawobs5b+1),(DPrawexp5tad+1)/(DPrawobs5tad+1)], axis=0)
np.savetxt('/storage/home/lua137/work/TADcalling/script/1DenrichmentmatTADvsboundary.txt',enrichmat,delimiter='\t')
names=['Gaps','Singletonsb','SingletonsTAD','level1b','level1TAD','level2b','level2TAD','level3b','level3TAD','level4b','level4TAD','level5b','level5TAD']
enrichmattable = pandas.DataFrame(data=enrichmat,columns=background.keys(),index=names)
g = sns.clustermap(enrichmattable,vmin=0,vmax=4,row_cluster=False)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
fi = g.fig
fi.savefig('/storage/home/lua137/work/TADcalling/1DenrichmentheatmapbyTADvsboundary'+celltype+'.png',dpi=400)
plt.close()

