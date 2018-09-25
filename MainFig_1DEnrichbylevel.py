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
def statepair_enrichment(state1, state2, ideas_percentage, pos1, pos2, bdict):
    cooccur = (ideas_percentage[pos1,state1]*ideas_percentage[pos2,state2] + ideas_percentage[pos2,state1]*ideas_percentage[pos1,state2])/2
    bg = bdict[state1]*bdict[state2]####(ideas_percentage[pos1,state1]+ideas_percentage[pos2,state1])/2 * (ideas_percentage[pos1,state2]+ideas_percentage[pos2,state2])/2
    if bg == 0:
        enrichscore = 0
    else:
        enrichscore = cooccur/bg
    return enrichscore
def state_enrichment(ideas, boundary):
    pos1 = boundary
    pos1 = pos1.astype(int)
    ideas1 = ideas[pos1,:]
    return ideas1.sum(axis=0)

def coverage(domain,ll):
        domain=domain.astype(int)
        for o in range(0,np.shape(domain)[0]):
                for p in range(domain[o,0],(domain[o,1]+1)):
                        ll[p]=1
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
DPrawexp1 = np.array([0.0]*len(background.keys()))
DPrawexp2 = np.array([0.0]*len(background.keys()))
DPrawexp3 = np.array([0.0]*len(background.keys()))
DPrawexp4 = np.array([0.0]*len(background.keys()))
DPrawexp5 = np.array([0.0]*len(background.keys()))
DPrawexpsolo = np.array([0.0]*len(background.keys()))
DPrawexpnob = np.array([0.0]*len(background.keys()))


DPrawout1=open('/storage/home/lua137/work/TADcalling/tadlevelboundary1.txt','wa')
DPrawout2=open('/storage/home/lua137/work/TADcalling/tadlevelboundary2.txt','wa')
DPrawout3=open('/storage/home/lua137/work/TADcalling/tadlevelboundary3.txt','wa')
DPrawout4=open('/storage/home/lua137/work/TADcalling/tadlevelboundary4.txt','wa')
DPrawout5=open('/storage/home/lua137/work/TADcalling/tadlevelboundary5.txt','wa')
DPrawoutsolo=open('/storage/home/lua137/work/TADcalling/tadlevelboundarysolo.txt','wa')

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
        #print statedict
        dpTAD_raw = pandas.read_table('/storage/home/lua137/work/TADcalling/dpruns/'+celltype+'/10kb/'+'dp_raw_pen0.1_newest.chr'+str(chrnum),sep='\t',header=None)
        dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
	dpTAD_rawfr = getlevel(dpTAD_rawa)
        dpTAD_raw1a = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,3]>1),0:2]
	dpTAD_rawn1a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=1),0:2]
        dpTAD_raw2a = dpTAD_rawfr[dpTAD_rawfr[:,2]==2,0:2]
	dpTAD_rawn2a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=2),0:2]
        dpTAD_raw3a = dpTAD_rawfr[dpTAD_rawfr[:,2]==3,0:2]
	dpTAD_rawn3a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=3),0:2]
        dpTAD_raw4a = dpTAD_rawfr[dpTAD_rawfr[:,2]==4,0:2]
	dpTAD_rawn4a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=4),0:2]
        dpTAD_raw5a = dpTAD_rawfr[dpTAD_rawfr[:,2]>=5,0:2]
	dpTAD_rawn5a = dpTAD_rawfr[(dpTAD_rawfr[:,2]!=5),0:2]
	dpTAD_rawsoloa = dpTAD_rawfr[(dpTAD_rawfr[:,2]==1)&(dpTAD_rawfr[:,3]==1),0:2]

        dpTAD_rawb = np.unique(dpTAD_rawa.flatten())
        dpTAD_raw1b = np.unique(dpTAD_raw1a.flatten())
        dpTAD_raw2b = np.unique(dpTAD_raw2a.flatten())
        dpTAD_raw3b = np.unique(dpTAD_raw3a.flatten())
        dpTAD_raw4b = np.unique(dpTAD_raw4a.flatten())
        dpTAD_raw5b = np.unique(dpTAD_raw5a.flatten())
	
        #dpTAD_rawn1b = np.unique(dpTAD_rawn1a.flatten())
        #dpTAD_rawn2b = np.unique(dpTAD_rawn2a.flatten())
        #dpTAD_rawn3b = np.unique(dpTAD_rawn3a.flatten())
        #dpTAD_rawn4b = np.unique(dpTAD_rawn4a.flatten())
        #dpTAD_rawn5b = np.unique(dpTAD_rawn5a.flatten())

	dpTAD_rawsolob = np.unique(dpTAD_rawsoloa.flatten())
	#print dpTAD_rawa
	#dpTAD_rawnob = np.setdiff1d(np.array(range(0,chrs_length[chrnum-1]/res)),dpTAD_rawb)
	tadarea = np.asarray(coverage(dpTAD_rawa,[0]*(chrs_length[chrnum-1]/res+1)))
	dpTAD_rawnob = np.where(tadarea==0)[0]
	#print tadarea
	print chrs_length[chrnum-1]/res+1,len(dpTAD_rawnob)
	#dpTAD_rawo1b = np.setdiff1d(dpTAD_raw1b,dpTAD_rawn1b)
        #dpTAD_rawo2b = np.setdiff1d(dpTAD_raw2b,dpTAD_rawn2b)
        #dpTAD_rawo3b = np.setdiff1d(dpTAD_raw3b,dpTAD_rawn3b)
        #dpTAD_rawo4b = np.setdiff1d(dpTAD_raw4b,dpTAD_rawn4b)
        #dpTAD_rawo5b = np.setdiff1d(dpTAD_raw5b,dpTAD_rawn5b)
        #print len(np.unique(dpTAD_rawb)),len(dpTAD_rawsolob),len(dpTAD_raw1b),len(dpTAD_raw2b),len(dpTAD_raw3b),len(dpTAD_raw4b),len(dpTAD_raw5b),len(np.intersect1d(dpTAD_raw1b,dpTAD_raw2b)),len(np.intersect1d(dpTAD_raw3b,dpTAD_raw2b)),len(np.intersect1d(dpTAD_raw3b,dpTAD_raw4b)),len(np.intersect1d(dpTAD_raw4b,dpTAD_raw5b)),len(dpTAD_rawo1b),len(dpTAD_rawo2b),len(dpTAD_rawo3b),len(dpTAD_rawo4b),len(dpTAD_rawo5b)

        for i in dpTAD_raw1b:
                DPrawout1.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw2b:
                DPrawout2.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw3b:
                DPrawout3.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw4b:
                DPrawout4.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_raw5b:
                DPrawout5.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')
        for i in dpTAD_rawsolob:
                DPrawoutsolo.write('chr'+str(chrnum)+'\t'+str(int(i)*res)+'\t'+str(int(i+1)*res)+'\t'+str(sum(dpTAD_rawa[:,0]==i))+'\t'+str(sum(dpTAD_rawa[:,1]==i))+'\n')

        DPrawexp += state_enrichment(mat, dpTAD_rawb)
        DPrawexp1 += state_enrichment(mat, dpTAD_raw1b)
        DPrawexp2 += state_enrichment(mat, dpTAD_raw2b)
        DPrawexp3 += state_enrichment(mat, dpTAD_raw3b)
        DPrawexp4 += state_enrichment(mat, dpTAD_raw4b)
        DPrawexp5 += state_enrichment(mat, dpTAD_raw5b)
        DPrawexpsolo += state_enrichment(mat, dpTAD_rawsolob)
	DPrawexpnob += state_enrichment(mat, dpTAD_rawnob)

DPrawobs = np.array([background[k]*np.sum(DPrawexp) for k in background.keys()])
DPrawobs1 = np.array([background[k]*np.sum(DPrawexp1) for k in background.keys()])
DPrawobs2 = np.array([background[k]*np.sum(DPrawexp2) for k in background.keys()])
DPrawobs3 = np.array([background[k]*np.sum(DPrawexp3) for k in background.keys()])
DPrawobs4 = np.array([background[k]*np.sum(DPrawexp4) for k in background.keys()])
DPrawobs5 = np.array([background[k]*np.sum(DPrawexp5) for k in background.keys()])
DPrawobssolo = np.array([background[k]*np.sum(DPrawexpsolo) for k in background.keys()])
DPrawobsnob = np.array([background[k]*np.sum(DPrawexpnob) for k in background.keys()])
#for k in range(0,len(background.keys())):
#	print DPrawexpsolo[k],background.keys()[k]
enrichmat = np.stack([(DPrawexpnob+1)/(DPrawobsnob+1),(DPrawexpsolo+1)/(DPrawobssolo+1),(DPrawexp1+1)/(DPrawobs1+1),(DPrawexp2+1)/(DPrawobs2+1),(DPrawexp3+1)/(DPrawobs3+1),(DPrawexp4+1)/(DPrawobs4+1),(DPrawexp5+1)/(DPrawobs5+1)], axis=0)
#np.savetxt('/storage/home/lua137/work/TADcalling/script/1Denrichmentmat.txt',enrichmat,delimiter='\t')
'''
names=['Gaps','Singletons','level1','level2','level3','level4','level>=5']
enrichmattable = pandas.DataFrame(data=enrichmat,columns=background.keys(),index=names)
sns.set(font_scale=2.5,font='sans-serif')
sns.axes_style({'font.family': ['sans-serif'],
                'font.sans-serif': ['Arial']
                })
sns.set(color_codes=True)
#cmap = "Blues" #sns.cubehelix_palette(light=2, as_cmap=True)#sns.palplot(sns.light_palette("navy", reverse=True))
g = sns.clustermap(enrichmattable,vmin=0,vmax=8,row_cluster=False)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
fi = g.fig
fi.savefig('/storage/home/lua137/work/TADcalling/1Denrichmentheatmapbylevel'+celltype+'.png',dpi=400)
plt.close()
'''
'''
ind = np.arange(len(background.keys()))  # the x locations for the groups
width = 0.1       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, (DIexp+1)/(DIobs+1), width, color='k')
rects2 = ax.bar(ind + width, (Arrowexp+1)/(Arrowobs+1), width, color='y')
rects3 = ax.bar(ind + width*2, (rGMAPexp+1)/(rGMAPobs+1), width, color='m')
rects4 = ax.bar(ind + width*5, (DPrawexp+1)/(DPrawobs+1), width, color='r')
rects5 = ax.bar(ind + width*6, (DPcorexp+1)/(DPcorobs+1), width, color='g')
rects6 = ax.bar(ind + width*3, (TADtreeexp+1)/(TADtreeobs+1), width, color='b')
rects7 = ax.bar(ind + width*4, (ICfinderexp+1)/(ICfinderobs+1), width, color='c')


# add some text for labels, title and axes ticks
ax.set_ylabel('Enrichment')
ax.set_xticks(ind + width*3)
ax.set_xticklabels(background.keys())

ax.legend((rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0],rects7[0]), names)
plt.xticks(rotation=90)
plt.savefig('/storage/home/lua137/work/TADcalling/stateEnrichment.png',dpi=400)
plt.close()
'''
