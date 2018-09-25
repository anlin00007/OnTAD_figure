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
import pyBigWig
def readf (filepath):
	dpTAD_raw = pandas.read_table(filepath,sep='\t',header=None)
        dpTAD_rawa1 = dpTAD_raw.loc[dpTAD_raw[5]>0,:].values[:,0:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[:,2]>0,0:2]-1
        dpTAD_rawb = np.unique(dpTAD_rawa.flatten())
        return (dpTAD_rawa,dpTAD_rawb)
def coverage(domain,ll):
        domain=domain.astype(int)
        for o in range(0,np.shape(domain)[0]):
                for p in xrange(domain[o,0],(domain[o,1]+1)):
                        ll[p]=1
        return ll
def computeMatrix(bwfile,boundarylist,chrn,winsize,res,chrs_l):
    mm = np.zeros((len(boundarylist),2*winsize+1))
    blist = boundarylist[(boundarylist>winsize)&(boundarylist<(chrs_l-winsize*res)/res)]
    for i in range(0,len(blist)):
        for j in range(0,winsize+1):
            mm[i,winsize-j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]-j)*res),int((blist[i]-j+1)*res))))/res)
            mm[i,winsize+j]=max(0,np.nansum(np.array(bwfile.values(chrn, int((blist[i]+j)*res),int((blist[i]+j+1)*res))))/res)
    return mm
def compute_jaccard_index(set_1, set_2):
    n = len(np.intersect1d(set_1,set_2))
    return n / float(len(set_1) + len(set_2) - n)
