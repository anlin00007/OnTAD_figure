source('/storage/home/lua137/work/TADcalling/script/ari.R')
require(zoo)
require(data.table)
res=10000
dis=200
getAUC<-function(ari)
{	AUC=c()
	for (i in 1:dim(ari$pr)[1]){
		o = 0
		x <- c(0,as.numeric(ari$rc[i,]),1)
		y <- c(1,as.numeric(ari$pr[i,]),0)
		for(i in 1:(length(x)-1)){o = o+ (x[i+1]-x[i])*(y[i+1]+y[i])/2}
		AUC <- c(AUC,o)
                #id <- order(x)
                #AUC <- c(AUC,sum(diff(x[id])*rollmean(y[id],2)))
	}
	AUC[is.na(AUC)] <- 0
	return(as.numeric(AUC))
}

dpAUC = rep(0,dis+1)
dpAUCnew = rep(0,dis+1)
dpAUCnew2 = rep(0,dis+1)
ArrowAUC = rep(0,dis+1)
DIAUC = rep(0,dis+1)
rGMAPAUC = rep(0,dis+1)
ICAUC = rep(0,dis+1)
TADtreeAUC = rep(0,dis+1)

n=0
'%ni%' <- Negate('%in%')
for (chrnum in 2:22){
   if (chrnum %ni% c(1,2,9)){
	n = n + 1
	hmat <- fread(paste('/storage/home/lua137/groupz/HiC/Gm12878/10kb/chr',as.character(chrnum),'.matrix',sep=''))
	hmat <- as.matrix(hmat)
	DItad = read.table(paste('/storage/home/lua137/work/TADcalling/DI_TAD/hg19/GM12878/10kb/GM12878_10kb_chr',as.character(chrnum),'.add.DI.out.7col.final',sep=''),sep='\t')
	DI=DItad[,2:3]/res + 1
	#print (head(DItad))
	#print (head(DI))
	Arrowhead = read.table(paste('/storage/home/lua137/work/TADcalling/juicer/Arrowhead.Gm12878.10kb.KR.chr',as.character(chrnum),sep=''),sep='\t')
        Arrow=Arrowhead[,2:3]/res + 1
        #print (head(Arrowhead))
        #print (head(Arrow[order(Arrow[,1]),]))
	rG = read.table(paste('/storage/home/lua137/work/TADcalling/rGMAP/GM12878_combined_10000_chr',as.character(chrnum),'.rGMAPTAD',sep=''),sep='\t')
        rGMAP=rG[,1:2]/res + 1
	fpath=paste('/storage/home/lua137/work/TADcalling/TADtree/final_alg/10kb/Gm12878/chr',chrnum,sep='')
	filenames = list.files(fpath,full.names=TRUE,pattern='N')
	finalname = paste(fpath,'/N',length(filenames),'.txt',sep='')
	TADtree = read.table(finalname,sep='\t',header=T)
	TADtreeTAD = TADtree[,c(2,3)]
	ICFinder = read.table(paste('/storage/home/lua137/work/TADcalling/IC-Finder/IC-Finder/Gm12878/chr',as.character(chrnum),'.domain',sep=''),sep='\t')
        ICTAD=ICFinder
	dpTAD_raw = read.table(paste('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/dp_raw_pen0.1_newest.chr',as.character(chrnum),sep=''),sep='\t')
        dpTAD_rawa1 = dpTAD_raw[dpTAD_raw[,6]>0,1:3]
        dpTAD_rawa = dpTAD_rawa1[dpTAD_rawa1[,3]>0,1:2]
	#dpTAD_rawnew = read.table(paste('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/','dp_raw_pen0.1_pknorm.chr',as.character(chrnum),sep=''),sep='\t')
        #dpTAD_rawnewa1 = dpTAD_rawnew[dpTAD_rawnew[,6]>0,1:3]
        #dpTAD_rawnewa = dpTAD_rawnewa1[dpTAD_rawnewa1[,3]>0,1:2]
	#dpTAD_rawnew2 = read.table(paste('/storage/home/lua137/work/TADcalling/dpruns/Gm12878/10kb/','HitHiCraw_pen0.1_max200_pknorm_chr',as.character(chrnum),'_data1.tad',sep=''),sep='\t')
        #dpTAD_rawnew2a1 = dpTAD_rawnew2[dpTAD_rawnew2[,6]>0,1:3]
        #dpTAD_rawnew2a = dpTAD_rawnew2a1[dpTAD_rawnew2a1[,3]>0,1:2]
	#print(head(dpTAD_raw))
	#print(head(dpTAD_rawa1))
	#print(head(dpTAD_rawa))
	dpari = getadjr2(hmat,dpTAD_rawa)
	#dpnewari = getadjr2(hmat,dpTAD_rawnewa)
	#dpnew2ari = getadjr2(hmat,dpTAD_rawnew2a)
	Arrowari = getadjr2(hmat,Arrow)
	DIari = getadjr2(hmat,DI[complete.cases(DI),])
	rGMAPari = getadjr2(hmat,rGMAP)
	ICari = getadjr2(hmat,ICTAD)
	TADtreeari = getadjr2(hmat,TADtreeTAD)
	#print (dpari$adjr2)
	#print (dpari$pr)
	#print (dpari$rc)
	#print (length(dpAUC))
	#print (dpari$adjr2)
	dpAUC = dpAUC+ dpari$adjr2[,1]#getAUC(dpari)
	#dpAUCnew = dpAUCnew+ dpnewari$adjr2[,1]
        #dpAUCnew2 = dpAUCnew2+ dpnew2ari$adjr2[,1]
	#print (dpAUC)
	ArrowAUC = ArrowAUC + Arrowari$adjr2[,1]#getAUC(Arrowari)
	DIAUC = DIAUC + DIari$adjr2[,1]#getAUC(DIari)
	rGMAPAUC = rGMAPAUC + rGMAPari$adjr2[,1]#getAUC(rGMAPari)
	ICAUC = ICAUC + ICari$adjr2[,1]#getAUC(ICari)
	TADtreeAUC = TADtreeAUC + TADtreeari$adjr2[,1]
	#print (ArrowAUC[19:21])
	#print (dpAUC[19:21])
	#pdf(paste('testAUC1050',chrnum,'.pdf',sep=''))
	#plot((1:dim(dpari$pr)[2])/10, (dpari$pr[20,]+1e-06)/(dpari$oh[20,]+1e-06), type='l', xlab="Quantile",ylab="Precision",col='red',xlim=c(0,1),ylim=c(0,5))
	#lines((1:dim(dpari$pr)[2])/10, (dpari$pr[100,]+1e-06)/(dpari$oh[100,]+1e-06), type='l',col='green')
	#lines((1:dim(dpari$pr)[2])/10, (dpari$pr[150,]+1e-06)/(dpari$oh[150,]+1e-06), type='l',col='blue')
	#lines((1:dim(dpari$pr)[2])/10, (rGMAPari$pr[20,]+1e-06)/(rGMAPari$oh[20,]+1e-06), type='c',col='red')
	#lines((1:dim(dpari$pr)[2])/10, (rGMAPari$pr[100,]+1e-06)/(rGMAPari$oh[100,]+1e-06), type='c',col='green')
        #lines((1:dim(dpari$pr)[2])/10, (rGMAPari$pr[150,]+1e-06)/(rGMAPari$oh[150,]+1e-06), type='c',col='blue')
	#lines((1:dim(dpari$pr)[2])/10, (Arrowari$pr[20,]+1e-06)/(Arrowari$oh[20,]+1e-06), type='p',col='red')
	#lines((1:dim(dpari$pr)[2])/10, (Arrowari$pr[100,]+1e-06)/(Arrowari$oh[100,]+1e-06), type='p',col='green')
        #lines((1:dim(dpari$pr)[2])/10, (Arrowari$pr[150,]+1e-06)/(Arrowari$oh[150,]+1e-06), type='p',col='blue')
	#plot(1:dim(dpari$pr)[1], dpAUC, type='l', xlab="Distance",ylab="AUC",col='red')
	#lines(1:dim(dpari$pr)[1],ArrowAUC, type='l',col='green')
	#lines(1:dim(dpari$pr)[1],DIAUC, type='l',col='blue')
	#lines(1:dim(dpari$pr)[1],rGMAPAUC, type='l',col='yellow')
	#lines(1:dim(dpari$pr)[1],ICAUC, type='l',col='orange')
	#legend(1, 0.1, legend=c("DP", "Arrow","DI","rGMAP","ICFinder"),col=c("red","green", "blue","yellow","orange"), lty=1:2, cex=0.8)
	#dev.off()
}
}
out = rbind(dpAUC,ArrowAUC,DIAUC,rGMAPAUC,TADtreeAUC)/n
write.table(out,file='testadjr2alltest.txt',row.names = F, col.names = F, quote = F,sep = '\t')
#pdf('testadjr2alltest.pdf')
#plot(1:(dis+1), dpAUC/n, type='l', xlab="Distance",ylab="AdjRsquare",col='magenta',ylim=c(0,0.8),xlim=c(0,150))
#lines(1:(dis+1), dpAUCnew/n, type='l',col='green')
#lines(1:(dis+1),ArrowAUC/n, type='l',col='red')
#lines(1:(dis+1),DIAUC/n, type='l',col='blue')
#lines(1:(dis+1),rGMAPAUC/n, type='l',col='cyan')
#lines(1:(dis+1),ICAUC/n, type='l',col='yellow')
#lines(1:(dis+1),TADtreeAUC/n, type='l',col='black')
#legend(100, 0.7, legend=c("Our method", "Arrow","DIcaller","rGMAP","ICFinder","TADtree"),col=c("magenta","red", "blue","cyan","yellow","black"), cex=0.8,lty=2)
#dev.off()
