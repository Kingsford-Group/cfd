import numpy as np
import csv
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import glob
import sys
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import scipy.stats
#from venn import venn

# read in gene quant file
def readQuantFile(filename):

    labels = []
    expdata = []# {}
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        for line in freader:
            #if float(line[2]) > 5:
                #expdata[line[0]] = float(line[2])
            expdata.append(float(line[2]))
            labels.append(line[0])
    return expdata,labels

def readGeneList(expdata,filename):

    genevals = []
    pvals = []
    geneorder=[]
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        #i=1
        for line in freader:
            #genevals.append(expdata[line[0]])
            pvals.append(float(line[1]))
            geneorder.append(line[0])
            #geneorder.append(i)
            #i +=1
    return genevals,geneorder

def readGenePctFile(filename):
    
    genelabels = []
    pctdata = []
    samplelabels = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if len(genelabels) == 0:
                genelabels = line[1:]
            else:
                pctdata.append([float(x) for x in line[1:]])
                fparts = line[0].split('/')
                samplelabels.append(fparts[-2])
    return pctdata,genelabels,samplelabels

def readPValFile(filename):
    
    pvals = []
    genes = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        for line in freader:
            pvals.append(float(line[1]))
            genes.append(line[0])
    return pvals,genes

def readGOfile(filename):

    data = pd.read_csv(filename)
    return data

def readNoiseSimFile(filename):
    
    noisesim = {}
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter=' ')
        for line in freader:
            if len(line) < 4:
                continue
            pctfound = (float(line[-2])-1)/50
            fparts = line[-1].split('_')
            noiselevel = '0.'+fparts[-2].split('.')[-1]
            if noiselevel in noisesim:
                noisesim[noiselevel].append(pctfound)
            else:
                noisesim[noiselevel] = [pctfound]
    return noisesim

def readPrevGenes(filename):

    hkgenes = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            hkgenes.append(line[0].strip())
    return hkgenes

# plot a histogram of gene expression values
def plotHistogram(expvals,genevals,geneorder,filename):

    p = sns.distplot(expvals.values(), bins=50,kde=False)
    genedata = {'geneexp':genevals,'yvals':[10]*len(genevals),'ordering':geneorder}
    sns.scatterplot(x='geneexp',y='yvals',hue='ordering',data=genedata,legend="full")#,color='red',marker='*',s=10)
    p.set_yscale("log")
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote figure to',filename)

def plotHeatmap(expvals,orderidx,samplenames,filename):

    # need data to be in numpy array with rows as samples, columns as expression values
    # reorder columns by ordering, then plot heatmap
    #print(len(orderidx))
    reordereddata = expvals[:,orderidx]

    plt.figure(figsize=(20,10))
    #print(reordereddata[1:5,165:175])
    #sys.exit()
    fig = sns.heatmap(reordereddata[:,:500],yticklabels=samplenames,xticklabels=[],cmap='Blues')
    fig.figure.axes[-1].set_ylabel('Fraction of genes in sample with higher expression', size=16)
    cbar = fig.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    plt.tick_params(labelsize=16)
    plt.axvline(168,0,31,color='red')
    fig.set_xlabel('top 500 genes ordered by p-value',fontsize=16)
    #fig.xaxis.set_label_position('top') 
    fig.set_ylabel('tissue sample',fontsize=18)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote heatmap figure to',filename)

def plotThreshData(threshdata,filename):

    datatuples = threshdata.items()
    labels = [[x[0]]*len(x[1]) for x in datatuples if x[0] in ['None','60','80']]
    labels = [x for l in labels for x in l[:149]]
    #tdata = []
    #for i,dataset in enumerate(datatuples):
        #dataset[1].sort()
        #datatuples[i][1] = dataset[1]
    tdata = [x for l in datatuples for x in l[1][:149] if l[0] in ['None','60','80']]
    threshdict = {'labels':labels,'pvals':tdata}
    fig,ax = plt.subplots()
    fig = sns.stripplot(x='labels',y='pvals',data=threshdict,order=['None','60','80'],palette=['#8dd3c7','#8dd3c7','#fb8072'])
    plt.ylim(0,0.0045) # full data
    #plt.ylim(0,0.00015) # zoomed in
    fig.set_xticklabels(['Poor (all data)', 'Medium (>60%)', 'High (>80%)'])
    legendelems = [matplotlib.lines.Line2D([0],[0],marker='o',markerfacecolor='#8dd3c7',label='does not pass multiple hypothesis testing',color='w'),matplotlib.lines.Line2D([0],[0],marker='o',markerfacecolor='#fb8072',label='passes multiple hypothesis testing',color='w')]
    ax.legend(handles=legendelems)
    #ax.legend(legendcols,['does not pass multiple hypothesis testing','passes multiple hypothesis testing'])
    fig.set_xlabel('Data quality (mapping percentage threshold)')
    fig.set_ylabel('p-values')
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote threshold figure to',filename)

def plotPvalHist(threshdata,filename):

    #threshdata is dictionary with keys as noise level strings ('None', '50','60', etc), and entries are lists of pvalues

    fig = sns.distplot(threshdata['None'],kde=False,color='#8dd3c7',hist=True)
    sns.distplot(threshdata['60'],kde=False,color='#bebada',hist=True)
    sns.distplot(threshdata['80'],kde=False,color='#fb8072',hist=True)
    fig.set_xlabel('p-values',fontsize=14)
    fig.set_ylabel('count',fontsize=14)
    plt.tick_params(labelsize=14)
    plt.xlim(0,0.9)
    fig.legend(['Poor (all data)','Medium (> 60%)','High (> 80%)'],title='Data quality',fontsize=14)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote threshold figure to',filename)


def plotNullData(truedata,nulldata,filename):

    #truedata = truedata[:149]
    labels = [['High Quality Samples']*len(truedata)]
    #nulllabels = [['null'+str(i+1)]*len(x) for i,x in enumerate(nulldata)]
    nulllabels = [['Random Subsets']*len(x) for i,x in enumerate(nulldata)]
    setlabels = [[i]*len(x) for i,x in enumerate(nulldata)]
    setlabels = labels+setlabels
    labels = labels+nulllabels
    setlabels = [x for l in setlabels for x in l]
    labels = [x for l in labels for x in l] # add [:149] for just top pvals
    for nullm in nulldata:
        fig = sns.distplot(nullm,kde=False,color='#8dd3c7',hist=True,hist_kws={"alpha":0.2})
    sns.distplot(truedata,kde=False,color='#fb8072',hist=True)
    fig.set_xlabel('p-values',fontsize=14)
    fig.set_ylabel('count',fontsize=14)
    plt.tick_params(labelsize=14)
    plt.xlim(0,0.9)
    fig.legend(['Random subsets','High quality (> 80%)'],fontsize=14)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote null figure to',filename)    

def plotVennDiag(geneset1, geneset2, labels, filename):

    v = venn2([set(geneset1),set(geneset2)],set_labels=labels)
    v.get_patch_by_id('10').set_color('#a6cee3') # human body map color
    v.get_patch_by_id('01').set_color('#b2df8a') # gtex color
    v.get_patch_by_id('11').set_color('grey') # overlap color - make this prettier
    #venn2(subsets=(64,45,104),set_labels=('HumanBodyMap','GTEx'))
    #v.get_patch_by_id('10').set_color('#fb8072') # chang et al color
    #v.get_patch_by_id('01').set_color('#8dd3c7') # eisenberg and levanon color
    #v.get_patch_by_id('11').set_color('grey')

    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote venn diagram figure to',filename)

def plotComparisons(compdata,filename):

    color2011 = '#ca0020'
    color2013 = '#f4a582'
    color2017 = '#0571b0'

    plt.figure(figsize=(6,6))
    #plt.rcParams["axes.labelsize"] = 14
    #totals = compdata.sum()
    fig = sns.barplot(x=['HBM','GTEx'],y=[compdata[x].sum() for x in ['HBM','GTEx']],color='grey') # none
    sns.barplot(x=['HBM','GTEx'],y=[compdata[x][:-1].sum() for x in ['HBM','GTEx']],color=color2011) # 2011
    sns.barplot(x=['HBM','GTEx'],y=[compdata[x][:-2].sum() for x in ['HBM','GTEx']],color=color2013) # 2013
    sns.barplot(x=['HBM','GTEx'],y=[compdata[x][:-3].sum() for x in ['HBM','GTEx']],color=color2013,linewidth=0) # 2011&2013
    sns.barplot(x=['HBM','GTEx'],y=[compdata[x][:-4].sum() for x in ['HBM','GTEx']],color=color2011,linewidth=0) # 2011&2017
    allbar = sns.barplot(x=['HBM','GTEx'],y=[compdata[x][0] for x in ['HBM','GTEx']],color='black',linewidth=0) #all 3
    allbar.patches[6].set_hatch('//')
    allbar.patches[6].set_edgecolor(color2011)
    allbar.patches[7].set_hatch('//')
    allbar.patches[7].set_edgecolor(color2011)
    allbar.patches[8].set_hatch('//')
    allbar.patches[8].set_edgecolor(color2017)
    allbar.patches[9].set_hatch('//')
    allbar.patches[9].set_edgecolor(color2017)
    matplotlib.rcParams['hatch.linewidth'] = 3
    neitherbar = plt.Rectangle((0,0),1,1,fc="grey", edgecolor = 'none')
    bar2011 = plt.Rectangle((0,0),1,1,fc=color2011,  edgecolor = 'none')
    bar2013 = plt.Rectangle((0,0),1,1,fc=color2013,  edgecolor = 'none')
    bothbar1113 = plt.Rectangle((0,0),1,1,fc=color2013,  edgecolor = color2011,linewidth=0,hatch='//')
    bothbar1117 = plt.Rectangle((0,0),1,1,fc=color2017,  edgecolor = color2011,linewidth=0,hatch='//')
    allbar = plt.Rectangle((0,0),1,1,fc="black", edgecolor = 'none')

    plt.legend([neitherbar,bar2011,bar2013,bothbar1113,bothbar1117,allbar],['None','2011 study only','2013 study only','2011 and 2013','2011 and 2017','all 3 previous studies'],fontsize=12)
    plt.ylim(0,235)
    #fig = sns.barplot(x="study",y="numberofgenes",hue="dataset",data=compdata,palette=['#a6cee3','#b2df8a','grey'])
    fig.set_xlabel('')
    fig.set_xticklabels(['HumanBodyMap','GTEx'],fontsize=14)
    fig.set_ylabel('Number of genes in other HK gene lists',fontsize=14)
    plt.tick_params(labelsize=14)
    #plt.xticks(rotation=30)
    #fig.set_xticklabels(['Chang et al., 2011','Eisenberg and\n Levanon, 2013','Caracausi et al., 2017'],fontsize=14)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote comparison figure to',filename)

def plotTopHKGenesComp(compdata,filename):

    #plt.figure(figsize=(15,15))
    fig = sns.barplot(x="study",y="numberofgenes",hue="dataset",data=compdata,palette=['#a6cee3','#b2df8a','grey'])
    fig.set_xlabel('')
    fig.set_ylabel('Size of overlap in gene sets',fontsize=14)
    plt.xticks(rotation=15)
    fig.set_xticklabels(['Chang et al., 2011','Eisenberg and\n Levanon, 2013','Caracausi et al., 2017'],fontsize=14)
    plt.tick_params(labelsize=14)
    plt.setp(fig.get_legend().get_texts(), fontsize='14') # for legend text
    plt.setp(fig.get_legend().get_title(), fontsize='14') # for legend title
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote comparison figure to',filename)


def plotVennDiagPrevHK(hkgenelists, labels, filename):

    color2011 = '#0571b0'
    color2013 = '#f4a582'
    color2017 = '#ca0020'

    v = venn3([set(x) for x in hkgenelists],set_labels=labels)
    v.get_patch_by_id('100').set_color(color2017)
    v.get_patch_by_id('001').set_color(color2011)
    v.get_patch_by_id('010').set_color(color2013)
    v.get_patch_by_id('101').set_color('grey')
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote venn diagram figure to',filename)

def plotNoiseSimResults(noisesim, filename):

    datalist = []
    for key,vallist in noisesim.iteritems():
        for val in vallist:
            datalist.append((key,val))
    noisedf = pd.DataFrame(datalist,columns=['Noise fraction','Accuracy'])
    fig = sns.barplot(x="Noise fraction",y="Accuracy",data=noisedf,order=['0.0','0.05','0.1','0.15','0.2','0.25'],palette=['#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494'])
    fig.set_xlabel('Fraction of noisy samples',fontsize=14)
    fig.set_ylabel('Accuracy',fontsize=14)
    plt.tick_params(labelsize=14)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('wrote simulation of noise data to',filename)

def main(infile,genefile,pvalfile,threshfiles,nullfiles,gofiles,prevresultscomp,noisefile,outseed):

    if len(infile) > 0 and len(genefile)>0:
        if '*' in infile:
            infilelist = glob.glob(infile)
            fulldata = []
            for inf in infilelist:
                expdata,labels = readQuantFile(inf)
                fulldata.append(expdata)
                genevals,neworder = readGeneList(expdata,genefile)
        else:
            fulldata,labels,samples = readGenePctFile(infile)
            genevals,neworder = readGeneList([],genefile)
            samples = ['thyroid 1','testis 1','ovary 1','leukocyte 1','skeletal muscle 1','prostate 1','lymph node 1','lung 1','adipose 1','adrenal gland 1','brain 1','breast 1','colon 1','kidney 1','heart 1','liver 1','adipose 2','adrenal gland 2', 'brain 2','breast 2', 'colon 2', 'kidney 2','heart 2','liver 2','lung 2','lymph node 2','prostate 2','skeletal muscle 2','leukocyte 2','ovary 2','testis 2','thyroid 2']
            reordersamples = [0,31,1,30,2,29,3,28,4,27,5,26,6,25,7,24,8,16,9,17,10,18,11,19,12,20,13,21,14,22,15,23]
            samples = [samples[x] for x in reordersamples]
        #print(infilelist)
        reordering = [0]*len(labels)
        nomatchidx = len(neworder)
        for i,x in enumerate(labels):
            try:
                reordering[i] = neworder.index(x)
            except:
                #need to put at the end
                reordering[i] = nomatchidx
                nomatchidx += 1
        #reordering = [neworder.index(x) for x in labels]
        idx = np.argsort(reordering)
        fulldata = np.array(fulldata)
        fulldata = fulldata[reordersamples,:]

        #plotHistogram(expdata,genevals,pvals,outfile)
        plotHeatmap(fulldata,idx,samples,outseed+'_top500genes_heatmap.pdf')

    if len(threshfiles) > 0:
        # plot boxplots of pvalues from each threshold value
        threshfilelist = glob.glob(threshfiles)
        threshdata = {}
        for f in threshfilelist:
            fparts = f.split('_')
            threshlabel = fparts[-2][6:]
            threshdata[threshlabel],_ = readPValFile(f)
        #plotThreshData(threshdata,outseed+'_pvalues_bythreshold.png')
        plotPvalHist(threshdata,outseed+'_pvalues_histogram.pdf')

    if len(nullfiles) > 0:

        nullfilelist = glob.glob(nullfiles)
        nulldata = []
        for f in nullfilelist:
            nullpvals,_ = readPValFile(f)
            nulldata.append(nullpvals)
        truedata,_ = readPValFile(pvalfile[0])
        plotNullData(truedata,nulldata,outseed+'_pvalues_nullmodels_histogram.pdf')
        
    if len(gofiles) > 0:
        for fname in gofiles:
            if 'HBM' in fname:
                hbmdata = readGOfile(fname)
                hbmlabel = ['HBM']*len(hbmdata)
                hbmdata['Data source'] = hbmlabel
            elif 'GTEx' in fname:
                gtexdata = readGOfile(fname)
                gtexlabel = ['GTEx']*len(gtexdata)
                gtexdata['Data source'] = gtexlabel
        allgodata = pd.concat([hbmdata,gtexdata],axis=0,ignore_index=True)
        #allgodata.reset_index()
        plotGOresults(allgodata,outseed+'_goresults.png')
        
    if len(pvalfile) > 1:
        for f in pvalfile:
            if 'hbm' in f:
                _,hbmgenes = readPValFile(f)
            elif 'gtex' in f:
                _,gtexgenes = readPValFile(f)
            else:
                print('dont recognize filename for pval file')
        #plotVennDiag(hbmgenes,gtexgenes,('HumanBodyMap','GTEx'),outseed+'_venndiag.png')
        plotVennDiag(hbmgenes,gtexgenes,('',''),outseed+'_venndiag.pdf')

    if len(prevresultscomp) > 0:
        prevresults_topgenes = pd.read_csv(prevresultscomp+'.txt',delimiter='\t')
        prevresults_all = pd.read_csv(prevresultscomp+'_fulllists.txt',delimiter='\t')
        plotTopHKGenesComp(prevresults_topgenes,outseed+'_prevresults_topgenes_barplot.pdf')
        plotComparisons(prevresults_all,outseed+'_prevresults_barplot.pdf')
        #plotComparisons(prevresults,outseed+'_prevresults_barplot.png')
        hkgenelists = [['RPL41','EEF1A1','RPL23A','TPT1','ACTB','RPS23','RPS3A','RPS10','RPS18','RPL27A'],['C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29'],['ACTG1','RPS18','POM121C','MRPL18','TOMM5','YTFDH1','TPT1','RPS27']]
        #labels=('Chang et al., 2011','Eisenberg and Levanon, 2013', 'Caracausi et al., 2017')
        genes2011 = readPrevGenes('datafiles/changetal2011_hkgenes.txt')
        genes2013 = readPrevGenes('datafiles/eisenberglevanon2013_hkgenes.txt')
        genes2017 = readPrevGenes('datafiles/caracausietal2017_hkgenes.txt')
        #plotVennDiag(genes2011,genes2013, ('Chang et al., 2011','Eisenberg and Levanon, 2013'), outseed+'_prevresults_venndiag.png')
        plotVennDiagPrevHK(hkgenelists,(),outseed+'_all3prevresults_venndiag.pdf')
        #('Chang et al., 2011','Eisenberg and Levanon, 2013','Caracausi et al., 2017')

    if len(noisefile) > 0:
        noisedata = readNoiseSimFile(noisefile)
        #statistical testing on noisedata
        pairsdone = []
        statresults = []
        for key1,accvals1 in noisedata.iteritems():
            for key2,accvals2 in noisedata.iteritems():
                if key1==key2: continue
                if (key1,key2) in pairsdone or (key2,key1) in pairsdone: continue
                _,pval = scipy.stats.ks_2samp(accvals1,accvals2)
                statresults.append((key1,key2,pval))
                pairsdone.append((key1,key2))
        for test in statresults:
            print(test)
        plotNoiseSimResults(noisedata,outseed+'_noisesimulation_barplot.pdf')

        #noisepvals = {}
        #for f in glob.glob('outputs/simdata/noisetests/*_1_allpvals.txt'):
        #    fparts = f.split('_')
        #    noiselabel = fparts[-3][10:]
        #    noisepvals[noiselabel],_ = readPValFile(f)
        #print(noisepvals.keys())
        #plotPvalHist(noisepvals,outseed+'_pvalues_histogram.png')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-f',type=str,default='',help='Salmon quant file of expression values')
    parser.add_argument('-g',type=str,default='',help='File containing list of top n genes')
    parser.add_argument('-p',nargs='+',default=[],help='true pval file')
    parser.add_argument('-t',type=str,default='',help='thresh files')
    parser.add_argument('-n',type=str,default='',help='null model pval files')
    parser.add_argument('-go',nargs='+',default=[],help='Files with GO enrichment results on HBM and GTEx')
    parser.add_argument('-pr',type=str,default='',help='File with info on previous HK results')
    parser.add_argument('-nf',type=str,default='',help='File containing noise results')

    parser.add_argument('-o',type=str,help='Name of output seed for figure')

    args = parser.parse_args()
    main(args.f, args.g, args.p, args.t, args.n, args.go, args.pr, args.nf, args.o)
