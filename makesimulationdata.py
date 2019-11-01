import csv
import argparse
import numpy as np
import random
import sys
# need to create an mxn matrix

#params needed: m (num features), n (num conditions), value range, "high" median, "low" median, "high" variance, "low" variance, proportion of high/high and low/low, addition of noisy data (proportion, noise level)

# each feature needs to be drawn from a specific distribution with given median and variance. for this, assume all samples have same distribution

def generateMatrix(numfeats,numsamples,highmean,highvar,lowmean,lowvar,weight):
    
    simmat = []
    numhigh = int(weight*numfeats)
    for x in range(numhigh):
        simmat.append(generateRandFeatVals(numsamples,highmean,highvar))
    for x in range(numfeats - numhigh):
        simmat.append(generateRandFeatVals(numsamples,lowmean,lowvar))
    return simmat

def generateRandFeatVals(length,mean,var):
    #note: this values will be in [0,1]
    return [random.gauss(mean,var) for x in range(length)]
    

def generateNoisyFeatVals(numnoisy):

    return [random.uniform(0,1) for x in range(numnoisy)]

def generateMatrixNoisy(numfeats,numsamples,weightnoise,weighttrue,highmean,highvar,lowmean,lowvar):

    simmat = []
    numnoise = int(weightnoise*numsamples)
    numsamples = numsamples- numnoise
    numtrue = int(weighttrue*numfeats)
    print('should be',numtrue,'true high median, low variance features')
    numother = numfeats - numtrue
    #for x in range(numnoise): # this is wrong, don't want noisy features, want some noisy/crappy SAMPLES
    #    simmat.append([random.uniform(0,1) for x in range(numsamples)])
    for x in range(numtrue):
        simmat.append(generateRandFeatVals(numsamples,highmean+0.1,lowvar/2)+generateNoisyFeatVals(numnoise))
    for x in range(numother):
        feattype = random.randint(1,3)
        if feattype == 1:
            simmat.append(generateRandFeatVals(numsamples,highmean,highvar)+generateNoisyFeatVals(numnoise))
        elif feattype == 2:
            simmat.append(generateRandFeatVals(numsamples,lowmean,highvar)+generateNoisyFeatVals(numnoise))
        else:
            simmat.append(generateRandFeatVals(numsamples,lowmean,lowvar)+generateNoisyFeatVals(numnoise))
    return simmat

def writeToFile(matrix,filename):
    
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for idx,featvals in enumerate(matrix):
            fwriter.writerow([idx]+featvals)
    print('wrote simulated results to',filename)

def checkMedAndVar(simmat):

    medians = []
    variances = []
    for feat in simmat:
        medians.append(np.median(feat))
        variances.append(np.var(feat))
    print(medians)
    print(variances)

def main(numfeats,numsamples,highmean,lowmean,highvar,lowvar,weight,noisy,allparams,outfile):

    if allparams:
        sizeparams = [[10000,1000],[20000,10000],[20000,20]]
        meanparams = [[0.85,0.15],[0.7,0.3]]
        varparams = [[0.5,0.1],[0.4,0.2]]
        weightparams = [0.25,0.5,0.75]
        for sizes in sizeparams:
            for means in meanparams:
                for varset in varparams:
                    simmat = generateMatrix(sizes[0],sizes[1],means[0],varset[0],means[1],varset[1],0.5)
                    outfilename = outfile+'_'+str(sizes[0])+'feats_'+str(sizes[1])+'samps_'+str(means[0])+'mean'+str(means[1])+'_'+str(varset[0])+'var'+str(varset[1])+'_weight0.5.txt'
                    writeToFile(simmat,outfilename)

        for weight in weightparams:
            sizes = sizeparams[0]
            means = meanparams[0]
            varset = varparams[0]
            simmat = generateMatrix(sizes[0],sizes[1],means[0],varset[0],means[1],varset[1],0.5)
            outfilename = outfile+'_'+str(sizes[0])+'feats_'+str(sizes[1])+'samps_'+str(means[0])+'mean'+str(means[1])+'_'+str(varset[0])+'var'+str(varset[1])+'_weight'+str(weight)+'.txt'
            writeToFile(simmat,outfilename)
   
    elif noisy == True:
        trueweight = 0.005
        simmat = generateMatrixNoisy(numfeats,numsamples,weight,trueweight ,highmean,highvar,lowmean,lowvar)
        #checkMedAndVar(simmat)
        #sys.exit()
        writeToFile(simmat,outfile)
     
    else:
        print('generating matrix with',numfeats,'features and',numsamples,'samples')
        simmat = generateMatrix(numfeats,numsamples,highmean,highvar,lowmean,lowvar,weight)
        #checkMedAndVar(simmat)
        writeToFile(simmat,outfile)

   

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-nf',type=int,default=10000,help='number of features to simulate')
    parser.add_argument('-ns',type=int,default=1000,help='number of samples to simulate')
    parser.add_argument('-hm',type=float,default=0.85,help='mean for high values')
    parser.add_argument('-lm',type=float,default=0.15,help='mean for low values')
    parser.add_argument('-hv',type=float,default=0.6,help='value for high variances')
    parser.add_argument('-lv',type=float,default=0.1,help='value for low variances')
    parser.add_argument('-w',type=float,default=0.5,help='proportion of high med/high var features or if generating noise, fraction of noisy features')
    parser.add_argument('-n',type=bool,default=False,help='True if generating noisy examples instead of mean/var examples')
    parser.add_argument('-allparams',type=bool,default=False,help='use True to run all parameter settings')
    parser.add_argument('-o',type=str,help='filename for output file')

    args = parser.parse_args()
    main(args.nf,args.ns,args.hm,args.lm,args.hv,args.lv,args.w,args.n,args.allparams,args.o)
