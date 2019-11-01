package main

import(
	"flag"
	"fmt"
	"strings"
	"os"
	"os/exec"
	"bufio"
	"encoding/csv"
	"strconv"
	"log"
	"github.com/gonum/stat"
	"sort"
	"math"
	"path/filepath"
	"runtime"
//	"encoding/binary"
)


func main() {

	runtime.GOMAXPROCS(1)

	inputfiles := flag.String("i","","list of input files, separated by commas or one input file containing all data as columns")
	inputtype := flag.String("t","","type of input - f if file containing list of filenames, s if filestring for glob")
	outfile := flag.String("o","","filename for outputs")

	flag.Parse()
	inputlist := strings.Split(*inputfiles,",")
	
	//read in all files from list, or use glob if filestring is input
	var filelist []string
	if len(inputlist) == 1 && *inputtype == "s" {
		var err error
		filelist,err = filepath.Glob(inputlist[0])
		if err != nil { fmt.Println(err) }
	} else if len(inputlist) == 1 && *inputtype == "f" {
		// read file containing filelist
		file, err := os.Open(inputlist[0])
		if err != nil {
			fmt.Print(err)
			os.Exit(1)
		}
		defer file.Close()
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			filelist = append(filelist, scanner.Text())
		}
	} else {
		filelist = inputlist
	}
	numcond := len(filelist)
	var namelist []string
	var datalists [][]float64
	if numcond < 2 {
		namelist,datalists = readCombinedFile(filelist[0]) 
	} else {
		datalists = make([][]float64,numcond)
		for i,filename := range filelist {
			namelist,datalists[i] = readSalmonFiles(filename)
		}
	}

	//remove transcripts with all zero values - too many pointless hypotheses tested
	datalists_nnz := make([][]float64,len(datalists))
	namelist_nnz := make([]string, len(namelist))
	idx := 0
	for i := 0; i < len(namelist); i++ {
		allzero := true
		for i2 := 0; i2 < len(datalists); i2++ {
			if datalists[i2][i] > 0.0 { allzero = false }
			if len(datalists_nnz[i2]) == 0 { datalists_nnz[i2] = make([]float64, len(datalists[i2])) }
		}
		if !allzero {
			namelist_nnz[idx] = namelist[i]
			for i3 := 0; i3 < len(datalists); i3++ {
				datalists_nnz[i3][idx] = datalists[i3][i]
			}
			idx++
		}
	}
	namelist_nnz = namelist_nnz[:idx]
	for i:=0; i < len(datalists); i++ {
		datalists_nnz[i] = datalists_nnz[i][:idx]
	}
	//fmt.Println("features with nonzeros:",namelist_nnz,len(namelist),len(namelist_nnz))

	//convert raw input values to percentages
	datalists_nnz = convertToPercentage(datalists_nnz)

	//write these to file
	//writeDataToFile(filelist,namelist_nnz,datalists_nnz,"geneexp_as_percentages.txt")
	//os.Exit(1)

	//fmt.Println("as pct:",datalists)
	//compute variance of each feature
	variances := computeFeatureVariance(datalists_nnz)
	//fmt.Println("variances:",variances[:10])

	//compute importance of each feature (median?)
	medianranks := computeFeatureImportance(datalists_nnz)
	//fmt.Println("medians:",medianranks[:10])
	
	chisqvals := combinePvals(variances,medianranks)

	pvals := chisqToPval(chisqvals)

	// multiple hypthesis correction!
	bhidx,pvalsidx := multHypTestBH(pvals)

	sortedvar := make([]float64,len(variances))
	sortedmed := make([]float64,len(medianranks))
	//sortedpvals := make([]float64,len(pvals))
	sortednames := make([]string,len(namelist_nnz))

	for oidx,sidx := range pvalsidx {
		sortedvar[oidx] = variances[sidx]
		sortedmed[oidx] = medianranks[sidx]
		//sortedpvals[oidx] = pvals[sidx]
		sortednames[oidx] = namelist_nnz[sidx]
	}

	writeOutputToFile(sortednames[:bhidx+1], pvals[:bhidx+1], *outfile) 
	//writeOutputToFile(sortednames, pvals, *outfile) 
	//fmt.Println(bhidx)
}


func readCombinedFile(filename string) ([]string, [][]float64) {

	var labels []string
	var data [][]float64

	if file,err := os.Open(filename); err == nil {
		defer file.Close()

		//scanner := bufio.NewScanner(file)
		//if scanner.Scan() {scanner.Text() } // skips header line
		scanner := csv.NewReader(file)
		scanner.Comma = '\t'
		readfile,_ := scanner.ReadAll()
		for _,line := range readfile {
			//line := strings.Fields(scanner.Text())
			labels = append(labels,string(line[0]))
			// if data is empty, create data to have (len(line)-1) slices
			// for each entry in line, convert to float and append to appropriate slice
			if len(data) == 0 {
				data = make([][]float64,len(line)-1)
			}
			for idx,val := range line[1:] {
				floatval,interr := strconv.ParseFloat(val,64)
				if interr != nil { log.Fatal(err) }
				data[idx] = append(data[idx],floatval)
			}
		}
	} else {
		fmt.Println("Could not open file, please check filename again")
		os.Exit(1)
	}
	return labels,data
}

func readSalmonFiles(filename string) ([]string, []float64) {

	var txnames []string
	var tpmvals []float64

	if file, err := os.Open(filename); err == nil {
		defer file.Close()
		
		scanner := bufio.NewScanner(file)
		if scanner.Scan() { scanner.Text() } // skip first line (header)
		for scanner.Scan() {
			line := strings.Fields(scanner.Text())
			txnames = append(txnames,line[0])
			tpm,interr := strconv.ParseFloat(line[3],64)
			if interr != nil { log.Fatal(err) }
			tpmvals = append(tpmvals,tpm)
		}
	} else {
		fmt.Println("Could not open file, please check filename again")
		os.Exit(1)
	}
	return txnames,tpmvals
}

func convertToPercentage(datalists [][]float64) [][]float64 {

	//numfeatures := len(datalists[0])
	numconditions := len(datalists)
	perclist := make([][]float64,numconditions)
	for idx := 0; idx < numconditions; idx++ {
		perclist[idx] = valsToPct(datalists[idx])
	}
	return perclist
}

type Slice struct {
	sort.Float64Slice
	idx []int
}


func valsToPct(vals []float64) []float64 {

	pcts := make([]float64, len(vals))	
	orig_vals := make([]float64, len(vals))
	copy(orig_vals,vals)
	valstosort := &Slice{Float64Slice: sort.Float64Slice(vals), idx:make([]int,len(vals))}
	for i := range valstosort.idx {
		valstosort.idx[i] = i
	}
	sort.Sort(valstosort)

	for i := 0; i < len(vals); i++ {
		if i > 0 && valstosort.Float64Slice[i] == valstosort.Float64Slice[i-1] {
			// if value is repeated, use same (max) percentage for all instances
			pcts[valstosort.idx[i]] = pcts[valstosort.idx[i-1]]
		} else {
			pcts[valstosort.idx[i]] = float64(len(vals) - i)/float64(len(vals))
		}
	}
	return pcts
}

func (s Slice) Swap(i, j int) {
	s.Float64Slice.Swap(i, j)
	s.idx[i], s.idx[j] = s.idx[j], s.idx[i]
}


func computeFeatureVariance(datalists [][]float64) []float64 {

	numfeatures := len(datalists[0])
	numconditions := len(datalists)
	featslice := make([]float64,numconditions)
	variances := make([]float64,numfeatures)
	// for each feature, make a slice of the feature values across conditions, compute variance of slice
	for idx := 0; idx < numfeatures; idx++ {
		for condnum := 0; condnum < numconditions; condnum++ {
			featslice[condnum] = datalists[condnum][idx]
		}
		variances[idx] = stat.Variance(featslice,nil)
	}
	//os.Exit(1)

	//redo using prev method
	varranks := valsToPct(variances)
	minvar := 1.0
	for i:= 0; i < len(varranks); i++ {
		varranks[i] = 1.0 - varranks[i]
		if varranks[i] > 0 && varranks[i] < minvar { minvar = varranks[i] }
	}
	// remove zeros
	minvar = minvar/2.0
	for i := 0; i < len(varranks); i++ {
		if varranks[i] == 0 { varranks[i] = minvar }
	}
	return varranks
}

func computeFeatureImportance(datalists [][]float64) []float64 {

	numfeatures := len(datalists[0])
	numconditions := len(datalists)
	featslice := make([]float64,numconditions)
	medians := make([]float64,numfeatures)
	// for each feature, make a slice of the feature values across conditions, compute variance of slice
	for idx := 0; idx < numfeatures; idx++ {
		for condnum := 0; condnum < numconditions; condnum++ {
			featslice[condnum] = datalists[condnum][idx]
		}
		sort.Float64s(featslice)
		if len(featslice)%2 != 0 {
			medians[idx] = featslice[len(featslice)/2]
		} else {
			medians[idx] = (featslice[len(featslice)/2] + featslice[len(featslice)/2 + 1])/2.0
		}
	}
	//os.Exit(1)
	//fmt.Println(medians)
	medranks := valsToPct(medians)
	minval := 1.0
	for i:= 0; i < len(medranks); i++ {
		medranks[i] = 1.0 - medranks[i]
		if medranks[i] > 0 && medranks[i] < minval { minval = medranks[i] }
	}
	// remove zeros
	minval = minval/2.0
	for i := 0; i < len(medranks); i++ {
		if medranks[i] == 0 { medranks[i] = minval }
	}
	//fmt.Println(medranks)
	return medranks
}

func multHypTestBH(pvals []float64) (int,[]int) {
        // multiple hypothesis correction through Benjamini-Hochberg procedure, at level 0.05

        alpha := 0.05
        numtests := len(pvals)
        // order lists by p-values
        //sort.Slice(allbdyvis, func(i,j int) bool {return allbdyvis[i].pval < allbdyvis[j].pval})
	//orig_vals := make([]float64, len(vals))
	//copy(orig_vals,vals)
	pvalssort := &Slice{Float64Slice: sort.Float64Slice(pvals), idx:make([]int,len(pvals))}
	for i := range pvalssort.idx {
		pvalssort.idx[i] = i
	}
	sort.Sort(pvalssort)

	imax := -1
        //sort.Float64s(pvals)
        // define thresholds
        thresh := make([]float64, numtests)
        for i := 0; i < numtests; i++ {
                thresh[i] = float64(i)*alpha/float64(numtests)
        }
        // find greatest index where pval < thresh
        
        //fmt.Println(len(pvals), len(thresh))
        for i := numtests-1; i >= 0; i-- {
                //fmt.Println(i)
                if pvalssort.Float64Slice[i] < thresh[i] {
                        imax = i
                        break
                }
        }
        return imax,pvalssort.idx
}

func writeOutputToFile(txnames []string, pvals []float64, outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
	labelline := []string{"label","p-value"}
        //fmt.Println(strings.Join(labelline, "\t"))
        line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,line1+"\n")

	for idx,vals := range txnames {
                strvals := make([]string, 2)
                strvals[0] = vals
                strvals[1] = strconv.FormatFloat(pvals[idx],'g',-1,64)
                newline := strings.Join(strvals, "\t")
                fmt.Fprintf(w,newline+"\n")
        }
        w.Flush()
        f.Close()
        fmt.Println("Wrote output values to", outfile)
}


func combinePvals(variances []float64, medranks []float64) []float64 {

	chisq := make([]float64, len(variances))

	// compute number of things with both lower val in variances AND lower val in medranks
	/*for i := 0; i < len(pvals); i++ {
		count := 0
		for j := 0; j < len(variances); j++ {
			if variances[j] <= variances[i] && medranks[j] <= medranks[i] { count++ }
		}
		pvals[i] = float64(count)/float64(len(pvals))
	}*/

	// using fisher's method
	for i := 0; i < len(variances); i++ {
		chisq[i] = -2.0*(math.Log(variances[i]) + math.Log(medranks[i]))
		// these are chi-squared values
	}

	return chisq
}

func chisqToPval(chi2vals []float64) []float64 {

	pvals := make([]float64,len(chi2vals))
	for i := 0; i < len(chi2vals); i++{
		chi2 := chi2vals[i]
		pvals[i] = math.Exp(-chi2/2.0) * (1.0 + chi2/2.0)
		//fmt.Println(chi2,pvals[i])
	}
	return pvals
}

func writeDataToFile(filelist []string,namelist_nnz []string,datalists_nnz [][]float64,outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
	labelline := append([]string{"sample"},namelist_nnz...)
	//labelline := []string{alllabels}
        //fmt.Println(strings.Join(labelline, "\t"))
        line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,line1+"\n")

	for idx,samplename := range filelist {
                strvals := make([]string, len(datalists_nnz[0])+1)
                strvals[0] = samplename
		for idx2,dataval := range datalists_nnz[idx] {
			strvals[idx2+1] = strconv.FormatFloat(dataval,'g',-1,64)
		}
                newline := strings.Join(strvals, "\t")
                fmt.Fprintf(w,newline+"\n")
        }
        w.Flush()
        f.Close()
        fmt.Println("Wrote data values to", outfile)
}
