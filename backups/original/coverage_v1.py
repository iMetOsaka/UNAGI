import os, sys
import pandas as pd
def main(argv):
        inputPath=argv[0]        
        inf = open(inputPath,"r")
        lists = [[] for _ in range(1,18)]
        #chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17 = ([] for i in range(1,18))
        k1= 'NC_001133.9'
        k2=0
        for l in inf:   
            L = l.strip().split("\t")
            if L[0] == k1:
                lists[k2].append(L)
            else:
                    k2+=1
                    k1 = L[0]         
                    lists[k2].append(L)

            
        genes_6={}
            
        for chromo in lists:
            coverage = []
            for sublist in chromo:
                    coverage.append(int(sublist[1]))
                
    
            curr = 0
            range_end = len(chromo)-2
            for k in range(0,range_end):
                if coverage[k] > 6 and k > curr:
                    pos = k 
                    for i in range(pos,range_end):
                        if coverage[i]>5:
                            continue
                        else:
                            curr = i-1
                            start = chromo[0][0] + "," + str(pos)
                            print(start,curr)
                            break
#initialisation
if __name__ == "__main__":
	main(sys.argv[1:])
                   