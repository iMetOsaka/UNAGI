import os, sys
def main(argv):    
    path=argv[0]
    path_ann=argv[1] #genes_model 
    inf = open(path,"r")
    three =[]
    for l in inf:
        L = l.strip().split("\t")
        three.append(L)
    #clustering
    a1=[]
    clustered=[]
    for i in range(len(three)-1):
        if int(three[i+1][1])- int(three[i][1]) < 30:
            a1.append(three[i])
        else:
            if a1 == []:
                clustered.append(three[i])
            else:
                a1.append(three[i])
                n = int(len(a1)/2)
                clustered.append(a1[n])
                a1 =[]
    #ends to chromosomes
    lists_3 = [[] for _ in range(1,18)]
    k1= 'NC_001133.9'
    k2=0
    for i in clustered:
                if str(i[0]) == k1:
                    lists_3[k2].append(int(i[1]))
                else:
                        k2+=1
                        k1 = str(i[0])         
                        lists_3[k2].append(int(i[1])) 
    #loading the annotation
    inf = open(path_ann,"r")
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
                        
    #finding peaks and new annotation
    for ch in range(0,17):
        #finding the ends in the middle of genes
        for i in lists[ch]:
            a3 =[]
            a3.append(int(i[1])-30)#bc we will add 30 to all start sites later        
            for k in range(int(i[1]),int(i[2])):
                if k in lists_3[ch] and int(i[2])- k > 100:
                    a3.append(k)
            a3.append(int(i[2]))
            for k in range(0, int(len(a3))-1):            
                print((str(i[0]),a3[k]+30,a3[k+1]))
#initialisation
if __name__ == "__main__":
	main(sys.argv[1:])