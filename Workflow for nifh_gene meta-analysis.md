## This is the workflow for nifh_gene meta-analysis

### Use Usearch to deal with the high-squencing data

``` 
### The fastq_mergepairs command merges (assembles) paired-end reads to create consensus sequences and, optionally, consensus quality scores

usearch -fastq_mergepairs *_R1.fastq -relabel @ -fastq_maxdiffs 8 -fastqout merged.fq

###Cut barcode xbp +  xbp in left and  xbp in right (according to the primer length)
usearch -fastx_truncate merged.fq -stripleft 29 -stripright 18 -fastqout stripped.fq

###Performs quality filtering 
usearch -fastq_filter filtered1.fa -fastq_trunclen 405 -fastq_maxee 0.8 -fastaout filtered.fa

### Find the set of unique sequences in an input file, also called dereplication.
usearch -fastx_uniques filtered.fa -fastaout uniques.fa -sizeout -relabel Uniq

###Uses the UNOISE algorithm to perform denoising (error-correction) of amplicon reads and creat Zotu (similar with the ASV)
usearch -unoise3 uniques.fa -zotus zotus.fa

###Genenate the OTU table (The table require further treatment)
usearch -otutab merged.fq -otus zotus.fa -otutabout otutab.txt -mapout map.txt

```
### Use FrameBot to translate nucleic acid sequences into protein sequences and the sequences encoding proteins that did not match the nifH protein sequence or that contained termination codons were discarded.
that contained termination codons were discarded.
```
/usr/bin/java -jar /project/framebot/RDPTools-2.0.3/FrameBot.jar framebot -o nifH_test /project/ChenZ/Database/nifH_test.index zotus.fa
```

### Generate new OTU table
```
seqtk seq -U nifH_test_corr_prot.fasta > 11.fasta

usearch -closed_ref 11.fasta -db /project/ChenZ/Database/22nifh_modified/all_nifh.faa -otutabout otutab_new.txt -id 0.97 -tabbedout closed.txt 

awk '{print $1 "\t" $4}' closed.txt > test1.txt

sed 's/;tax.*$//' test1.txt > test2.txt

awk '{print$2}' test2.txt > picked.txt

seqkit grep -f picked.txt /project/ChenZ/Database/22nifh_modified/all_nifh.faa > new_otus.fa

tsv-utils annotation otutab.txt test2.txt > replace.txt

cat replace.txt | sed '/\*/'d > replace_del.txt

awk 'NR==1{print}' otutab.txt > id.txt

sed  "s/#OTU ID/#OTU\tOTU2/g" id.txt > id2.txt

cat id2.txt replace_del.txt > new1.txt

awk '{$1="";print $0}' OFS="\t" new1.txt > new2.txt

cat new2.txt | sed 's/^[ \t]\+//' | sed '/\-/'d > new_otutab.txt
```

### Use R software to merge the all Otu Tables 
```{r}
library(plyr)
a = list.files("input2")            ##input2 
dir = paste("./input2/",a,sep="")
n = length(dir)

### input all table files

for (i in 1:n){
   spe<-read.table(file=dir[i],header=T,check.names=FALSE,sep="\t")
   spe=aggregate(spe[,-1],list(spe$OTU2),sum)
   aa=dir[i]
   bb=as.numeric(nchar(aa))-4
   cc=substring(aa,10,bb)                
   write.csv(spe,file=paste("./input/result_",cc,".csv",sep=""),row.names=TRUE)     
 }


a = list.files("input")
dir = paste("./input/",a,sep="")
n = length(dir)
merge.data = read.csv(file = dir[1],header=TRUE,sep=",")
merge.data=merge.data[,-1]
for (i in 2:n){
  new.data = read.csv(file = dir[i], header=TRUE,sep=",")
  new.data=new.data[,-1]
  merge.data = merge(merge.data,new.data,by="Group.1",all=TRUE)
}

aa=as.data.frame(merge.data)
aa[is.na(aa)] <- 0
bb=as.data.frame.matrix(aa)
bb=arrange(bb,bb[,1])                      
cc=bb[,1]

## Output the merged OTU tables and the represent sequences

write.csv(bb,file = "./input/merge2.csv",row.names=F) 

write.table(cc,"./input/pick_last.txt",row.names=F,col.names=F,quote=FALSE)
```

