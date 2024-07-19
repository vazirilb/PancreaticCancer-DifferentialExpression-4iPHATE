#Pipeline used for RNASeq analysis

#QC (trim_galore)
trim_galore --phred33 --paired -q 28 --retain_unpaired -stringency 5 -e 0.1 -r1 35 -r2 35 file1.fqgz file2.fqgz

#Alignment (STAR)
STAR --runMode genomeGenerate --genomeDir STAR_GRCh38.p12_OH149.genome --genomeFastaFiles GRCh38.p12.genome.fa --sjdbGTFfile gencode.v30.basic.annotation.gtf --sjdbOverhang 149

STAR --genomeDir STAR_GRCh38.p12_OH149.genome --runThreadN 12 --readFilesIn trimmed1.fqgz trimmed2.fqgz --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --outSAMunmapped Within

#Transcript quantification (Salmon)
gffread -w gencode.v30.transcripts.salmon.fa -g GRCh38.p12.genome.fa gencode.v30.basic.annotation.gtf
salmon quant -t gencode.v30.transcripts.salmon.fa --gencode -l A --gcBias -a aligned_file.bam


#load into R, annotation             
library(tximport)
library(edgeR)           
library(biomaRt)
txdb=loadDb(file=gencode.v30.basic.annotation.TxDb)
k=keys(txdb,keytype="TXNAME")
tx2gene=select(txdb,k,"GENEID","TXNAME")
txi=tximport(file1_quant.sf,type="salmon",tx2gene=tx2gene,ignoreTxVersion=TRUE)
cts=round(txi$counts,0)

ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
uni.genes=unique(gsub("\\..*","",rownames(txi$abundance)))
anno.df=getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id","gene_biotype")
              ,filters="ensembl_gene_id",values=uni.genes,mart=ensembl)

dg.list=DGEList(counts=cts                         
                ,genes=anno.df[match(gsub("\\..*","",rownames(cts)),anno.df$ensembl_gene_id),])

dg.list=dg.list[rownames(subset(dg.list$genes,gene_biotype=="protein_coding")),]
dg.list=dg.list[rownames(subset(dg.list$genes,!chromosome_name%in%c("MT","Y"))),]                             
keep.ix=apply(dg.list$counts[,dg.list$samples$serum=="High"],1,FUN=function(x){sum(x >= 20)>1})
keep.ix=keep.ix & apply(dg.list$counts,1,sum)>100 & apply(dg.list$counts,1,max)>50
dg.list=dg.list[keep.ix,,keep.lib.sizes=FALSE]

#output counts matrix
write.csv(dg.list$counts)
