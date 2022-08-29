Oyster Genome Population Analysis
================
Jon Puritz
6/1/2022

## Initial Setup

Clone Repo

``` bash
git clone https://github.com/The-Eastern-Oyster-Genome-Project/2021_Genome_Submissions.git
```

Change into repo directory

``` bash
cd 2021_Genome_Submissions
mkdir -p Masked
mkdir -p Population_Genomics
cd Masked
```

Download masked genome

Download raw sequence files

## Run variant calling

Run modified dDocent pipeline **but make sure to change nunber of
processors to match your machine and your email address**

``` bash
bash ../scripts/dDocent_ngs ../other_files/config
```

Wait for code to finsh (will take a few days)

## Filter variant calls **Note** 40 is the number of processors

``` bash
cd raw.vcf
bash ../../../scripts/PVCF_filter.sh MASKED 40 
```

## Structural Variant Calling

``` bash
cd Masked
conda activate HMASK
ls *-RGmd.bam | sed 's/-RGmd.bam//g'| parallel -j 40 "delly call -x repeats.bed -o {}.bcf -g reference.fasta {}-RGmd.bam"
delly merge -o sites.bcf *.bcf
ls *-RGmd.bam | sed 's/-RGmd.bam//g'| parallel "delly call -x repeats.bed -o {}.geno.bcf -g reference.fasta -v sites.bcf {}-RGmd.bam"
bcftools merge --threads 40 -m id -O b -o merged.bcf *_*.geno.bcf

bcftools view --threads 40 -S ../other_files/filter.set merged.bcf -o no.inbred.sv.bcf -O b
bcftools index --threads 40 no.inbred.sv.bcf

delly filter -f germline -o germline.sv.bcf no.inbred.sv.bcf

bcftools view --threads 40 germline.sv.bcf -i 'FILTER == "PASS" & SVTYPE == "INV" ' | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' | mawk '{print $1 "\t" $2-1 "\t" $3-1 "\t" $4}' > delly.inversions.masked.bed

bedtools merge -i delly.inversions.masked.bed > delly.merged.inversions.masked.bed

bcftools view -f "PASS" -i 'F_MISSING==0'  germline.sv.bcf |bcftools +setGT -- -t q -n . -i 'FORMAT/FT="LowQual"'  | bcftools +fill-tags| bcftools view -i 'F_MISSING==0' |bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\n' | mawk '{print $1 "\t" $2-1 "\t" $3-1 "\t" $4 "\t" $5 "\t" $3 - $2}' > sv.full.nomissing.bed

bcftools view -f "PASS" -i 'F_MISSING==0'  germline.sv.bcf | bcftools +setGT -- -t q -n . -i 'FORMAT/FT="LowQual"' | bcftools +fill-tags| bcftools view -i 'F_MISSING==0' | bcftools query -f '%ID\n' > filtered.sv.sites

bcftools view --threads 40 -i "ID=@filtered.sv.sites" merged.bcf -O b -o filtered.sv.bcf

bcftools view -f "PASS" -i 'SVTYPE == "DEL" || SVTYPE == "DUP" ' filtered.sv.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t[%RDCN\t]\n' > masked.cnv.tab

bcftools view -f "PASS" -i 'SVTYPE == "DEL" || SVTYPE == "DUP" ' filtered.sv.bcf | bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' | mawk '{print $1 "\t" $2-1 "\t" $3-1 "\t" $4}' > delly.cnv.masked.bed



cp delly.inversions.masked.bed delly.merged.inversions.masked.bed masked.cnv.tab delly.cnv.masked.bed sv.full.bed sv.full.nomissing.bed ../Population_Genomics



mkdir -p SV_calls

mv *.bcf ./SV_calls
cd ..
#
#conda install dicey
#dicey chop reference.fasta
#bwa mem -t 16 reference.fasta read1.fq.gz read2.fq.gz | samtools sort -@ 8 -o srt.bam -
#samtools index srt.bam
#dicey mappability2 srt.bam 
#gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 
#
#ls *-RGmd.bam | sed 's/-RGmd.bam//g'| parallel "delly call -x repeats.bed -o {}.cnv.bcf -m map.fa.gz -g reference.fasta -v sites.bcf {}-RGmd.bam"
#
#delly merge -e -p -o sites.cnv.bcf -m 1000 -n 100000 *.cnv.bcf
#
#ls *-RGmd.bam | sed 's/-RGmd.bam//g'| parallel -j 16 "delly cnv -u -v sites.cnv.bcf -o {}.cnv.geno.bcf -m map.fa.gz -g reference.fasta -l sites.bcf #{}-RGmd.bam"
#bcftools merge -m id -O b -o merged.cnv.geno.bcf *.cnv.geno.bcf
#delly classify -o filtered.cnv.bcf -f germline -q 30 -x 0.25 merged.cnv.geno.bcf
```

## Setup Analyses

``` bash
cd Population_Genomics
ln -s ../Masked/raw.vcf/filtered/*.vcf.gz .

mkdir -p ./Output/Delly
cp delly.inversions.masked.bed delly.merged.inversions.masked.bed masked.cnv.tab delly.cnv.masked.bed sv.full.bed sv.full.nomissing.bed ./Output/Delly
```

### Setup R Environment, Colors, and plotting set

``` r
library("dplyr")
library("ggplot2")
library(R.utils)
library(gghighlight)
library(ggman)
library(ggtext)
library(patchwork)
library(plotrix)
library(plyr)
library(qqman)
library(qvalue)
library(reshape2)
library(tidyr)
library(zoo)
library(infer)
options(dplyr.summarise.inform = FALSE)
library(bigsnpr)
library("wesanderson")
library("directlabels")
library(OutFLANK)
library(pcadapt)
library(adegenet)
library(poppr)
library(vcfR)
library(maps)
library(mapdata)
library(fiftystater)
library("rgdal") 
library("plyr")
library(RColorBrewer)
library(matrixStats)
library(data.table)
library(ggdendro)
library(ggridges)
library(hierfstat)
library(stringr)
```

### Setup

``` r
pops <- read.table("./other_files/population_coordinates.full.tab", header = TRUE)
popind <- read.table("./other_files/popmap", header=TRUE)

popmap <-join(pops,popind,type = "inner")
```

    ## Joining by: POP

``` r
popmap <- popmap[order(popmap$IND),]
popmap$Type <- factor(popmap$Type, levels=c("Inbred","Selected","Wild"))

popmap$POP <- factor(popmap$POP, levels=c("HI","SM","UMFS","NEH","NG", "HG","HC","CS","DEBY" ,"CLP","HC-VA","LOLA","CL","SL","OBOYS","LM"))

popmap.ni <- droplevels(subset(popmap, POP != "NG" & POP !="HG" & POP != "LM"))

p_noinbred <- c("#FF7F00","#FDBF6F","#fda56f","#871a1a","#E31A1c","#FB9A99","#8a584a" ,"#33A02C","#B2DF8A","#86a395","#1F78B4","#A6CEE3","#A6CEE3")

full_c <- c("#FF7F00","#FDBF6F","#fda56f","#871a1a","#707070","Grey","#E31A1c","#FB9A99","#8a584a" ,"#33A02C","#B2DF8A","#86a395","#1F78B4","#A6CEE3","#A6CEE3","#6A3D9A")

p1 <- c("#1F78B4","#33A02C","#FB9A99","#FB9A99","#E31A1c","#B2DF8A","Grey","#FF7F00","#B2Df8A","#E31A1C","#707070","#A6CEE3","#A6CEE3","#FDBF6F","#FDBF6F")

popmap.wild <- droplevels(subset(popmap.ni, POP != "DEBY" & POP !="LOLA" & POP != "NEH" & POP != "OBOYS" & POP != "UMFS"))

p_wild <- c("#FF7F00","#FDBF6F","#E31A1c","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3")
p_wild_map <- c("#1F78B4","#33A02C", "#FB9A99", "#E31A1c","#B2DF8A","#FF7F00","#A6CEE3","#FDBF6F" )

popmap.wildAE <- droplevels(subset(popmap.wild, POP != "CL" & POP !="SL" ))
p_wildAE <- c("#FF7F00","#FDBF6F","#E31A1c","#FB9A99","#33A02C","#B2DF8A")

gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.ni$IND
gp.popmap$pop <- popmap.ni$POP

target <- c("CL", "SL", "OBOYS","LOLA","DEBY","HC-VA","CLP","CS","HC","NEH","UMFS","SM","HI")
target <- rev(target)
gp.popmap <-gp.popmap %>% arrange(factor(pop, levels = target))


# Points to plot for each population and the population name
pop_coordinates <- pops
pop_coordinates$Lat <- as.character(pop_coordinates$Lat)
pop_coordinates$Lat <- as.numeric(pop_coordinates$Lat)
pop_coordinates$Long <- as.character(pop_coordinates$Long)
pop_coordinates$Long <- as.numeric(pop_coordinates$Long)

pop_coordinates <- pop_coordinates[order(pop_coordinates$POP),]

pop_coordinates.D <- subset(pop_coordinates, POP == "LOLA" | POP== "NEH" | POP =="OBOYS" | POP =="DEBY" | POP =="UMFS")

pop_coordinates.I <- subset(pop_coordinates, POP == "NG" | POP =="HG")

prefix.I <- "./Population_Genomics/Intermediate_Files/"
prefix.O <- "./Population_Genomics/Output/"
prefix.OE <- "./Population_Genomics/Output/Extra/"
prefix.OF <- "./Population_Genomics/Output/Figures/"
prefix.OFS <- "./Population_Genomics/Output/Figures/Supplemental/"
prefix <- "./Population_Genomics/"
NCORES <- nb_cores()

pal1 <- wes_palette("Zissou1")
pal2 <- wes_palette("GrandBudapest1")
pal3 <- wes_palette("GrandBudapest2")
pal4 <- wes_palette("Royal2")
col_pal <- c(pal2[1],pal1[1], pal2[3],pal3[2],pal3[1],pal4[5],pal3[3],pal3[4],pal4[1],pal4[3])
chroms <- c("1","2","3","4","5","6","7","8","9","10")
ncbis <- c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1","NC_035786.1","NC_035787.1","NC_035788.1","NC_035789.1")
```

# Sequencing and Mapping Statistics

``` bash
cd Masked
for i in `cat namelist` 
do 
READS=$(mawk '/total reads:/' ./logfiles/$i.trim.log | head -1 | cut -f3 -d" ")
READS=$(($READS * 2))
FILREADS=$(mawk '/reads passed filter/' ./logfiles/$i.trim.log | cut -f4 -d " ")
paste <(echo -e $i "\t" $READS "\t" $FILREADS) <(samtools view -@64 -F 0x100 -F 0x800 -F0x4 -c $i-RGmd.bam) <(samtools view -@64 -f 0x400 -c $i-RGmd.bam) <(samtools view -@64 -F 0x100 -F 0x800 -c $i.F.bam); done > sequencing.mapping.counts

cat <(echo -e "SAMPLE\tRAW_READS\tTRIMMED_READS\tTOTAL_MAPPINGS\tDUPLICATES_DETECTED\tMAPPINGS_POST_FILTERING\tPOP") <(paste sequencing.mapping.counts <(cut -f2 ../other_files/popmap)) > sequencing.mapping.counts.txt

mv sequencing.mapping.counts sequencing.mapping.counts.txt ../../Population_Genomics


cat bamlist.list | parallel -j 8 "samtools coverage {} > {}.coverage"
rm seq.mean.coverage
for i in *.bam.coverage; do mawk '!/NC_007175.2/ && !/#/' $i | mawk -v x=${i%.F.bam.coverage} '{sum=sum+$6} END {print sum/NR "\t"x}'; done >> seq.mean.coverage
mv seq.mean.coverage ../Population_Genomics
ls *-RGmd.bam | parallel -j 8 "samtools coverage {} > {}.coverage.full"
rm seq.mean.coverage.full
for i in *_*.coverage.full; do mawk '!/NC_007175.2/ && !/#/' $i | mawk -v x=${i%.coverage.full} '{sum=sum+$6} END {print sum/NR "\t"x}'; done >> seq.mean.coverage.full
mv seq.mean.coverage.full ../Population_Genomics
mkdir -p Coverage
mv *.coverage *.coverage.full ./Coverage/
cd ../Population_Genomics
paste <(mawk '{new=$1*684741128/(684741128-106558067); print new "\t" $2}' seq.mean.coverage | cut -f1) <(mawk '{new=$1*684741128/(684741128-106558067); print new "\t" $2}' seq.mean.coverage.full) <(cut -f2 ../other_files/popmap) > total.seq.coverage
cat <(echo -e "Mean_Filtered_Coverage\tMean_Unfiltered_Coverage\tSample\tPOP")  total.seq.coverage > seq.mean.coverage.txt
mkdir -p ./Output/Coverage_Statistics/
cp seq.mean.coverage.txt sequencing.mapping.counts.txt ./Output/Coverage_Statistics/
```

``` r
seq.map.results <-read.table(paste(prefix.O,"Coverage_Statistics/","sequencing.mapping.counts.txt", sep=""), sep="\t", header=T)

seq.table <-rbind(seq.map.results[,2:6] %>% dplyr::summarise(dplyr::across(everything(), mean)), seq.map.results[,2:6] %>% dplyr::summarise(dplyr::across(everything(), std.error)))
row.names(seq.table) <- c("Mean","Stanard Error")
seq.table
```

    ##               RAW_READS TRIMMED_READS TOTAL_MAPPINGS DUPLICATES_DETECTED
    ## Mean           79071288      76812419       74318475           4954759.9
    ## Stanard Error   2036262       1968176        1904303            132522.5
    ##               MAPPINGS_POST_FILTERING
    ## Mean                         54201810
    ## Stanard Error                 1368984

``` r
seq.cov.results <-read.table(paste(prefix,"seq.mean.coverage.txt", sep=""), sep="\t", header=T)
cov.table <- rbind(seq.cov.results[,1:2] %>% dplyr::summarise(dplyr::across(everything(), mean)), seq.cov.results[,1:2] %>% dplyr::summarise(dplyr::across(everything(), std.error)))

row.names(cov.table) <- c("Mean","Stanard Error")
cov.table
```

    ##               Mean_Filtered_Coverage Mean_Unfiltered_Coverage
    ## Mean                      80.4515689               87.4132300
    ## Stanard Error              0.3228894                0.2277883

``` r
group_by(seq.cov.results, POP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Mean_Filtered_Coverage, na.rm = TRUE),
    sd = sd(Mean_Filtered_Coverage, na.rm = TRUE),
    se=std.error(Mean_Filtered_Coverage, na.rm = TRUE) )
```

    ## # A tibble: 16 × 5
    ##    POP   count  mean    sd    se
    ##    <fct> <int> <dbl> <dbl> <dbl>
    ##  1 CL        6  78.0 1.39  0.566
    ##  2 CLP       6  81.2 0.960 0.392
    ##  3 CS        6  80.4 0.656 0.268
    ##  4 DEBY      6  79.9 0.700 0.286
    ##  5 HC        6  80.8 1.78  0.728
    ##  6 HC-VA     6  81.0 0.865 0.353
    ##  7 HG        3  94.4 3.61  2.09 
    ##  8 HI        6  81.3 1.33  0.542
    ##  9 LM        5  78.2 2.30  1.03 
    ## 10 LOLA      6  79.3 1.06  0.432
    ## 11 NEH       6  81.5 0.897 0.366
    ## 12 NG        4  79.9 1.40  0.701
    ## 13 OBOYS     6  78.6 1.11  0.455
    ## 14 SL        6  79.0 1.10  0.447
    ## 15 SM        6  80.0 0.514 0.210
    ## 16 UMFS      6  80.1 0.901 0.368

``` r
group_by(seq.cov.results, POP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Mean_Unfiltered_Coverage, na.rm = TRUE),
    sd = sd(Mean_Unfiltered_Coverage, na.rm = TRUE),
    se=std.error(Mean_Unfiltered_Coverage, na.rm = TRUE) )
```

    ## # A tibble: 16 × 5
    ##    POP   count  mean    sd    se
    ##    <fct> <int> <dbl> <dbl> <dbl>
    ##  1 CL        6  85.9 1.12  0.459
    ##  2 CLP       6  87.9 0.745 0.304
    ##  3 CS        6  87.2 0.472 0.193
    ##  4 DEBY      6  87.1 0.518 0.212
    ##  5 HC        6  87.7 1.40  0.570
    ##  6 HC-VA     6  87.8 0.683 0.279
    ##  7 HG        3  97.3 2.14  1.23 
    ##  8 HI        6  87.9 0.947 0.387
    ##  9 LM        5  86.0 1.80  0.804
    ## 10 LOLA      6  86.6 0.849 0.347
    ## 11 NEH       6  88.0 0.691 0.282
    ## 12 NG        4  86.6 1.21  0.603
    ## 13 OBOYS     6  86.3 0.930 0.380
    ## 14 SL        6  86.5 0.857 0.350
    ## 15 SM        6  87.1 0.358 0.146
    ## 16 UMFS      6  87.1 0.688 0.281

# Map Creation

``` r
state <- map_data("state", boundary = FALSE)
w2hr <- map_data("world2Hires","USA")
states <- map_data("state")
cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#999999","#F0E442" , "#CC79A7")
pal <- wes_palette("Zissou1", 8, type = "continuous")

coastline = readOGR(dsn="./other_files/map_data", layer="us_medium_shoreline")
```

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "/home/jpuritz/2022_Population_Genomics/other_files/map_data", layer: "us_medium_shoreline"
    ## with 59023 features
    ## It has 8 fields

``` r
states= readOGR(dsn="./other_files/map_data", layer="state_bounds")
```

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "/home/jpuritz/2022_Population_Genomics/other_files/map_data", layer: "state_bounds"
    ## with 19 features
    ## It has 2 fields
    ## Integer64 fields read as strings:  ID

``` r
coast <- spTransform(coastline, CRS("+proj=longlat +ellps=clrk66"))
states <- spTransform(states, CRS("+proj=longlat +ellps=clrk66"))
coast@data$id = rownames(coastline@data)
states@data$id = rownames(states@data)
coastline.points = fortify(coastline, region="id")
states.points = fortify(states, region="id")
```

``` r
coast.df = join(coastline.points, coast@data, by="id")
states.df = join(states.points, coast@data, by="id")
```

``` r
png(filename=paste(prefix.OE,"PopulationMap.png", sep=""), type="cairo",units="px", width=2500, height=3500, res=200, bg="transparent")
ggplot()+
geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black" ,size=0.25, alpha=0.5) +
geom_path(data = states.df, aes(x=long, y = lat, group = group), color="black",size=0.25, alpha=0.5)+
  geom_point(data=pop_coordinates[which(pop_coordinates$POP != "LM"),], aes(x=Long,y=Lat, fill=POP), size=8,color="grey", shape=21,alpha=0.9)+
  scale_fill_manual(values=p1, name="Population")+
  scale_color_manual(values=p1, name="Population")+
  geom_label_repel(data=pop_coordinates[which(pop_coordinates$POP != "LM"),], aes(x=Long, y=Lat, label=POP, color=POP), fill="#f2f2f2", min.segment.length = 0,segment.colour = 'grey40',size=12, segment.size=0.5, box.padding=1,fontface="bold", alpha=1) + xlab("Longitude") + ylab("Latitude") + guides(fill="none", color="none")+
  xlab("Longitude") + 
  ylab("Latitude") +
  #coord_equal() +
  coord_map(xlim = c(-95.5, -67.5),  ylim = c(24.25, 47.25), projection = "mercator")+
  theme_bw() +
  theme(text = element_text(size=24),plot.margin = margin(0, 0, 0, 0, "cm"))
```

    ## Warning: Removed 7 rows containing missing values (geom_point).

    ## Warning: Removed 7 rows containing missing values (geom_label_repel).

``` r
dev.off()
```

    ## png 
    ##   2

``` r
finished_mapnolm <-ggplot()+
geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black" ,size=0.25, alpha=0.5) +
geom_path(data = states.df, aes(x=long, y = lat, group = group), color="black",size=0.25, alpha=0.5)+
  geom_point(data=pop_coordinates[which(pop_coordinates$POP != "LM" & pop_coordinates$Lat != "NA"),], aes(x=Long,y=Lat, fill=POP), size=5,color="grey", shape=21,alpha=0.9)+
  scale_fill_manual(values=p_wild_map, name="Population")+
  scale_color_manual(values=p_wild_map, name="Population")+
  geom_label_repel(data=pop_coordinates[which(pop_coordinates$POP != "LM" & pop_coordinates$Lat != "NA"),], aes(x=Long, y=Lat, label=POP, color=POP), fill="#f2f2f2", min.segment.length = 0,segment.colour = 'grey40',size=8, segment.size=0.5, box.padding=1,fontface="bold", alpha=1) + xlab("Longitude") + ylab("Latitude") + guides(fill="none", color="none")+
  xlab("Longitude") + 
  ylab("Latitude") +
  coord_map(xlim = c(-97.5, -65),  ylim = c(22.75, 48.5), projection = "mercator")+
  theme_bw() +
  theme(text = element_text(size=12))

finished_mapnolm_small <-ggplot()+
geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black" ,size=0.25, alpha=0.5) +
geom_path(data = states.df, aes(x=long, y = lat, group = group), color="black",size=0.25, alpha=0.5)+
  geom_point(data=pop_coordinates[which(pop_coordinates$POP != "LM" & pop_coordinates$Lat != "NA"),], aes(x=Long,y=Lat, fill=POP), size=5,color="grey", shape=21,alpha=0.9)+
  scale_fill_manual(values=p_wild_map, name="Population")+
  scale_color_manual(values=p_wild_map, name="Population")+
  geom_label_repel(data=pop_coordinates[which(pop_coordinates$POP != "LM" & pop_coordinates$Lat != "NA"),], aes(x=Long, y=Lat, label=POP, color=POP), fill="#f2f2f2", min.segment.length = 0,segment.colour = 'grey40',size=3.5, segment.size=0.5, box.padding=1,fontface="bold", alpha=1) + xlab("Longitude") + ylab("Latitude") + guides(fill="none", color="none")+
  xlab("Longitude") + 
  ylab("Latitude") +
  coord_map(xlim = c(-97.5, -65),  ylim = c(22.75, 48.5), projection = "mercator")+
  theme_bw() +
  theme(text = element_text(size=12))
```

``` r
wild_ae.coor <- droplevels(pop_coordinates[which(pop_coordinates != "LM" & pop_coordinates != "SL" & pop_coordinates != "CL" & pop_coordinates != "OBOYS" & pop_coordinates != "LOLA" & pop_coordinates != "NEH" & pop_coordinates != "DEBY" &  pop_coordinates != "UMFS" &  pop_coordinates != "NG" &  pop_coordinates != "HG"),])

wild_ae.coor <- wild_ae.coor[complete.cases(wild_ae.coor),]

p_wildAE.color <- c( "#33A02C","#FB9A99","#E31A1c", "#B2DF8A", "#FF7F00", "#FDBF6F")


finished_AE <-ggplot()+
#geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#CCCC99E6", color="black" ,size=0.25) +
geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black" ,size=0.25, alpha=0.5) +
geom_path(data = states.df, aes(x=long, y = lat, group = group), color="black",size=0.25, alpha=0.5)+
  geom_point(data=wild_ae.coor, aes(x=Long,y=Lat, fill=POP), size=5,color="grey", shape=21,alpha=0.9)+
  scale_fill_manual(values=p_wildAE.color, name="Population")+
  scale_color_manual(values=p_wildAE.color, name="Population")+
  geom_label_repel(data=wild_ae.coor, aes(x=Long, y=Lat, label=POP, color=POP), min.segment.length = 0,segment.colour = 'grey40',size=3.5, segment.size=0.5, box.padding=1, fill="grey95", fontface = "bold", alpha=1) + 
  xlab("Longitude") + ylab("Latitude")+ guides(fill="none", color="none")+
  coord_map(xlim = c(-85, -65),  ylim = c(37, 47), projection = "mercator")+
  theme_bw() +
  theme(text = element_text(size=12),plot.margin = margin(0, 0, 0, 0, "cm"))
```

# Linkage

``` bash
cd Population_Genomics

VCF='SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz'

RESDIR='./Intermediate_Files/LD/AtlanticEastern/'

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 499990 --ld-window-bp 500000 --out $RESDIR'geno_ld_window_499950-500000_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 99990 --ld-window-bp 100000 --out $RESDIR'geno_ld_window_99990-100000_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 49990 --ld-window-bp 50000 --out $RESDIR'geno_ld_window_49990-50000_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 9990 --ld-window-bp 10000 --out $RESDIR'geno_ld_window_9990-10000_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 4990 --ld-window-bp 5000 --out $RESDIR'geno_ld_window_4990-5000_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 490 --ld-window-bp 500 --out $RESDIR'geno_ld_window_490-500_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 190 --ld-window-bp 200 --out $RESDIR'geno_ld_window_190-200_AtlanticEastern' &

vcftools --gzvcf $VCF --keep ../other_files/wildAE.pop --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out $RESDIR'geno_ld_window_50-50_AtlanticEastern' 

RESDIR='./Intermediate_Files/LD/selected/'
vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 499950 --ld-window-bp 500000 --out $RESDIR'geno_ld_window_499950-500000_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 99990 --ld-window-bp 100000 --out $RESDIR'geno_ld_window_99990-100000_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 49990 --ld-window-bp 50000 --out $RESDIR'geno_ld_window_49990-50000_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 9990 --ld-window-bp 10000 --out $RESDIR'geno_ld_window_9990-10000_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 4990 --ld-window-bp 5000 --out $RESDIR'geno_ld_window_4990-5000_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 490 --ld-window-bp 500 --out $RESDIR'geno_ld_window_490-500_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 190 --ld-window-bp 200 --out $RESDIR'geno_ld_window_190-200_selected' &

vcftools --gzvcf $VCF --keep ../other_files/domestic.pop --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out $RESDIR'geno_ld_window_50-50_selected' 
```

``` r
### Code by KEL and BDW and modified slightly by JBP

load.ld.data <- function (path = NULL) {
    if (!dir.exists(path)) {
        stop("ERROR: The path that you have provided is not a directory")
    } else {
        # Full paths from working directory to file
        file.paths <- paste(path, list.files(path), sep = "")
    }
    # Just file names from path
    file.names <- unlist(lapply(strsplit(file.paths, "/"), function(x) return(x[5])))
    
    # Initialize array that will hold all window sizes
    windows <- c()
    # Initialize list that will accumulate data.frames that will be named
    ld.vars <- list()
    for (i in 1:length(file.paths))
    {
        # current.window - range of bp window provided in the file names
        current.window <- strsplit(file.names[i], split = "_")[[1]][4]
        current.df <- read.csv(file.paths[which(grepl(paste("_", current.window, sep=""), file.paths, fixed=TRUE) == TRUE)], 
                               sep="\t",
                               header=TRUE,
                               stringsAsFactors=FALSE)
        # window.name - minimum window size in the file name to use as name of data.frame
        window.name <- strsplit(current.window, split = "-")[[1]][1]
        # Add current.window to array of all windows
        windows <- c(windows, current.window)
        # Add current.df to list of all named data.frames
        ld.vars[[window.name]] <- current.df
    }
    
    return(ld.vars)
}

# get the LD decay plots for all of the chromosomes in the ld.datasets

decay.plots <- function (path = NULL, file.suffix = NULL, plot.title = NULL, file.prefix= NULL) {
    ld.datasets <- load.ld.data(path = path)

    CHR <- c("NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", 
             "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", 
             "NC_035788.1", "NC_035789.1")
    
    col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                  '#ffffff', '#000000')
    
    chr.num <- length(unique(ld.datasets[[1]]$CHR))
    for (i in length(ld.datasets)) {
        if (length(unique(ld.datasets[[i]]$CHR)) != chr.num) {
            print(names(ld.datasets[i]))
            #stop("The number of chromosomes in the datasets do not match")
        } else {

        }
    }
    
    chr.index <- rep(NA, chr.num)
    for (i in 1:chr.num){
        chr.index[i] <- substr(unique(ld.datasets[[1]]$CHR)[i], start = 9, stop = 9)
    }
    chr.index <- sort(chr.index)

    ld.stats <- data.frame(matrix(NA, length(ld.datasets), chr.num))
    ld.error <- data.frame(matrix(NA, length(ld.datasets), chr.num))
    ld.stats[,1] <- as.numeric(names(ld.datasets))
    ld.error[,1] <- as.numeric(names(ld.datasets))
    for (i in 1:length(ld.datasets)) {
        ld.datasets[[i]]$R.2 <- as.numeric(ld.datasets[[i]]$R.2)
        for (j in 1:chr.num) {
            ld.stats[i, j+1] <- summary(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])[4]
            ld.error[i, j+1] <- sd(na.omit(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])) / sqrt(length(na.omit(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])))
        }
    }
    
    ld.stats <- ld.stats[order(ld.stats[, 1]), ]
    ld.error <- ld.error[order(ld.error[, 1]), ]

    png(paste(prefix.OFS, file.prefix, file.suffix, ".png", sep=""), type="cairo",
        width = 8,
        height = 8,
        units = "in",
        res = 300)
    par(bty="l")
    for (i in 2:ncol(ld.stats)) {
        if (i == 2){
            plot(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                 pch = i,
                 col = col_vector[i],
                 cex = 2,
                 xlab="Distance (bp)", 
                 ylab=bquote("LD"~(R^2)~"\U00B1 Std. Error"),
                 ylim=c(0, 0.3),
                 type = "b",
                 log = "x",
                 bty="l")
                 arrows(ld.stats[,1], ld.stats[,i]-ld.error[,i], ld.stats[,1], ld.stats[,i]+ld.error[,i], col=col_vector[i], length=0.05, angle=90, code=3)
        } else {
            points(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                   pch = i,
                   col = col_vector[i],
                   cex = 2,
                   type = "b")
            arrows(ld.stats[,1], ld.stats[,i]-ld.error[,i], ld.stats[,1], ld.stats[,i]+ld.error[,i], col=col_vector[i], length=0.05, angle=90, code=3)
        }
        legend("topright", legend = CHR, pch = 2:(length(CHR)+1), col = col_vector[2:(length(CHR)+1)], box.lty=0)
    }
    dev.off()
}
```

``` r
# Get plot for selection Atlantic populations
 decay.plots(path="Population_Genomics/Intermediate_Files/LD/AtlanticEastern/",
             file.suffix = "AtlanticEastern", file.prefix= "Figure.S3A.")
```

    ## png 
    ##   2

``` r
 decay.plots(path="Population_Genomics/Intermediate_Files/LD/selected/",
             file.suffix = "selected", file.prefix= "Figure.S3B.")
```

    ## png 
    ##   2

# PCA

``` bash
bcftools view  --threads 40 -S ../other_files/filter.set ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz |  bcftools +fill-tags| bcftools view -i 'F_MISSING==0 & MAF > 0.05'  -O z -o ./Population_Genomics/VCF/SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz

plink2 --vcf ./Population_Genomics/SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz --make-bed --out ./Population_Genomics/PLINK/Masked.NoInbred --allow-extra-chr
```

``` r
rds.m.ni <- snp_readBed2(paste(prefix,"Masked.NoInbred.bed", sep=""), backingfile = tempfile(),ncores = NCORES)
obj.bed <- bed(paste(prefix,"Masked.NoInbred.bed", sep=""))
```

``` r
masked_all.m.ni<- snp_attach(rds.m.ni)
G.m.ni <- masked_all.m.ni$genotypes
CHR.m.ni <- masked_all.m.ni$map$chromosome
POS.m.ni <- masked_all.m.ni$map$physical.pos
```

``` r
svd.m.ni <- snp_autoSVD(G.m.ni, as.integer(factor(CHR.m.ni)), POS.m.ni, ncores = NCORES,min.mac = 7, size=10)
save(svd.m.ni,  file = "LD_thinned_NI_data.RData")
```

Load FST Data if needed

``` r
# To load the data again
load("LD_thinned_NI_data.RData")
```

``` r
ind.keep <- attr(svd.m.ni, "subset")

rssq <- bigsnpr:::prod_and_rowSumsSq(
  obj.bed, ind_row = rows_along(obj.bed), ind_col = ind.keep,
  center = svd.m.ni$center, scale = svd.m.ni$scale, V = matrix(0, length(ind.keep), 0))[[2]]

sum(rssq)  # squared frobenius norm of genotype matrix for subset (= total variance explained)
```

    ## [1] 84971119

``` r
var_exp <- svd.m.ni$d^2 / sum(rssq)
signif(var_exp, 2)
```

    ##  [1] 0.072 0.030 0.026 0.024 0.022 0.020 0.017 0.017 0.016 0.014

``` r
round(cumsum(var_exp), 3)
```

    ##  [1] 0.072 0.103 0.128 0.152 0.174 0.194 0.211 0.227 0.243 0.258

``` r
dffp <- as.data.frame(svd.m.ni$u[,1:6])
colnames(dffp) <- c("PC1","PC2","PC3","PC4","PC5","PC6")
dffp$POP <- popmap.ni$POP
dffp$Type <- popmap.ni$Type
```

``` r
png(filename=paste(prefix.OE,"PCA.allpops.1.2.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_noI_12<- ggplot(data = dffp,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6)+
scale_shape_manual(values=c(23,21),name="Origin") +
scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") + #ggtitle("PCA of all Populaitons/Lines") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method = list(cex = 2, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_noI_12
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OE,"PCA.allpops.3.4.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_noI_34<- ggplot(data = dffp,aes( x=PC3, y=PC4,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6)+
scale_shape_manual(values=c(23,21),name="Origin") +
scale_fill_manual(values=p_noinbred) + scale_color_manual(values=p_noinbred) + 
#guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =6)))+
guides(fill = "none", alpha="none",size="none" , shape= "none", color ="none")+
geom_dl(alpha=1,method = list(cex = 2, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_noI_34
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OFS,"Figure.S1.PCAs.png", sep=""), type="cairo",units="px", width=6000, height=3000, res=200, bg="transparent")
S1A <- pca_noI_12 + inset_element(finished_mapnolm, left = 0.5, bottom = 0.5, right = 0.975, top = 0.975, align_to = 'full')
S1A  + pca_noI_34 + plot_layout(guides="collect")
dev.off()
```

    ## png 
    ##   2

#### Setup FST Calculations later

``` r
filename <- paste(prefix.I,"NoInbredLociLocationsAfterThinning.maskedPCs.txt", sep="")
write(paste(CHR.m.ni[attr(svd.m.ni, "subset")],
            POS.m.ni[attr(svd.m.ni, "subset")], sep='\t'),filename)
```

### Wild only samples

``` r
rds.m.w <- snp_readBed2(paste(prefix,"Masked.Wild.bed", sep=""), backingfile = tempfile(),ncores = NCORES)
obj.bed <- bed(paste(prefix,"Masked.Wild.bed", sep=""))
```

``` r
masked_all.m.w <- snp_attach(rds.m.w)
G.m.w <- masked_all.m.w$genotypes
CHR.m.w <- masked_all.m.w$map$chromosome
POS.m.w <- masked_all.m.w$map$physical.pos
```

``` r
svd.m.w <- snp_autoSVD(G.m.w, as.integer(factor(CHR.m.w)), POS.m.w, ncores = NCORES, min.mac =4,size=10)

### Save data as calculations take time
save(svd.m.w,  file = "LD_thinned_w_data.RData")
```

Load FST Data if needed

``` r
# To load the data again
load("LD_thinned_w_data.RData")
```

``` r
dffp <- as.data.frame(svd.m.w$u[,1:6])
colnames(dffp) <- c("PC1","PC2","PC3","PC4","PC5","PC6")
dffp$POP <- popmap.wild$POP
dffp$Type <- popmap.wild$Type
```

``` r
png(filename=paste(prefix.OE,"PCA.wild.1.2.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_wild_12<- ggplot(data = dffp,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wild, name="Population/Line") + scale_color_manual(values=p_wild , name="Population/Line") + ggtitle("PCA of all Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wild_12
dev.off()
```

    ## png 
    ##   2

``` r
ind.keep <- attr(svd.m.w, "subset")

rssq <- bigsnpr:::prod_and_rowSumsSq(
  obj.bed, ind_row = rows_along(obj.bed), ind_col = ind.keep,
  center = svd.m.w$center, scale = svd.m.w$scale, V = matrix(0, length(ind.keep), 0))[[2]]

sum(rssq)  # squared frobenius norm of genotype matrix for subset (= total variance explained)
```

    ## [1] 43425602

``` r
var_exp <- svd.m.w$d^2 / sum(rssq)
signif(var_exp, 2)
```

    ##  [1] 0.077 0.028 0.026 0.023 0.022 0.021 0.021 0.021 0.021 0.021

``` r
round(cumsum(var_exp), 3)
```

    ##  [1] 0.077 0.104 0.130 0.153 0.175 0.197 0.218 0.239 0.260 0.281

``` r
png(filename=paste(prefix.OE,"PCA.wild.3.4.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_wild_34<- ggplot(data = dffp,aes( x=PC3, y=PC4,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wild, name="Population/Line") + scale_color_manual(values=p_wild , name="Population/Line") + ggtitle("PCA of all Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wild_34
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OE,"PCA.wild.5.6.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_wild_56<- ggplot(data = dffp,aes( x=PC5, y=PC6,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wild, name="Population/Line") + scale_color_manual(values=p_wild , name="Population/Line") + ggtitle("PCA of all Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wild_56
dev.off()
```

    ## png 
    ##   2

``` r
filename <- paste(prefix.I,"WildLociLocationsAfterThinning.maskedPCs.txt", sep="")
write(paste(CHR.m.w[attr(svd.m.w, "subset")],
            POS.m.w[attr(svd.m.w, "subset")], sep='\t'),filename)
```

``` bash
bcftools view --threads 40 -S ./other_files/wildAE.pop ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz | bcftools +fill-tags|bcftools view -i 'F_MISSING==0 & MAF > 0.05' -O z -o ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz

bcftools index --threads 40 ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz

plink2 --vcf ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz --make-bed --out ./Population_Genomics/MASKED.wildAE --allow-extra-chr
```

``` r
rds.m.wae <- snp_readBed2(paste(prefix,"MASKED.wildAE.bed", sep=""), backingfile = tempfile(),ncores = NCORES)
obj.bed <- bed(paste(prefix,"MASKED.wildAE.bed", sep=""))
```

``` r
masked_all.m.wae <- snp_attach(rds.m.wae)
G.m.wae <- masked_all.m.wae$genotypes
CHR.m.wae <- masked_all.m.wae$map$chromosome
POS.m.wae <- masked_all.m.wae$map$physical.pos
```

``` r
svd.m.wildae <- snp_autoSVD(G.m.wae, as.integer(factor(CHR.m.wae)), POS.m.wae, ncores = NCORES, min.mac =4,size=10)
#### Save data as calculations take time
save(svd.m.wildae,  file = "LD_thinned_wae_data.RData")
```

Load FST Data if needed

``` r
#### To load the data again
load("LD_thinned_wae_data.RData")
```

``` r
filename <- paste(prefix.I,"WildAELociLocationsAfterThinningmaskedPCs.txt", sep="")
write(paste(CHR.m.wae[attr(svd.m.wildae, "subset")],
            POS.m.wae[attr(svd.m.wildae, "subset")], sep='\t'),
      filename)
```

``` r
dffp <- as.data.frame(svd.m.wildae$u[,1:6])
colnames(dffp) <- c("PC1","PC2","PC3","PC4","PC5","PC6")
dffp$POP <- popmap.wildAE$POP
dffp$Type <- popmap.wildAE$Type
```

``` r
png(filename=paste(prefix.OE,"PCA.wildAE.1.2.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")

pca_wildAE_12<- ggplot(data = dffp,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wildAE, name="Population/Line") + scale_color_manual(values=p_wildAE , name="Population/Line") + ggtitle("PCA of all Atlantic Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none", shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wildAE_12

dev.off()
```

    ## png 
    ##   2

``` r
shape_n <- "Thinned Genome Wide SNPs"
shapes <- c( "Thinned Genome Wide SNPs" = 21 )


dffp$Type <- factor(shape_n)[1]
pca_wildAE_12p <- ggplot(data = dffp,aes( x=PC1, y=PC2,fill=POP,color=POP, shape=Type, size=3, alpha=0.75)) + 
geom_point(alpha=0.75, size = 4, color="black")+
scale_fill_manual(values=p_wildAE, name="Population") + 
scale_color_manual(values=p_wildAE , name="Population") +
scale_shape_manual(values=shapes , name="Genomic Architecture") +
xlab("PC1")+ ylab("PC2") +
guides(alpha="none",size="none", fill = guide_legend(override.aes = list(size =4, shape = c(21,21,21,21,21,21),alpha = 0.75, order =1)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=12)) + theme(legend.position = "right")




ind.keep <- attr(svd.m.wildae, "subset")

rssq <- bigsnpr:::prod_and_rowSumsSq(
  obj.bed, ind_row = rows_along(obj.bed), ind_col = ind.keep,
  center = svd.m.wildae$center, scale = svd.m.wildae$scale, V = matrix(0, length(ind.keep), 0))[[2]]

sum(rssq)  # squared frobenius norm of genotype matrix for subset (= total variance explained)
```

    ## [1] 22448860

``` r
var_exp <- svd.m.wildae$d^2 / sum(rssq)
signif(var_exp, 2)
```

    ##  [1] 0.039 0.037 0.032 0.032 0.030 0.029 0.029 0.029 0.029 0.029

``` r
round(cumsum(var_exp), 3)
```

    ##  [1] 0.039 0.076 0.108 0.140 0.170 0.199 0.228 0.256 0.285 0.313

``` r
png(filename=paste(prefix.OE,"PCA.wildAE.3.4.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_wildAE_34<- ggplot(data = dffp,aes( x=PC3, y=PC4,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wildAE, name="Population/Line") + scale_color_manual(values=p_wildAE , name="Population/Line") + ggtitle("PCA of all Atlantic Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none"  , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wildAE_34
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OE,"PCA.wildAE.5.6.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_wildAE_56<- ggplot(data = dffp,aes( x=PC5, y=PC6,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6, shape = 21)+
scale_fill_manual(values=p_wildAE, name="Population/Line") + scale_color_manual(values=p_wildAE , name="Population/Line") + ggtitle("PCA of all Atlantic Wild Populaitons") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,21,21,21,21),alpha = 0.95)), alpha="none",size="none"  , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))

pca_wildAE_56
dev.off()
```

    ## png 
    ##   2

# FST and Outlier Calculations

``` r
rds.comp.m.random <- snp_readBed("./Population_Genomics/Random.50k.masked.bed", backingfile = tempfile())
rds.wild.m.random <- snp_readBed("./Population_Genomics/Random.50k.masked.wild.bed", backingfile = tempfile())
rds.wildAE.m.random <- snp_readBed("./Population_Genomics/Random.50k.Masked.WildAE.bed", backingfile = tempfile())

masked_comp.m.r <- snp_attach(rds.comp.m.random)
G.comp.m.r <- masked_comp.m.r$genotypes
CHR.comp.m.r <- masked_comp.m.r$map$chromosome
POS.comp.m.r <- masked_comp.m.r$map$physical.pos

masked_wild.m.r <- snp_attach(rds.wild.m.random)
G.wild.m.r <- masked_wild.m.r$genotypes
CHR.wild.m.r <- masked_wild.m.r$map$chromosome
POS.wild.m.r <- masked_wild.m.r$map$physical.pos

masked_wildAE.m.r <- snp_attach(rds.wildAE.m.random)
G.wildAE.m.r <- masked_wildAE.m.r$genotypes
CHR.wildAE.m.r <- masked_wildAE.m.r$map$chromosome
POS.wildAE.m.r <- masked_wildAE.m.r$map$physical.pos
```

``` r
my_fst.comp.m <- MakeDiploidFSTMat(G.m.ni[], locusNames = paste0(CHR.m.ni,"_", POS.m.ni), popNames = popmap.ni$POP)

my_fst.comp.m.r <- MakeDiploidFSTMat(G.comp.m.r[], locusNames = paste0(CHR.comp.m.r,"_", POS.comp.m.r), popNames = popmap.ni$POP)

my_fst.wild.m <- MakeDiploidFSTMat(G.m.w[], locusNames = paste0(CHR.m.w,"_", POS.m.w), popNames = popmap.wild$POP)

my_fst.wild.m.r <- MakeDiploidFSTMat(G.wild.m.r[], locusNames = paste0(CHR.wild.m.r,"_", POS.wild.m.r), popNames = popmap.wild$POP)

my_fst.wildae.m <- MakeDiploidFSTMat(G.m.wae[], locusNames = paste0(CHR.m.wae,"_", POS.m.wae), popNames = popmap.wildAE$POP)

my_fst.wildae.m.r <- MakeDiploidFSTMat(G.wildAE.m.r[], locusNames = paste0(CHR.wildAE.m.r,"_", POS.wildAE.m.r), popNames = popmap.wildAE$POP)

save(my_fst.comp.m, my_fst.comp.m.r, my_fst.wildae.m, my_fst.wildae.m.r,my_fst.wild.m,my_fst.wild.m.r,  file = "FST_data.RData")
```

Load FST Data if needed

``` r
# To load the data again
load("FST_data.RData")
```

``` r
#### Evaluating OutFLANK with trimmed SNPs ####
out_trim.m <- OutFLANK(my_fst.comp.m.r, NumberOfSamples=13, qthreshold = 0.05, Hmin = 0.1)
str(out_trim.m)
```

    ## List of 6
    ##  $ FSTbar               : num 0.0793
    ##  $ FSTNoCorrbar         : num 0.154
    ##  $ dfInferred           : num 9.4
    ##  $ numberLowFstOutliers : int 0
    ##  $ numberHighFstOutliers: int 2
    ##  $ results              :'data.frame':   50000 obs. of  15 variables:
    ##   ..$ LocusName       : Factor w/ 50000 levels "NC_035780.1_10002843",..: 2118 3590 3606 4711 4769 4779 4851 4970 82 461 ...
    ##   ..$ He              : num [1:50000] 0.142 0.1841 0.0973 0.1311 0.1311 ...
    ##   ..$ FST             : num [1:50000] 0.29756 0.06207 0.00509 0.11325 0.03096 ...
    ##   ..$ T1              : num [1:50000] 0.021724 0.005769 0.000249 0.007532 0.002048 ...
    ##   ..$ T2              : num [1:50000] 0.073 0.0929 0.049 0.0665 0.0662 ...
    ##   ..$ FSTNoCorr       : num [1:50000] 0.3561 0.1264 0.0836 0.1888 0.1198 ...
    ##   ..$ T1NoCorr        : num [1:50000] 0.026 0.01175 0.0041 0.01255 0.00792 ...
    ##   ..$ T2NoCorr        : num [1:50000] 0.073 0.0929 0.049 0.0665 0.0662 ...
    ##   ..$ meanAlleleFreq  : num [1:50000] 0.923 0.897 0.949 0.929 0.929 ...
    ##   ..$ indexOrder      : int [1:50000] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ GoodH           : Factor w/ 2 levels "goodH","lowH": 1 1 2 1 1 1 1 1 1 1 ...
    ##   ..$ qvalues         : num [1:50000] 0.736 0.907 NA 0.907 0.907 ...
    ##   ..$ pvalues         : num [1:50000] 0.024 0.799 NA 0.541 0.717 ...
    ##   ..$ pvaluesRightTail: num [1:50000] 0.012 0.601 NA 0.271 0.642 ...
    ##   ..$ OutlierFlag     : logi [1:50000] FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
head(out_trim.m$results)
```

    ##           LocusName        He          FST            T1         T2  FSTNoCorr
    ## 1 NC_035780.1_33577 0.1420118  0.297560976  0.0217236467 0.07300570 0.35609756
    ## 2 NC_035780.1_50755 0.1840894  0.062068966  0.0057692308 0.09294872 0.12643678
    ## 3 NC_035780.1_51155 0.0973044  0.005090909  0.0002492877 0.04896724 0.08363636
    ## 4 NC_035780.1_71645 0.1310815  0.113253012  0.0075320513 0.06650641 0.18875502
    ## 5 NC_035780.1_77993 0.1310815  0.030955585  0.0020477208 0.06615028 0.11978466
    ## 6 NC_035780.1_79381 0.2422748 -0.033846154 -0.0041132479 0.12152778 0.04175824
    ##      T1NoCorr   T2NoCorr meanAlleleFreq indexOrder GoodH   qvalues    pvalues
    ## 1 0.025997151 0.07300570      0.9230769          1 goodH 0.7356229 0.02400976
    ## 2 0.011752137 0.09294872      0.8974359          2 goodH 0.9072544 0.79887825
    ## 3 0.004095442 0.04896724      0.9487179          3  lowH        NA         NA
    ## 4 0.012553419 0.06650641      0.9294872          4 goodH 0.9072544 0.54103434
    ## 5 0.007923789 0.06615028      0.9294872          5 goodH 0.9072544 0.71650900
    ## 6 0.005074786 0.12152778      0.8589744          6 goodH 0.9203879 0.03082145
    ##   pvaluesRightTail OutlierFlag
    ## 1       0.01200488       FALSE
    ## 2       0.60056087       FALSE
    ## 3               NA       FALSE
    ## 4       0.27051717       FALSE
    ## 5       0.64174550       FALSE
    ## 6       0.98458927       FALSE

``` r
OutFLANKResultsPlotter(out_trim.m, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.2, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

``` r
hist(out_trim.m$results$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-2.png)<!-- -->

``` r
out_trim.m.w <- OutFLANK(my_fst.wild.m.r, NumberOfSamples=8, qthreshold = 0.05, Hmin = 0.1)
str(out_trim.m.w)
```

    ## List of 6
    ##  $ FSTbar               : num 0.0472
    ##  $ FSTNoCorrbar         : num 0.127
    ##  $ dfInferred           : num 6.5
    ##  $ numberLowFstOutliers : int 0
    ##  $ numberHighFstOutliers: int 11
    ##  $ results              :'data.frame':   50000 obs. of  15 variables:
    ##   ..$ LocusName       : Factor w/ 50000 levels "NC_035780.1_1000349",..: 2163 2868 4764 4866 4915 823 827 926 946 1081 ...
    ##   ..$ He              : num [1:50000] 0.497 0.117 0.264 0.5 0.264 ...
    ##   ..$ FST             : num [1:50000] 0.3642 0.0938 0.0684 0.0282 0.0529 ...
    ##   ..$ T1              : num [1:50000] 0.09549 0.00561 0.00918 0.00714 0.00709 ...
    ##   ..$ T2              : num [1:50000] 0.2622 0.0598 0.1342 0.253 0.1342 ...
    ##   ..$ FSTNoCorr       : num [1:50000] 0.4238 0.1577 0.1396 0.0941 0.1396 ...
    ##   ..$ T1NoCorr        : num [1:50000] 0.11111 0.00942 0.01873 0.02381 0.01873 ...
    ##   ..$ T2NoCorr        : num [1:50000] 0.2622 0.0598 0.1342 0.253 0.1342 ...
    ##   ..$ meanAlleleFreq  : num [1:50000] 0.542 0.938 0.844 0.5 0.844 ...
    ##   ..$ indexOrder      : int [1:50000] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ GoodH           : Factor w/ 2 levels "goodH","lowH": 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ qvalues         : num [1:50000] 0.385 0.987 0.987 0.987 0.987 ...
    ##   ..$ pvalues         : num [1:50000] 0.00394 0.55385 0.7175 0.74876 0.7175 ...
    ##   ..$ pvaluesRightTail: num [1:50000] 0.00197 0.27693 0.35875 0.62562 0.35875 ...
    ##   ..$ OutlierFlag     : logi [1:50000] FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
head(out_trim.m.w$results)
```

    ##            LocusName        He        FST          T1         T2  FSTNoCorr
    ## 1  NC_035780.1_33932 0.4965278 0.36423841 0.095486111 0.26215278 0.42384106
    ## 2  NC_035780.1_42599 0.1171875 0.09377593 0.005605159 0.05977183 0.15767635
    ## 3  NC_035780.1_80203 0.2636719 0.06839187 0.009176587 0.13417659 0.13955638
    ## 4  NC_035780.1_88994 0.5000000 0.02823529 0.007142857 0.25297619 0.09411765
    ## 5  NC_035780.1_93609 0.2636719 0.05286506 0.007093254 0.13417659 0.13955638
    ## 6 NC_035780.1_229531 0.1866319 0.05445026 0.005158730 0.09474206 0.12041885
    ##      T1NoCorr   T2NoCorr meanAlleleFreq indexOrder GoodH   qvalues     pvalues
    ## 1 0.111111111 0.26215278      0.5416667          1 goodH 0.3853534 0.003936431
    ## 2 0.009424603 0.05977183      0.9375000          2 goodH 0.9870541 0.553852524
    ## 3 0.018725198 0.13417659      0.8437500          3 goodH 0.9870541 0.717498420
    ## 4 0.023809524 0.25297619      0.5000000          4 goodH 0.9870541 0.748764643
    ## 5 0.018725198 0.13417659      0.8437500          5 goodH 0.9870541 0.717498420
    ## 6 0.011408730 0.09474206      0.8958333          6 goodH 0.9870541 0.923479064
    ##   pvaluesRightTail OutlierFlag
    ## 1      0.001968216       FALSE
    ## 2      0.276926262       FALSE
    ## 3      0.358749210       FALSE
    ## 4      0.625617679       FALSE
    ## 5      0.358749210       FALSE
    ## 6      0.461739532       FALSE

``` r
OutFLANKResultsPlotter(out_trim.m.w, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.2, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-3.png)<!-- -->

``` r
hist(out_trim.m.w$results$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-4.png)<!-- -->

``` r
out_trim.m.wae <- OutFLANK(my_fst.wildae.m.r, NumberOfSamples=6, qthreshold = 0.05, Hmin = 0.1)
str(out_trim.m.wae)
```

    ## List of 6
    ##  $ FSTbar               : num 0.0152
    ##  $ FSTNoCorrbar         : num 0.0962
    ##  $ dfInferred           : num 6.27
    ##  $ numberLowFstOutliers : int 0
    ##  $ numberHighFstOutliers: int 3
    ##  $ results              :'data.frame':   50000 obs. of  15 variables:
    ##   ..$ LocusName       : Factor w/ 50000 levels "NC_035780.1_10004670",..: 1931 4769 4770 4878 4985 54 436 803 838 839 ...
    ##   ..$ He              : num [1:50000] 0.198 0.33 0.346 0.313 0.129 ...
    ##   ..$ FST             : num [1:50000] 0.0642 -0.0195 0.0709 0.0809 -0.0213 ...
    ##   ..$ T1              : num [1:50000] 0.00648 -0.00324 0.0125 0.01296 -0.00139 ...
    ##   ..$ T2              : num [1:50000] 0.1009 0.1662 0.1764 0.1602 0.0653 ...
    ##   ..$ FSTNoCorr       : num [1:50000] 0.1284 0.046 0.1207 0.1329 0.0603 ...
    ##   ..$ T1NoCorr        : num [1:50000] 0.01296 0.00764 0.0213 0.0213 0.00394 ...
    ##   ..$ T2NoCorr        : num [1:50000] 0.1009 0.1662 0.1764 0.1602 0.0653 ...
    ##   ..$ meanAlleleFreq  : num [1:50000] 0.889 0.208 0.778 0.806 0.931 ...
    ##   ..$ indexOrder      : int [1:50000] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ GoodH           : Factor w/ 2 levels "goodH","lowH": 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ qvalues         : num [1:50000] 0.978 1 0.978 0.978 0.985 ...
    ##   ..$ pvalues         : num [1:50000] 0.47 0.334 0.545 0.429 0.568 ...
    ##   ..$ pvaluesRightTail: num [1:50000] 0.235 0.833 0.273 0.215 0.716 ...
    ##   ..$ OutlierFlag     : logi [1:50000] FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
head(out_trim.m.wae$results)
```

    ##            LocusName        He          FST           T1         T2  FSTNoCorr
    ## 1  NC_035780.1_32455 0.1975309  0.064220183  0.006481481 0.10092593 0.12844037
    ## 2  NC_035780.1_78394 0.3298611 -0.019498607 -0.003240741 0.16620370 0.04596100
    ## 3  NC_035780.1_78556 0.3456790  0.070866142  0.012500000 0.17638889 0.12073491
    ## 4  NC_035780.1_90198 0.3132716  0.080924855  0.012962963 0.16018519 0.13294798
    ## 5  NC_035780.1_98948 0.1292438 -0.021276596 -0.001388889 0.06527778 0.06028369
    ## 6 NC_035780.1_103968 0.1527778 -0.005988024 -0.000462963 0.07731481 0.07185629
    ##      T1NoCorr   T2NoCorr meanAlleleFreq indexOrder GoodH   qvalues   pvalues
    ## 1 0.012962963 0.10092593     0.88888889          1 goodH 0.9781183 0.4695704
    ## 2 0.007638889 0.16620370     0.20833333          2 goodH 0.9996659 0.3344355
    ## 3 0.021296296 0.17638889     0.77777778          3 goodH 0.9781183 0.5454419
    ## 4 0.021296296 0.16018519     0.80555556          4 goodH 0.9781183 0.4292908
    ## 5 0.003935185 0.06527778     0.93055556          5 goodH 0.9852297 0.5678277
    ## 6 0.005555556 0.07731481     0.08333333          6 goodH 0.9781183 0.7654322
    ##   pvaluesRightTail OutlierFlag
    ## 1        0.2347852       FALSE
    ## 2        0.8327823       FALSE
    ## 3        0.2727209       FALSE
    ## 4        0.2146454       FALSE
    ## 5        0.7160861       FALSE
    ## 6        0.6172839       FALSE

``` r
OutFLANKResultsPlotter(out_trim.m.wae, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.2, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-5.png)<!-- -->

``` r
hist(out_trim.m.wae$results$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-53-6.png)<!-- -->

``` r
P1.m <- pOutlierFinderChiSqNoCorr(my_fst.comp.m, Fstbar = out_trim.m$FSTNoCorrbar, 
                                   dfInferred = out_trim.m$dfInferred, qthreshold = 0.0005, Hmin=0.1)

P1.m <- P1.m[order(as.numeric(rownames(P1.m))),,drop=FALSE]                  
hist(P1.m$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

``` r
P1.m.w <- pOutlierFinderChiSqNoCorr(my_fst.wild.m, Fstbar = out_trim.m.w$FSTNoCorrbar, 
                                   dfInferred = out_trim.m.w$dfInferred, qthreshold = 0.0005, Hmin=0.1)

P1.m.w <- P1.m.w[order(as.numeric(rownames(P1.m.w))),,drop=FALSE]                  
hist(P1.m.w$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-54-2.png)<!-- -->

``` r
P1.m.wae <- pOutlierFinderChiSqNoCorr(my_fst.wildae.m, Fstbar = out_trim.m.wae$FSTNoCorrbar, 
                                   dfInferred = out_trim.m.wae$dfInferred, qthreshold = 0.0005, Hmin=0.1)

P1.m.wae <- P1.m.wae[order(as.numeric(rownames(P1.m.wae))),,drop=FALSE]                  
hist(P1.m.wae$pvaluesRightTail)
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-54-3.png)<!-- -->

``` r
m.plot.df<- setNames(data.frame(matrix(ncol = 7, nrow = length(my_fst.comp.m$LocusName))), c("CHR", "CHRR", "BP", "SNP", "P", "Q", "Outlier"))
m.plot.df$CHR <-as.integer(as.factor(CHR.m.ni))
m.plot.df$CHRR <- CHR.m.ni
m.plot.df$SNP <-P1.m$LocusName
m.plot.df$BP <- POS.m.ni
m.plot.df$P <- P1.m$pvalues
m.plot.df$Q <- P1.m$qvalues
m.plot.df$F <- P1.m$FST
m.plot.df$Outlier <- P1.m$OutlierFlag
m.plot.df <- m.plot.df[which(m.plot.df$Outlier == TRUE | m.plot.df$Outlier == FALSE),]

m.w.plot.df<- setNames(data.frame(matrix(ncol = 7, nrow = length(my_fst.wild.m$LocusName))), c("CHR", "CHRR", "BP", "SNP", "P", "Q", "Outlier"))
m.w.plot.df$CHR <-as.integer(as.factor(CHR.m.w))
m.w.plot.df$CHRR <- CHR.m.w
m.w.plot.df$SNP <-P1.m.w$LocusName
m.w.plot.df$BP <- POS.m.w
m.w.plot.df$P <- P1.m.w$pvalues
m.w.plot.df$Q <- P1.m.w$qvalues
m.w.plot.df$F <- P1.m.w$FST
m.w.plot.df$Outlier <- P1.m.w$OutlierFlag
m.w.plot.df <- m.w.plot.df[which(m.w.plot.df$Outlier == TRUE | m.w.plot.df$Outlier == FALSE),]

m.wae.plot.df<- setNames(data.frame(matrix(ncol = 7, nrow = length(my_fst.wildae.m$LocusName))), c("CHR", "CHRR", "BP", "SNP", "P", "Q", "Outlier"))
m.wae.plot.df$CHR <-as.integer(as.factor(CHR.m.wae))
m.wae.plot.df$CHRR <- CHR.m.wae
m.wae.plot.df$SNP <-P1.m.wae$LocusName
m.wae.plot.df$BP <- POS.m.wae
m.wae.plot.df$P <- P1.m.wae$pvalues
m.wae.plot.df$Q <- P1.m.wae$qvalues
m.wae.plot.df$F <- P1.m.wae$FST
m.wae.plot.df$Outlier <- P1.m.wae$OutlierFlag
m.wae.plot.df <- m.wae.plot.df[which(m.wae.plot.df$Outlier == TRUE | m.wae.plot.df$Outlier == FALSE),]
```

``` r
withr::with_options(c(scipen = 10), write(paste(m.plot.df$CHRR,m.plot.df$BP-1,m.plot.df$BP,m.plot.df$F,m.plot.df$Outlier, sep='\t'),
      "./Population_Genomics/Masked_NoInbred_FST_Outflank.bed"))
withr::with_options(c(scipen = 10),write(paste(m.w.plot.df$CHRR,m.w.plot.df$BP-1,m.w.plot.df$BP,m.w.plot.df$F,m.w.plot.df$Outlier, sep='\t'),
      "./Population_Genomics/Masked_Wild_FST_Outflank.bed"))
withr::with_options(c(scipen = 10),write(paste(m.wae.plot.df$CHRR,m.wae.plot.df$BP-1,m.wae.plot.df$BP,m.wae.plot.df$F,m.wae.plot.df$Outlier, sep='\t'),
      "./Population_Genomics/Masked_WildAE_FST_Outflank.bed"))
```

``` bash
bedtools intersect -b <(sort -k 1,1 -k2,2n ./Population_Genomics/Masked_NoInbred_FST_Outflank.bed) -a <(sort -k 1,1 -k2,2n ./other_files/10kb.bed) -wa -wb -sorted  |bedtools merge -d -1 -c 7,8 -o mean,distinct | mawk '{ if ($0 ~ /TRUE/) {print $1 "\t" $2 "\t" $3 "\t" $4 "\tTRUE"} else {print $0}}' > ./Population_Genomics/masked.OF_fst.10kb.bed
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\tOUTLIER") ./Population_Genomics/masked.OF_fst.10kb.bed > ./Population_Genomics/masked.OF_fst.10kb.txt
mkdir -p ./Population_Genomics/Output/Fst
cp ./Population_Genomics/Masked_Wild_FST_Outflank.bed ./Population_Genomics/Masked_WildAE_FST_Outflank.bed ./Population_Genomics/Masked_NoInbred_FST_Outflank.bed ./Population_Genomics/Output/Fst/
```

# Pairwise Population Differentiation

# Wild vs. Selected Line Comparisons

### PI

``` bash
ln -s ./Masked/raw.vcf/filtered/SNP.MASKED.TRSdp5g75.nDNA.g1.maf01.FIL.vcf.gz ./Population_Genomics/

ls ./other_files/*.temp | sed 's/.temp//g' | sed 's:\./other_files/::g' | parallel "bcftools view --threads 4 -S ./other_files/{}.temp ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf01.FIL.vcf.gz | vcftools --vcf -  --window-pi 10000 --out ./Population_Genomics/Intermediate_Files/{}.pi"

bcftools view  --threads 40 -S ./other_files/AE.pop ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf01.max2alleles.FIL.vcf.gz |  bcftools +fill-tags| bcftools view -i 'F_MISSING==0 & MAF > 0.05'  -O z -o ./Population_Genomics/SNP.MASKED.AE.TRSdp5g75.nDNA.g1.maf01.max2alleles.nomissing.FIL.vcf.gz

ls ./other_files/*.temp | sed 's/.temp//g' | sed 's:\./other_files/::g' | parallel "bcftools view --threads 4 -S ./other_files/{}.temp ./Population_Genomics/SNP.MASKED.AE.TRSdp5g75.nDNA.g1.maf01.max2alleles.nomissing.FIL.vcf.gz | vcftools --vcf -  --window-pi 10000 --out ./Population_Genomics/Intermediate_Files/{}.AE.pi"

cd ./Population_Genomics/Intermediate_Files/

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.pi.windowed.pi)  <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.pi.windowed.pi)  > pop.wild.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.AE.pi.windowed.pi)  <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.AE.pi.windowed.pi)  > pop.wild.AE.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS2.pi.windowed.pi)  <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.pi.windowed.pi) > pop.selected.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.AE.pi.windowed.pi)   <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.AE.pi.windowed.pi) > pop.selected.AE.windowed.pi

bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.pi | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > mean.pop.wild.windowed.pi

bedtools merge -i <(mawk '!/CHR/' pop.wild.AE.windowed.pi | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > mean.pop.wild.AE.windowed.pi

bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.pi | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > mean.pop.selected.windowed.pi

bedtools merge -i <(mawk '!/CHR/' pop.selected.AE.windowed.pi | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > mean.pop.selected.AE.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.wild.windowed.pi > ../mean.pop.wild.windowed.pi.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.wild.AE.windowed.pi > ../mean.pop.wild.windowed.AE.pi.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.selected.windowed.pi > ../mean.pop.selected.windowed.pi.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.selected.AE.windowed.pi > ../mean.pop.selected.windowed.AE.pi.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS2.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.pi.windowed.pi) > pop.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.AE.pi.windowed.pi)  <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.AE.pi.windowed.pi) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.AE.pi.windowed.pi) > pop.windowed.AE.pi



cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.wild.windowed.pi mean.pop.selected.windowed.pi > mean.pop.windowed.pi

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tGROUP") mean.pop.wild.AE.windowed.pi mean.pop.selected.AE.windowed.pi > mean.pop.windowed.AE.pi

cd ../
mkdir -p ./Output/pi
cp mean.pop.windowed.AE.pi  mean.pop.windowed.pi mean.pop.domestic.windowed.pi.txt  mean.pop.selected.windowed.AE.pi.txt  mean.pop.selected.windowed.pi.txt  mean.pop.wild.windowed.AE.pi.txt  mean.pop.wild.windowed.pi.txt ./Output/pi
```

``` r
pi.wvd.dataframe<-read.table(paste(prefix.I,"pop.windowed.pi", sep=""), sep="\t", header=T)
pi.wvd.dataframe$GROUP <- factor(pi.wvd.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
pi.wvd.dataframe$POP <- factor(pi.wvd.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

filename <- paste(prefix.I,"pop.windowed.AE.pi", sep="")
pi.wvd.AE.dataframe<-read.table(filename, sep="\t", header=T)
pi.wvd.AE.dataframe$GROUP <- factor(pi.wvd.AE.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
pi.wvd.AE.dataframe$POP <- factor(pi.wvd.AE.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"))

filename <- paste(prefix.I,"mean.pop.windowed.pi", sep="")
pi.wvd.mean.dataframe<-read.table(filename, sep="\t", header=T)
pi.wvd.mean.dataframe$GROUP <- factor(pi.wvd.mean.dataframe$GROUP, levels=c("WILD", "SELECTED" ))

filename <- paste(prefix.I,"mean.pop.windowed.AE.pi", sep="")
pi.wvd.mean.AE.dataframe<-read.table(filename, sep="\t", header=T)
pi.wvd.mean.AE.dataframe$GROUP <- factor(pi.wvd.mean.AE.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
```

``` r
bd.pi <-ggplot(pi.wvd.dataframe, aes(x=PI,y = POP))+
  #geom_violin(aes(color=GROUP,fill= GROUP)) +
  #stat_summary(fun=mean, geom="point", shape=23, size=4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.005)+
  geom_boxplot(aes(fill=GROUP), width=0.75,outlier.shape = 23, outlier.color = "black")+
  geom_vline(xintercept =mean(pi.wvd.dataframe$PI),linetype = "dashed") +
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(pi.wvd.dataframe$POP)))+
  #ylab("Population/Line")+
  labs(x=expression(pi)) +
  ggtitle("Wild vs. Selected Nucleotide Diversity")+
  theme_classic() +
  theme(axis.title.y = element_blank()) 
```

``` r
bc.pi <-ggplot(pi.wvd.mean.dataframe, aes(x=PI,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP)) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = 23, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(pi.wvd.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Mean Wild vs. Selected Nucleotide Diversity")+
  theme_classic() +
  labs(x=expression(pi)) +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank())

bc.pi.ae <-ggplot(pi.wvd.mean.AE.dataframe, aes(x=PI,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP)) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = 23, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(pi.wvd.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Mean Wild vs. Selected Nucleotide Diversity using Atlantic Samples Only")+
  theme_classic() +
  labs(x=expression(pi)) +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank())
```

``` r
pi.wvd.dataframe %>% t_test(formula = PI ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df  p_value alternative estimate  lower_ci upper_ci
    ##       <dbl>   <dbl>    <dbl> <chr>          <dbl>     <dbl>    <dbl>
    ## 1      8.79 401217. 1.45e-18 two.sided   0.000103 0.0000802 0.000126

``` r
pi.wvd.dataframe %>% t_test(formula = PI ~ GROUP, alternative = "greater", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df  p_value alternative estimate  lower_ci upper_ci
    ##       <dbl>   <dbl>    <dbl> <chr>          <dbl>     <dbl>    <dbl>
    ## 1      8.79 401217. 7.25e-19 greater     0.000103 0.0000839      Inf

``` r
group_by(pi.wvd.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PI, na.rm = TRUE),
    sd = sd(PI, na.rm = TRUE),
    se=std.error(PI, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP     count    mean      sd         se
    ##   <fct>     <int>   <dbl>   <dbl>      <dbl>
    ## 1 WILD     299402 0.00538 0.00403 0.00000736
    ## 2 SELECTED 186199 0.00528 0.00394 0.00000914

``` r
group_by(pi.wvd.dataframe, POP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PI, na.rm = TRUE),
    sd = sd(PI, na.rm = TRUE),
    se=std.error(PI, na.rm = TRUE) )
```

    ## # A tibble: 13 × 5
    ##    POP   count    mean      sd        se
    ##    <fct> <int>   <dbl>   <dbl>     <dbl>
    ##  1 HI    37373 0.00524 0.00395 0.0000204
    ##  2 SM    37337 0.00503 0.00386 0.0000200
    ##  3 HC    37733 0.00530 0.00399 0.0000205
    ##  4 CS    37443 0.00533 0.00399 0.0000206
    ##  5 CLP   37397 0.00533 0.00399 0.0000207
    ##  6 HC-VA 37367 0.00532 0.00399 0.0000206
    ##  7 CL    37371 0.00574 0.00419 0.0000217
    ##  8 SL    37381 0.00574 0.00419 0.0000217
    ##  9 UMFS  37216 0.00483 0.00367 0.0000190
    ## 10 NEH   37033 0.00450 0.00346 0.0000180
    ## 11 DEBY  37362 0.00561 0.00409 0.0000212
    ## 12 LOLA  37376 0.00584 0.00421 0.0000218
    ## 13 OBOYS 37212 0.00558 0.00406 0.0000211

``` r
test2<- pi.wvd.dataframe[which(pi.wvd.dataframe$POP != "LOLA" & pi.wvd.dataframe$POP != "DEBY"),]

test2 %>% t_test(formula = PI ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df   p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>     <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      30.2 212317. 1.26e-199 two.sided   0.000406 0.000380 0.000433

``` r
test2 %>% t_test(formula = PI ~ GROUP, alternative = "greater", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df   p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>     <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      30.2 212317. 6.31e-200 greater     0.000406 0.000384      Inf

``` r
group_by(test2, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PI, na.rm = TRUE),
    sd = sd(PI, na.rm = TRUE),
    se=std.error(PI, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP     count    mean      sd         se
    ##   <fct>     <int>   <dbl>   <dbl>      <dbl>
    ## 1 WILD     299402 0.00538 0.00403 0.00000736
    ## 2 SELECTED 111461 0.00497 0.00376 0.0000113

``` r
pi.wvd.AE.dataframe %>% t_test(formula = PI ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df   p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>     <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      36.1 148453. 2.01e-283 two.sided   0.000463 0.000438 0.000488

``` r
pi.wvd.AE.dataframe %>% t_test(formula = PI ~ GROUP, alternative = "greater", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df   p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>     <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      36.1 148453. 1.00e-283 greater     0.000463 0.000442      Inf

``` r
group_by(pi.wvd.AE.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PI, na.rm = TRUE),
    sd = sd(PI, na.rm = TRUE),
    se=std.error(PI, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP     count    mean      sd         se
    ##   <fct>     <int>   <dbl>   <dbl>      <dbl>
    ## 1 WILD     185960 0.00430 0.00316 0.00000733
    ## 2 SELECTED  73632 0.00383 0.00286 0.0000105

``` r
filename <- paste(prefix,"mean.pop.wild.windowed.pi.txt", sep="")
wpi <- read.table(filename, header = TRUE)
filename <- paste(prefix,"mean.pop.selected.windowed.pi.txt", sep="")
dpi <- read.table(filename, header = TRUE)

wpi <- wpi %>% dplyr::rename(WILD_PI= PI)
dpi <- dpi %>% dplyr::rename(SELECTED_PI= PI)
w.d.pi <- join(wpi,dpi, by =c("CHROM", "BIN_START"),type="full")


w.d.pi$DIFF <- w.d.pi$WILD_PI - w.d.pi$SELECTED_PI
w.d.pi$BIN_START <- w.d.pi$BIN_START -1
```

``` r
SNP.pi<-c(1:(nrow(w.d.pi)))
mydf.pi<-data.frame(SNP.pi,w.d.pi)
mydf.pi$CHR <- mydf.pi$CHROM
mydf.pi %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf.pi
mydf.pi$CHR <- as.numeric(mydf.pi$CHR)  

m2.pi <-ggman(mydf.pi,chrom="CHR",bp="BIN_START",pvalue="WILD_PI",snp="SNP.pi",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm.pi <- m2.pi[[1]]
dfm.pi$chrom_alt <- factor(dfm.pi$chrom_alt, levels=c(0,1))

dfmsplit.pi <- split(dfm.pi, dfm.pi$chrom)
xbreaks <- sapply(dfmsplit.pi,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m2.pi.total <- ggplot(dfm.pi, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm.pi$index)+10)) +
  guides(colour = FALSE,alpha=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm.pi$marker+0.001,na.rm=TRUE)))+
  labs(x = "Chromosome") +
  ggtitle("Wild") +theme_classic()+
  labs(y=expression(pi)) +
  scale_color_manual(values=c("#D55E00", "#8f3e00")) 
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
m2.pi <-m2.pi.total



m3.pi <-ggman(mydf.pi,chrom="CHR",bp="BIN_START",pvalue="SELECTED_PI",snp="SNP.pi",logTransform=FALSE,sigLine=NA, ymax=1, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm3.pi <- m3.pi[[1]]
dfm3.pi$chrom_alt <- factor(dfm3.pi$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm3.pi, dfm3.pi$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m3.pi.total <- ggplot(dfm3.pi, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm3.pi$index)+10)) +
  guides(colour = FALSE,alpha=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm3.pi$marker+0.001,na.rm=TRUE)))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Selected") +theme_classic()+
  labs(y=expression(pi)) + 
  scale_color_manual(values=c("#009E73", "#00543d"))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
m3.pi <-m3.pi.total


md4.pi <-ggman(mydf.pi,chrom="CHR",bp="BIN_START",pvalue="DIFF",snp="SNP.pi",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf.pi, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", :
    ## NaNs produced

    ## Warning in ggman(mydf.pi, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", :
    ## `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")`
    ## instead.

``` r
dfm4.pi <- md4.pi[[1]]
dfm4.pi$chrom_alt <- factor(dfm4.pi$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm4.pi, dfm4.pi$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

md4.pi <- ggplot(dfm4.pi, aes(x= index, y=marker, colour = as.factor(chrom_alt),size=(0.1+ abs(marker))^2))+
  geom_point(alpha=0.85)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm4.pi$index)+10)) +
  guides(colour = FALSE,alpha=FALSE,size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(-0.006,0.007))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Wild - Selected") +theme_classic()+
  labs(y=expression("Difference in " ~ pi)) + 
  scale_color_manual(values=c("#7d5700", "#634500"))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
png(filename=paste(prefix.OF,"Figure.2.pi_comparison.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

layout <- "
AAADD
BBBDD
CCCEE
CCCFF
"
m2.pi+m3.pi+md4.pi+bd.pi+bc.pi+bc.pi.ae +plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')
```

    ## Warning: Removed 32 rows containing missing values (geom_point).

    ## Warning: Removed 90 rows containing missing values (geom_point).

    ## Warning: Removed 122 rows containing missing values (geom_point).

``` r
dev.off()
```

    ## png 
    ##   2

## Tajima’s D

``` bash
conda create -n vk python=3.8 numpy scipy
conda activate vk
pip install VCF-kit
```

``` bash
ls ./other_files/*.temp | sed 's/.temp//g' | sed 's:\./other_files/::g' | parallel "bcftools view --threads 4 -S ./other_files/{}.temp ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz| vk tajima 10000 10000 - > ./Population_Genomics/Intermediate_Files/{}.tajima"
```

``` bash
pwd
cd ./Population_Genomics/Intermediate_Files/
pwd
```

    ## /home/jpuritz/2022_Population_Genomics
    ## /home/jpuritz/2022_Population_Genomics/Population_Genomics/Intermediate_Files

``` bash
cd ./Population_Genomics/Intermediate_Files/
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.tajima) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.tajima) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.tajima) <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.tajima) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.tajima) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.tajima) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.tajima) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.tajima) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.tajima) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS2.tajima) <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.tajima) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.tajima) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.tajima) > ../pop.windowed.tajima

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.tajima) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.tajima) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.tajima) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.tajima) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC_VA.tajima) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.tajima)  <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.tajima) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.tajima)  > ../pop.wild.windowed.tajima

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.tajima) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.tajima) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.tajima) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS2.tajima)  <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.tajima) > ../pop.selected.windowed.tajima

cd ..

bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.tajima | sort -k1,1 -k2,2n ) -c 4,5,6,8 -d -1 -o mean,mean,mean,first > mean.pop.wild.windowed.tajima

bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.tajima | sort -k1,1 -k2,2n ) -c 4,5,6,8 -d -1 -o mean,mean,mean,first > mean.pop.selected.windowed.tajima

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tGROUP") mean.pop.wild.windowed.tajima mean.pop.selected.windowed.tajima > mean.pop.windowed.tajima

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tGROUP") mean.pop.wild.windowed.tajima > ./Intermediate_Files/mean.pop.wild.windowed.tajima.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tTAJD\tGROUP") mean.pop.selected.windowed.tajima > ./Intermediate_Files/mean.pop.selected.windowed.tajima.txt


mkdir -p ./Output/tajimasD

cp pop.windowed.tajima mean.pop.windowed.tajima mean.pop.wild.windowed.tajima.txt mean.pop.selected.windowed.tajima.txt ./Output/tajimasD
```

``` r
tajd.wvd.dataframe<-read.table(paste(prefix,"pop.windowed.tajima", sep=""), sep="\t", header=T)
summary(tajd.wvd.dataframe)
```

    ##          CHROM          BIN_START            BIN_END             N_SITES     
    ##  NC_035784.1: 93823   Min.   :        0   Min.   :    10000   Min.   :  1.0  
    ##  NC_035782.1: 63387   1st Qu.: 18480000   1st Qu.: 18490000   1st Qu.: 61.0  
    ##  NC_035780.1: 57737   Median : 36820000   Median : 36830000   Median :156.0  
    ##  NC_035783.1: 56871   Mean   : 38160118   Mean   : 38170118   Mean   :152.5  
    ##  NC_035781.1: 55770   3rd Qu.: 55060000   3rd Qu.: 55070000   3rd Qu.:229.0  
    ##  NC_035788.1: 44910   Max.   :104140000   Max.   :104150000   Max.   :631.0  
    ##  (Other)    :102655                                                          
    ##      N_SNPs           TAJD               POP              GROUP       
    ##  Min.   :  1.0   Min.   :-2.58444   HC     : 36957   SELECTED:182084  
    ##  1st Qu.: 37.0   1st Qu.:-0.06336   CS     : 36692   WILD    :293069  
    ##  Median :101.0   Median : 0.34627   CLP    : 36628                    
    ##  Mean   :104.1   Mean   : 0.32090   HC-VA  : 36609                    
    ##  3rd Qu.:157.0   3rd Qu.: 0.74076   HI     : 36590                    
    ##  Max.   :543.0   Max.   : 3.54663   DEBY   : 36564                    
    ##                                     (Other):255113

``` r
tajd.wvd.dataframe$GROUP <- factor(tajd.wvd.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
tajd.wvd.dataframe$POP <- factor(tajd.wvd.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

tajd.wvd.mean.dataframe<-read.table(paste(prefix,"mean.pop.windowed.tajima", sep=""), sep="\t", header=T)
tajd.wvd.mean.dataframe$GROUP <- factor(tajd.wvd.mean.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
summary(tajd.wvd.mean.dataframe)
```

    ##          CHROM         BIN_START            BIN_END             N_SITES   
    ##  NC_035784.1:14562   Min.   :        0   Min.   :    10000   Min.   :  1  
    ##  NC_035782.1: 9861   1st Qu.: 18460000   1st Qu.: 18470000   1st Qu.: 55  
    ##  NC_035780.1: 8976   Median : 36760000   Median : 36770000   Median :153  
    ##  NC_035783.1: 8801   Mean   : 38140951   Mean   : 38150951   Mean   :150  
    ##  NC_035781.1: 8645   3rd Qu.: 55050000   3rd Qu.: 55060000   3rd Qu.:227  
    ##  NC_035788.1: 7141   Max.   :104140000   Max.   :104150000   Max.   :631  
    ##  (Other)    :16338                                                        
    ##      N_SNPs            TAJD              GROUP      
    ##  Min.   :  1.00   Min.   :-1.6809   WILD    :37172  
    ##  1st Qu.: 35.25   1st Qu.: 0.1152   SELECTED:37152  
    ##  Median :100.25   Median : 0.3792                   
    ##  Mean   :101.01   Mean   : 0.3379                   
    ##  3rd Qu.:153.62   3rd Qu.: 0.6134                   
    ##  Max.   :470.12   Max.   : 2.2815                   
    ## 

``` r
bd <-ggplot(tajd.wvd.dataframe, aes(x=TAJD,y = POP))+
  #geom_violin(aes(color=GROUP,fill= GROUP)) +
  #stat_summary(fun=mean, geom="point", shape=23, size=4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.005)+
  geom_boxplot(aes(fill=GROUP), width=0.75,outlier.shape = 23, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(tajd.wvd.dataframe$POP)))+
  #ylab("Population/Line")+
  ggtitle("Wild vs. Selected Tajima's *D*")+
  theme_classic() +
  labs(x="Tajima's *D*") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank(), plot.title = element_markdown()) 

bd <- bd + guides(color = "none", fill ="none")
```

``` r
bc <-ggplot(tajd.wvd.mean.dataframe, aes(x=TAJD,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP)) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = 23, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(tajd.wvd.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. Selected Tajima's *D*")+
  theme_classic() +
  labs(x="Tajima's *D*") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())
```

``` r
tajd.wvd.dataframe %>% t_test(formula = TAJD ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -120. 368890.       0 two.sided     -0.245   -0.249   -0.241

``` r
tajd.wvd.dataframe %>% t_test(formula = TAJD ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -120. 368890.       0 less          -0.245     -Inf   -0.242

``` r
group_by(tajd.wvd.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(TAJD, na.rm = TRUE),
    sd = sd(TAJD, na.rm = TRUE),
    se=std.error(TAJD, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP     count  mean    sd      se
    ##   <fct>     <int> <dbl> <dbl>   <dbl>
    ## 1 WILD     293069 0.227 0.658 0.00122
    ## 2 SELECTED 182084 0.472 0.698 0.00163

``` r
group_by(tajd.wvd.dataframe, POP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(TAJD, na.rm = TRUE),
    sd = sd(TAJD, na.rm = TRUE),
    se=std.error(TAJD, na.rm = TRUE) )
```

    ## # A tibble: 13 × 5
    ##    POP   count  mean    sd      se
    ##    <fct> <int> <dbl> <dbl>   <dbl>
    ##  1 HI    36590 0.213 0.661 0.00346
    ##  2 SM    36543 0.230 0.703 0.00368
    ##  3 HC    36957 0.358 0.679 0.00353
    ##  4 CS    36692 0.104 0.652 0.00340
    ##  5 CLP   36628 0.142 0.654 0.00341
    ##  6 HC-VA 36609 0.143 0.652 0.00341
    ##  7 CL    36507 0.311 0.603 0.00315
    ##  8 SL    36543 0.314 0.605 0.00317
    ##  9 UMFS  36413 0.543 0.739 0.00388
    ## 10 NEH   36257 0.629 0.801 0.00421
    ## 11 DEBY  36564 0.323 0.617 0.00322
    ## 12 LOLA  36541 0.319 0.617 0.00323
    ## 13 OBOYS 36309 0.550 0.635 0.00333

``` r
filename <- paste(prefix.I,"mean.pop.wild.windowed.tajima.txt", sep="")
wtajd <- read.table(filename, header = TRUE)
filename <- paste(prefix.I,"mean.pop.selected.windowed.tajima.txt", sep="")
dtajd <- read.table(filename, header = TRUE)

wtajd <- wtajd %>% dplyr::rename(WILD_TAJD= TAJD)
dtajd <- dtajd %>% dplyr::rename(SELECTED_TAJD= TAJD)
w.d.tajd <- join(wtajd,dtajd, by =c("CHROM", "BIN_START"),type="full")


w.d.tajd$DIFF <- w.d.tajd$WILD_TAJD - w.d.tajd$SELECTED_TAJD
```

``` r
SNP<-c(1:(nrow(w.d.tajd)))
mydf<-data.frame(SNP,w.d.tajd)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf
mydf$CHR <- as.numeric(mydf$CHR)  

m2 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WILD_TAJD",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "WILD_TAJD", :
    ## NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- m2[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m2.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),max(dfm$marker+0.001,na.rm=TRUE)))+
  labs(x = "Chromosome") +
  ggtitle("Wild") +theme_classic()+
  labs(y="Tajima's *D*") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#D55E00", "#8f3e00")) 
m2 <-m2.total



m3 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="SELECTED_TAJD",snp="SNP",logTransform=FALSE,sigLine=NA, ymax=1, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- m3[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m3.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),max(dfm$marker+0.001,na.rm=TRUE)))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Selected") +theme_classic()+
  labs(y="Tajima's *D*") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#009E73", "#00543d"))

m3 <-m3.total


md4 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="DIFF",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", : NaNs
    ## produced

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", :
    ## `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")`
    ## instead.

``` r
dfm <- md4[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

md4 <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt),size=(0.1+ abs(marker))^2))+
  geom_point(alpha=0.85)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.005,na.rm=TRUE),max(dfm$marker+0.005,na.rm=TRUE)))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Wild - Selected") +theme_classic()+
  labs(y="Difference in Tajima's *D*") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#7d5700", "#634500"))


png(filename=paste(prefix.OFS,"Figure.S5.TajimasD.Comparison.W.S.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AAAD
BBBD
CCCE
CCCE
"
#######(m2/m3/md4 | bd/bc) +plot_layout(design=layout)
m2+m3+md4+bd+bc+plot_layout(design=layout) + plot_annotation(tag_levels = 'A')
```

    ## Warning: Removed 20 rows containing missing values (geom_point).

    ## Warning: Removed 20 rows containing missing values (geom_point).

``` r
dev.off()
```

    ## png 
    ##   2

## Heterozygosity

``` bash
ln -s ~/2021_Genome_Submissions/Masked/raw.vcf/filtered/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz ./Population_Genomics/

INDS=$(seq -s, 0 5 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/CL.MASKED.popstats &
INDS=$(seq -s, 6 11 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/CLP.MASKED.popstats &
INDS=$(seq -s, 12 17 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/CS.MASKED.popstats &
INDS=$(seq -s, 18 23 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/DEBY.MASKED.popstats &
INDS=$(seq -s, 24 29 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/HC.MASKED.popstats &
INDS=$(seq -s, 30 35 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/HC-VA.MASKED.popstats &
INDS=$(seq -s, 39 44 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/HI.MASKED.popstats &
INDS=$(seq -s, 45 49 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/LM.MASKED.popstats &
INDS=$(seq -s, 50 55 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/LOLA.MASKED.popstats
INDS=$(seq -s, 56 61 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/NEH.MASKED.popstats &
INDS=$(seq -s, 66 71 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/OBOYS.MASKED.popstats
INDS=$(seq -s, 72 77 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/SL.MASKED.popstats &
INDS=$(seq -s, 78 83 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/SM.MASKED.popstats
INDS=$(seq -s, 84 89 )
popStats -y GT --file <( bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf.gz) -t $INDS > ./Population_Genomics/Intermediate_Files/UMFS.MASKED.popstats &
```

``` bash
bedtools makewindows -g ./other_files/genome.file -w 10000 | mawk '!/NC_007175.2/' > ./other_files/10kb.bed
AWK1='{print $1 "\t" $2-1 "\t" $2 "\t" $3 "\t" $4}'
ls ./Population_Genomics/Intermediate_Files/*.MASKED.popstats | sed 's/.MASKED.popstats//g' | sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "cut -f1,2,4,5 ./Population_Genomics/Intermediate_Files/{}.MASKED.popstats | mawk '$AWK1' > ./Population_Genomics/Intermediate_Files/{}.het.bed"

ls ./Population_Genomics/Intermediate_Files/*.het.bed | sed 's/.het.bed//g'  |sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "bedtools intersect -a ./Population_Genomics/Intermediate_Files/{}.het.bed -b ./Population_Genomics/SNP.MASKED.AE.TRSdp5g75.nDNA.g1.maf05.max2alleles.bed > ./Population_Genomics/Intermediate_Files/{}.het.ae.bed"

ls ./Population_Genomics/Intermediate_Files/*.het.bed | sed 's/.het.bed//g'  |sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "bedtools intersect -b ./Population_Genomics/Intermediate_Files/{}.het.bed -a ./other_files/10kb.bed -wa -wb | bedtools merge -d -1 -c 7,8 -o mean > ./Population_Genomics/Intermediate_Files/{}.het.per10kb.bed"

ls ./Population_Genomics/Intermediate_Files/*.het.ae.bed | sed 's/.het.ae.bed//g'  |sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "bedtools intersect -b ./Population_Genomics/Intermediate_Files/{}.het.ae.bed -a ./other_files/10kb.bed -wa -wb | bedtools merge -d -1 -c 7,8 -o mean > ./Population_Genomics/Intermediate_Files/{}.het.ae.per10kb.bed"

cd ./Population_Genomics/Intermediate_Files/


cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tHET\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC-VA.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.het.per10kb.bed)  <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.het.per10kb.bed)  > pop.wild.windowed.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tHET\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC-VA.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.het.ae.per10kb.bed)   <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.het.ae.per10kb.bed)  > pop.wild.windowed.ae.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tHET\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS.het.per10kb.bed)  <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.het.per10kb.bed) > pop.selected.windowed.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tN_SITES\tN_SNPs\tHET\tPOP\tGROUP")  <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.het.ae.per10kb.bed) > pop.selected.windowed.ae.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCL\tWILD"}' CL.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tDEBY\tSELECTED"}' DEBY.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC-VA.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tLOLA\tSELECTED"}' LOLA.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tOBOYS\tSELECTED"}' OBOYS.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tSL\tWILD"}' SL.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.het.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.het.per10kb.bed) > ../pop.windowed.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tPOP\tGROUP") <(mawk '{if (NR>1)print $0 "\tCLP\tWILD"}' CLP.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tCS\tWILD"}' CS.het.ae.per10kb.bed)  <(mawk '{if (NR>1)print $0 "\tHC\tWILD"}' HC.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHC-VA\tWILD"}' HC-VA.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tHI\tWILD"}' HI.het.ae.per10kb.bed)  <(mawk '{if (NR>1)print $0 "\tNEH\tSELECTED"}' NEH.het.ae.per10kb.bed)  <(mawk '{if (NR>1)print $0 "\tSM\tWILD"}' SM.het.ae.per10kb.bed) <(mawk '{if (NR>1)print $0 "\tUMFS\tSELECTED"}' UMFS.het.ae.per10kb.bed) > ../pop.windowed.ae.het


bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.het| sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.wild.windowed.het
bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.het | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.selected.windowed.het

bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.ae.het| sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.wild.windowed.ae.het
bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.ae.het | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.selected.windowed.ae.het

cd .. 

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.wild.windowed.het mean.pop.selected.windowed.het > mean.pop.windowed.het

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.wild.windowed.ae.het mean.pop.selected.windowed.ae.het > mean.pop.windowed.ae.het


cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.wild.windowed.het > mean.pop.wild.windowed.het.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.selected.windowed.het > mean.pop.selected.windowed.het.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.wild.windowed.ae.het > mean.pop.wild.windowed.ae.het.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tEXP_HET\tOBS_HET\tGROUP") mean.pop.selected.windowed.ae.het > mean.pop.selected.windowed.ae.het.txt

mkdir -p ./Output/Heterozygosity
cp pop.windowed.het mean.pop.windowed.het pop.windowed.ae.het mean.pop.windowed.ae.het ./Output/Heterozygosity
```

``` r
filename <- paste(prefix,"pop.windowed.het", sep="")
het.wvd.dataframe<-read.table(filename, sep="\t", header=T)
summary(het.wvd.dataframe)
```

    ##          CHROM          BIN_START            BIN_END             EXP_HET      
    ##  NC_035784.1: 94666   Min.   :        0   Min.   :    10000   Min.   :0.0000  
    ##  NC_035782.1: 64116   1st Qu.: 18460000   1st Qu.: 18470000   1st Qu.:0.1846  
    ##  NC_035780.1: 58357   Median : 36760000   Median : 36770000   Median :0.2213  
    ##  NC_035783.1: 57213   Mean   : 38141147   Mean   : 38151147   Mean   :0.2165  
    ##  NC_035781.1: 56199   3rd Qu.: 55050000   3rd Qu.: 55060000   3rd Qu.:0.2527  
    ##  NC_035788.1: 46436   Max.   :104140000   Max.   :104150000   Max.   :0.5000  
    ##  (Other)    :106236                                                           
    ##     OBS_HET            POP              GROUP       
    ##  Min.   :0.0000   CL     : 37171   SELECTED:185855  
    ##  1st Qu.:0.1962   CLP    : 37171   WILD    :297368  
    ##  Median :0.2436   CS     : 37171                    
    ##  Mean   :0.2460   DEBY   : 37171                    
    ##  3rd Qu.:0.2896   HC     : 37171                    
    ##  Max.   :1.0000   HC-VA  : 37171                    
    ##                   (Other):260197

``` r
het.wvd.dataframe$GROUP <- factor(het.wvd.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
het.wvd.dataframe$POP <- factor(het.wvd.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

filename <- paste(prefix,"mean.pop.windowed.het", sep="")
het.wvd.mean.dataframe<-read.table(filename, sep="\t", header=T)
het.wvd.mean.dataframe$GROUP <- factor(het.wvd.mean.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
summary(het.wvd.mean.dataframe)
```

    ##          CHROM         BIN_START            BIN_END             EXP_HET      
    ##  NC_035784.1:14564   Min.   :        0   Min.   :    10000   Min.   :0.0000  
    ##  NC_035782.1: 9864   1st Qu.: 18460000   1st Qu.: 18470000   1st Qu.:0.1949  
    ##  NC_035780.1: 8978   Median : 36760000   Median : 36770000   Median :0.2167  
    ##  NC_035783.1: 8802   Mean   : 38141147   Mean   : 38151147   Mean   :0.2155  
    ##  NC_035781.1: 8646   3rd Qu.: 55047500   3rd Qu.: 55057500   3rd Qu.:0.2369  
    ##  NC_035788.1: 7144   Max.   :104140000   Max.   :104150000   Max.   :0.4917  
    ##  (Other)    :16344                                                           
    ##     OBS_HET            GROUP      
    ##  Min.   :0.0000   WILD    :37171  
    ##  1st Qu.:0.2182   SELECTED:37171  
    ##  Median :0.2444                   
    ##  Mean   :0.2470                   
    ##  3rd Qu.:0.2708                   
    ##  Max.   :0.9000                   
    ## 

``` r
filename <- paste(prefix,"pop.windowed.ae.het", sep="")
het.wvd.ae.dataframe<-read.table(filename, sep="\t", header=T)
summary(het.wvd.ae.dataframe)
```

    ##          CHROM         BIN_START            BIN_END             EXP_HET      
    ##  NC_035784.1:58088   Min.   :        0   Min.   :    10000   Min.   :0.0000  
    ##  NC_035782.1:39328   1st Qu.: 18470000   1st Qu.: 18480000   1st Qu.:0.2228  
    ##  NC_035780.1:35736   Median : 36780000   Median : 36790000   Median :0.2614  
    ##  NC_035783.1:35112   Mean   : 38148766   Mean   : 38158766   Mean   :0.2546  
    ##  NC_035781.1:34456   3rd Qu.: 55050000   3rd Qu.: 55060000   3rd Qu.:0.2931  
    ##  NC_035788.1:28232   Max.   :104140000   Max.   :104150000   Max.   :0.5000  
    ##  (Other)    :64584                                                           
    ##     OBS_HET            POP             GROUP       
    ##  Min.   :0.0000   CLP    :36942   SELECTED: 73884  
    ##  1st Qu.:0.2350   CS     :36942   WILD    :221652  
    ##  Median :0.2863   HC     :36942                    
    ##  Mean   :0.2885   HC-VA  :36942                    
    ##  3rd Qu.:0.3350   HI     :36942                    
    ##  Max.   :1.0000   NEH    :36942                    
    ##                   (Other):73884

``` r
het.wvd.ae.dataframe$GROUP <- factor(het.wvd.ae.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
het.wvd.ae.dataframe$POP <- factor(het.wvd.ae.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"))

filename <- paste(prefix,"mean.pop.windowed.ae.het", sep="")
het.wvd.mean.ae.dataframe<-read.table(filename, sep="\t", header=T)
het.wvd.mean.ae.dataframe$GROUP <- factor(het.wvd.mean.ae.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
summary(het.wvd.mean.ae.dataframe)
```

    ##          CHROM         BIN_START            BIN_END             EXP_HET      
    ##  NC_035784.1:14522   Min.   :        0   Min.   :    10000   Min.   :0.0000  
    ##  NC_035782.1: 9832   1st Qu.: 18470000   1st Qu.: 18480000   1st Qu.:0.2222  
    ##  NC_035780.1: 8934   Median : 36780000   Median : 36790000   Median :0.2529  
    ##  NC_035783.1: 8778   Mean   : 38148766   Mean   : 38158766   Mean   :0.2490  
    ##  NC_035781.1: 8614   3rd Qu.: 55050000   3rd Qu.: 55060000   3rd Qu.:0.2794  
    ##  NC_035788.1: 7058   Max.   :104140000   Max.   :104150000   Max.   :0.5000  
    ##  (Other)    :16146                                                           
    ##     OBS_HET            GROUP      
    ##  Min.   :0.0000   WILD    :36942  
    ##  1st Qu.:0.2500   SELECTED:36942  
    ##  Median :0.2875                   
    ##  Mean   :0.2898                   
    ##  3rd Qu.:0.3252                   
    ##  Max.   :1.0000                   
    ## 

``` r
bd <-ggplot(het.wvd.dataframe, aes(x=OBS_HET,y = POP))+
  #geom_violin(aes(color=GROUP,fill= GROUP)) +
  #stat_summary(fun=mean, geom="point", shape=23, size=4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.005)+
  geom_boxplot(aes(fill=GROUP), width=0.75,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(het.wvd.dataframe$POP)))+
  #ylab("Population/Line")+
  ggtitle("Wild vs. Selected Heterozygosity")+
  coord_cartesian(xlim=c(0,0.5))+
  theme_classic() +
  labs(x="Observed Heterozygosity") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank(), plot.title = element_markdown()) 


bd <- bd + guides(color = "none", fill ="none")

bd.ae <-ggplot(het.wvd.ae.dataframe, aes(x=OBS_HET,y = POP))+
  #geom_violin(aes(color=GROUP,fill= GROUP)) +
  #stat_summary(fun=mean, geom="point", shape=23, size=4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.005)+
  geom_boxplot(aes(fill=GROUP), width=0.75,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(het.wvd.dataframe$POP)))+
  #ylab("Population/Line")+
  ggtitle("Wild vs. Selected Heterozygosity")+
  coord_cartesian(xlim=c(0,0.5))+
  theme_classic() +
  labs(x="Observed Heterozygosity") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank(), plot.title = element_markdown()) 
```

``` r
bc <-ggplot(het.wvd.mean.dataframe, aes(x=OBS_HET,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(het.wvd.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. Selected Heterozygosity")+
  theme_classic() +
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  coord_cartesian(xlim=c(0,0.5))+
  labs(x="Observed Heterozygosity") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())

bc.ae <-ggplot(het.wvd.mean.ae.dataframe, aes(x=OBS_HET,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(het.wvd.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. Selected Heterozygosity using Atlantic Samples Only")+
  theme_classic() +
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  coord_cartesian(xlim=c(0,0.5))+
  labs(x="Observed Heterozygosity") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())
```

``` r
het.wvd.mean.dataframe %>% t_test(formula = OBS_HET ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -19.5 74319. 3.56e-84 two.sided   -0.00869 -0.00956 -0.00781

``` r
het.wvd.mean.dataframe %>% t_test(formula = OBS_HET ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -19.5 74319. 1.78e-84 less        -0.00869     -Inf -0.00795

``` r
group_by(het.wvd.mean.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(OBS_HET, na.rm = TRUE),
    sd = sd(OBS_HET, na.rm = TRUE),
    se=std.error(OBS_HET, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean     sd       se
    ##   <fct>    <int> <dbl>  <dbl>    <dbl>
    ## 1 WILD     37171 0.243 0.0603 0.000313
    ## 2 SELECTED 37171 0.251 0.0613 0.000318

``` r
group_by(het.wvd.dataframe, POP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(OBS_HET, na.rm = TRUE),
    sd = sd(OBS_HET, na.rm = TRUE),
    se=std.error(OBS_HET, na.rm = TRUE) )
```

    ## # A tibble: 13 × 5
    ##    POP   count  mean     sd       se
    ##    <fct> <int> <dbl>  <dbl>    <dbl>
    ##  1 HI    37171 0.243 0.0943 0.000489
    ##  2 SM    37171 0.236 0.0970 0.000503
    ##  3 HC    37171 0.246 0.0923 0.000479
    ##  4 CS    37171 0.246 0.0902 0.000468
    ##  5 CLP   37171 0.247 0.0925 0.000480
    ##  6 HC-VA 37171 0.246 0.0920 0.000477
    ##  7 CL    37171 0.238 0.0909 0.000471
    ##  8 SL    37171 0.239 0.0912 0.000473
    ##  9 UMFS  37171 0.252 0.106  0.000552
    ## 10 NEH   37171 0.238 0.111  0.000577
    ## 11 DEBY  37171 0.265 0.0937 0.000486
    ## 12 LOLA  37171 0.255 0.0943 0.000489
    ## 13 OBOYS 37171 0.246 0.100  0.000521

``` r
test<- het.wvd.dataframe[which(het.wvd.dataframe$POP != "LOLA" & het.wvd.dataframe$POP != "DEBY"),]

test %>% t_test(formula = OBS_HET ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -7.77 178740. 7.97e-15 two.sided   -0.00280 -0.00351 -0.00210

``` r
test %>% t_test(formula = OBS_HET ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic    t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>   <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -7.77 178740. 3.99e-15 less        -0.00280     -Inf -0.00221

``` r
group_by(test, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(OBS_HET, na.rm = TRUE),
    sd = sd(OBS_HET, na.rm = TRUE),
    se=std.error(OBS_HET, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP     count  mean     sd       se
    ##   <fct>     <int> <dbl>  <dbl>    <dbl>
    ## 1 WILD     297368 0.243 0.0927 0.000170
    ## 2 SELECTED 111513 0.245 0.106  0.000318

``` r
het.wvd.mean.ae.dataframe %>% t_test(formula = OBS_HET ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -9.12 67495. 7.40e-20 two.sided   -0.00530 -0.00644 -0.00416

``` r
het.wvd.mean.ae.dataframe %>% t_test(formula = OBS_HET ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df  p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>    <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -9.12 67495. 3.70e-20 less        -0.00530     -Inf -0.00434

``` r
group_by(het.wvd.mean.ae.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(OBS_HET, na.rm = TRUE),
    sd = sd(OBS_HET, na.rm = TRUE),
    se=std.error(OBS_HET, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean     sd       se
    ##   <fct>    <int> <dbl>  <dbl>    <dbl>
    ## 1 WILD     36942 0.287 0.0657 0.000342
    ## 2 SELECTED 36942 0.292 0.0903 0.000470

``` r
filename <- paste(prefix,"mean.pop.wild.windowed.het.txt", sep="")
whet <- read.table(filename, header = TRUE)
filename <- paste(prefix,"mean.pop.selected.windowed.het.txt", sep="")
dhet <- read.table(filename, header = TRUE)

whet <- whet %>% dplyr::rename(WILD_OBS_HET= OBS_HET)
dhet <- dhet %>% dplyr::rename(SELECTED_OBS_HET= OBS_HET)
w.d.het <- join(whet,dhet, by =c("CHROM", "BIN_START"),type="full")


w.d.het$DIFF <- w.d.het$WILD_OBS_HET - w.d.het$SELECTED_OBS_HET
```

``` r
SNP<-c(1:(nrow(w.d.het)))
mydf<-data.frame(SNP,w.d.het)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf
mydf$CHR <- as.numeric(mydf$CHR)  

m2 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WILD_OBS_HET",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- m2[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m2.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),0.9))+
  labs(x = "Chromosome") +
  ggtitle("Wild") +theme_classic()+
  labs(y="Observed Heterozygosity") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#D55E00", "#8f3e00")) 
m2 <-m2.total



m3 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="SELECTED_OBS_HET",snp="SNP",logTransform=FALSE,sigLine=NA, ymax=1, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- m3[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

m3.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),0.9))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Selected") +theme_classic()+
  labs(y="Observed Heterozygosity") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#009E73", "#00543d"))

m3 <-m3.total


md4 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="DIFF",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", : NaNs
    ## produced

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "DIFF", :
    ## `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")`
    ## instead.

``` r
dfm <- md4[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

md4 <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt),size=(0.1+ abs(marker))^2))+
  geom_point(alpha=0.85)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.4,0.4))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") +
  ggtitle("Wild - Domestic") +theme_classic()+
  labs(y="Difference in Observed Heterozygosity") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#7d5700", "#634500"))



png(filename=paste(prefix.OFS,"Figure.S4.Heterozygosity.Comparison.W.S.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AAADD
BBBDD
CCCEE
CCCFF
"
m2+m3+md4+bd+bc+bc.ae +plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')
dev.off()
```

    ## png 
    ##   2

## CNV Diversity

``` r
mask_cn1 <- read.table("./Population_Genomics/masked.cnv.tab",header = FALSE)

group <- popmap.ni$POP

mask_cn1 %>% 
  mutate(V1 = str_replace(V1, "NC_035780.1", "1")) %>% 
  mutate(V1 = str_replace(V1, "NC_035781.1", "2")) %>% 
  mutate(V1 = str_replace(V1, "NC_035782.1", "3")) %>% 
  mutate(V1 = str_replace(V1, "NC_035783.1", "4")) %>% 
  mutate(V1 = str_replace(V1, "NC_035784.1", "5")) %>% 
  mutate(V1 = str_replace(V1, "NC_035785.1", "6")) %>% 
  mutate(V1 = str_replace(V1, "NC_035786.1", "7")) %>% 
  mutate(V1 = str_replace(V1, "NC_035787.1", "8")) %>%
  mutate(V1 = str_replace(V1, "NC_035788.1", "9")) %>% 
  mutate(V1 = str_replace(V1, "NC_035789.1", "10"))  -> mask_cn1

mask_cnvs <- paste(mask_cn1$V1,mask_cn1$V2,mask_cn1$V3,sep="_")

mask_cnv_df <- t(mask_cn1[,4:93])

rownames(mask_cnv_df) <- popmap$IND
colnames(mask_cnv_df) <- mask_cnvs

mask_gen <-df2genind(mask_cnv_df,ploidy = 1)

pop(mask_gen) <- popmap$POP

MAF=0.05

mask_gen.ni <- popsub(mask_gen,c("CL","CLP","CS","DEBY","HC","HC-VA","HI","LOLA","NEH","OBOYS","SL","SM", "UMFS"))
m.freq.ni <- minorAllele(mask_gen.ni)
m.freq.AE <- minorAllele(popsub(mask_gen.ni,c("CLP","CS","HC","HC-VA","HI","SM","UMFS","NEH")))
m.freq.wild <- minorAllele(popsub(mask_gen.ni,c("CL","CLP","CS","HC","HC-VA","HI","SL","SM")))
m.freq.wildAE <- minorAllele(popsub(mask_gen.ni,c("CLP","CS","HC","HC-VA","HI","SM")))
m.freq.wildGM <- minorAllele(popsub(mask_gen.ni,c("CL","SL")))
m.freq.selected <- minorAllele(popsub(mask_gen.ni,c("DEBY","LOLA","NEH", "OBOYS","UMFS")))
m.freq.selectedAE <- minorAllele(popsub(mask_gen.ni,c("NEH","UMFS")))
m.freq.selectedGM <- minorAllele(popsub(mask_gen.ni,c("OBOYS")))
                                 
m.freq.df <- cbind.data.frame(m.freq.ni,m.freq.wild,m.freq.AE,m.freq.wildAE,m.freq.wildGM,m.freq.selected,m.freq.selectedAE,m.freq.selectedGM)
colnames(m.freq.df) <- c("Overall","Wild","Atlantic","WildAE", "WildGM", "Selected","SelectedAE","SelectedGM")

#write.table(m.freq.df, file='./Population_Genomics/test.tsv', quote=FALSE, sep='\t')

mask_gen.ae <- mask_gen.ni[,loc=which(m.freq.AE > MAF & m.freq.AE < 1-MAF)]
mask_gen.niw <- mask_gen.ni[,loc=which(m.freq.ni > MAF & m.freq.ni < 1-MAF)]
mask_gen.wildAE <- popsub(mask_gen.ni,c("CLP","CS","HC","HC-VA","HI","SM"))
mask_gen.wildAE <- mask_gen.wildAE[,loc=which(m.freq.wildAE > MAF & m.freq.wildAE < 1-MAF)]


cnv.ar <- allelic.richness(mask_gen.niw, diploid=FALSE)
cnv.ae.ar <- allelic.richness(mask_gen.ae, diploid=FALSE)

cnv.ar.wild <- cnv.ar$Ar[c("CL","CLP","CS","HC","HC-VA","HI","SL","SM")]
cnv.ar.ae.wild <- cnv.ae.ar$Ar[c("CLP","CS","HC","HC-VA","HI","SM")]

cnv.ar.selected <- cnv.ar$Ar[c("DEBY","LOLA","NEH", "OBOYS","UMFS")]
cnv.ar.ae.selected <- cnv.ae.ar$Ar[c("NEH", "UMFS")]

cnv.ar.wild.df <- gather(cnv.ar.wild,"POP", "AR")
cnv.ar.wild.df$Group <- rep("WILD",length(cnv.ar.wild.df$POP) )

cnv.ar.ae.wild.df <- gather(cnv.ar.ae.wild,"POP", "AR")
cnv.ar.ae.wild.df$Group <- rep("WILD",length(cnv.ar.ae.wild.df$POP) )

cnv.ar.selected.df <- gather(cnv.ar.selected,"POP", "AR")
cnv.ar.selected.df$Group <- rep("SELECTED",length(cnv.ar.selected.df$POP) )

cnv.ar.ae.selected.df <- gather(cnv.ar.ae.selected,"POP", "AR")
cnv.ar.ae.selected.df$Group <- rep("SELECTED",length(cnv.ar.ae.selected.df$POP) )

cnv.ar.df <- rbind.data.frame(cnv.ar.wild.df, cnv.ar.selected.df)
cnv.ar.ae.df <- rbind.data.frame(cnv.ar.ae.wild.df, cnv.ar.ae.selected.df)

pop.cnv.bs.df <- data.frame()
pop.cnv.ae.bs.df <- data.frame()

for (i in c("CL","CLP","CS","HC","HC-VA","HI","SL","SM")){
#assign(paste("cnv.bs",i, sep="."), basic.stats(popsub(mask_gen,i), diploid=FALSE))  
bs.temp <- basic.stats(popsub(mask_gen.niw,i), diploid=FALSE)
bs.temp$perloc$Loci <- locNames(popsub(mask_gen.niw,i))
bs.temp.df <- bs.temp$perloc[c("Loci", "Ht", "Hs")]
bs.temp.df$POP <- rep(i,length(bs.temp.df$Loci) )
bs.temp.df$Group <- rep("WILD",length(bs.temp.df$Loci) )
row.names(bs.temp.df) <- NULL

pop.cnv.bs.df <- rbind.data.frame(pop.cnv.bs.df, bs.temp.df)
print(i)
}
```

    ## [1] "CL"
    ## [1] "CLP"
    ## [1] "CS"
    ## [1] "HC"
    ## [1] "HC-VA"
    ## [1] "HI"
    ## [1] "SL"
    ## [1] "SM"

``` r
for (i in c("DEBY","LOLA","NEH", "OBOYS","UMFS")){
#assign(paste("cnv.bs",i, sep="."), basic.stats(popsub(mask_gen,i), diploid=FALSE))  
bs.temp <- basic.stats(popsub(mask_gen.niw,i), diploid=FALSE)
bs.temp$perloc$Loci <- locNames(popsub(mask_gen.niw,i))
bs.temp.df <- bs.temp$perloc[c("Loci", "Ht", "Hs")]
bs.temp.df$POP <- rep(i,length(bs.temp.df$Loci) )
bs.temp.df$Group <- rep("SELECTED",length(bs.temp.df$Loci) )
row.names(bs.temp.df) <- NULL

pop.cnv.bs.df <- rbind.data.frame(pop.cnv.bs.df, bs.temp.df)

}

for (i in c("CLP","CS","HC","HC-VA","HI","SM")){
#assign(paste("cnv.bs",i, sep="."), basic.stats(popsub(mask_gen,i), diploid=FALSE))  
bs.temp <- basic.stats(popsub(mask_gen.ae,i), diploid=FALSE)
bs.temp$perloc$Loci <- locNames(popsub(mask_gen.ae,i))
bs.temp.df <- bs.temp$perloc[c("Loci", "Ht", "Hs")]
bs.temp.df$POP <- rep(i,length(bs.temp.df$Loci) )
bs.temp.df$Group <- rep("WILD",length(bs.temp.df$Loci) )
row.names(bs.temp.df) <- NULL

pop.cnv.ae.bs.df <- rbind.data.frame(pop.cnv.ae.bs.df, bs.temp.df)

}

for (i in c("NEH", "UMFS")){
#assign(paste("cnv.bs",i, sep="."), basic.stats(popsub(mask_gen,i), diploid=FALSE))  
bs.temp <- basic.stats(popsub(mask_gen.ae,i), diploid=FALSE)
bs.temp$perloc$Loci <- locNames(popsub(mask_gen.ae,i))
bs.temp.df <- bs.temp$perloc[c("Loci", "Ht", "Hs")]
bs.temp.df$POP <- rep(i,length(bs.temp.df$Loci) )
bs.temp.df$Group <- rep("SELECTED",length(bs.temp.df$Loci) )
row.names(bs.temp.df) <- NULL

pop.cnv.ae.bs.df <- rbind.data.frame(pop.cnv.ae.bs.df, bs.temp.df)

}
```

``` r
pop.cnv.bs.df$Group <- factor(pop.cnv.bs.df$Group, levels=c("WILD", "SELECTED" ))
pop.cnv.bs.df$POP <- factor(pop.cnv.bs.df$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

cnv.ar.df$Group <- factor(cnv.ar.df$Group, levels=c("WILD", "SELECTED" ))
cnv.ar.df$POP <- factor(cnv.ar.df$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

cnv.total.df <- pop.cnv.bs.df
cnv.total.df$AR <- cnv.ar.df$AR

pop.cnv.ae.bs.df$Group <- factor(pop.cnv.ae.bs.df$Group, levels=c("WILD", "SELECTED" ))
pop.cnv.ae.bs.df$POP <- factor(pop.cnv.ae.bs.df$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"))

cnv.ar.ae.df$Group <- factor(cnv.ar.ae.df$Group, levels=c("WILD", "SELECTED" ))
cnv.ar.ae.df$POP <- factor(cnv.ar.ae.df$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"))

cnv.ae.total.df <- pop.cnv.ae.bs.df
cnv.ae.total.df$AR <- cnv.ar.ae.df$AR


Loc.names.df <- as.data.frame(str_split_fixed(cnv.total.df$Loci, "_", 3))
Loc.names.df%>% 
  mutate(V1 = str_replace(V1, "^1$","NC_035780.1")) %>% 
  mutate(V1 = str_replace(V1, "^2","NC_035781.1")) %>% 
  mutate(V1 = str_replace(V1, "^3","NC_035782.1")) %>% 
  mutate(V1 = str_replace(V1, "^4","NC_035783.1")) %>% 
  mutate(V1 = str_replace(V1, "^5","NC_035784.1")) %>% 
  mutate(V1 = str_replace(V1, "^6","NC_035785.1")) %>% 
  mutate(V1 = str_replace(V1, "^7","NC_035786.1")) %>% 
  mutate(V1 = str_replace(V1, "^8","NC_035787.1")) %>%
  mutate(V1 = str_replace(V1, "^9","NC_035788.1")) %>% 
  mutate(V1 = str_replace(V1, "^10","NC_035789.1"))  -> Loc.names.df
Loc.names.df$V2 <- as.numeric(strtoi(Loc.names.df$V2))

write(paste(Loc.names.df$V1,Loc.names.df$V2 - 1,Loc.names.df$V2,Loc.names.df$V3, cnv.total.df$Ht, cnv.total.df$AR, cnv.total.df$POP, cnv.total.df$Group, sep='\t'),"./Population_Genomics/Intermediate_Files/masked.cnv.div.bed")


Loc.names.ae.df <- as.data.frame(str_split_fixed(cnv.ae.total.df$Loci, "_", 3))
Loc.names.ae.df%>% 
  mutate(V1 = str_replace(V1, "^1$","NC_035780.1")) %>% 
  mutate(V1 = str_replace(V1, "^2","NC_035781.1")) %>% 
  mutate(V1 = str_replace(V1, "^3","NC_035782.1")) %>% 
  mutate(V1 = str_replace(V1, "^4","NC_035783.1")) %>% 
  mutate(V1 = str_replace(V1, "^5","NC_035784.1")) %>% 
  mutate(V1 = str_replace(V1, "^6","NC_035785.1")) %>% 
  mutate(V1 = str_replace(V1, "^7","NC_035786.1")) %>% 
  mutate(V1 = str_replace(V1, "^8","NC_035787.1")) %>%
  mutate(V1 = str_replace(V1, "^9","NC_035788.1")) %>% 
  mutate(V1 = str_replace(V1, "^10","NC_035789.1"))  -> Loc.names.ae.df
Loc.names.ae.df$V2 <- as.numeric(strtoi(Loc.names.ae.df$V2))

write(paste(Loc.names.ae.df$V1,Loc.names.ae.df$V2 - 1,Loc.names.ae.df$V2,Loc.names.ae.df$V3, cnv.ae.total.df$Ht, cnv.ae.total.df$AR, cnv.ae.total.df$POP, cnv.ae.total.df$Group, sep='\t'),"./Population_Genomics/Intermediate_Files/masked.cnv.div.ae.bed")
```

``` bash
for i in {"HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"}
do
mawk -v v=$i '$7 == v' ./Population_Genomics/Intermediate_Files/masked.cnv.div.bed > ./Population_Genomics/Intermediate_Files/$i.masked.cnv.div.bed
done

for i in {"HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"}
do
mawk -v v=$i '$7 == v' ./Population_Genomics/Intermediate_Files/masked.cnv.div.ae.bed > ./Population_Genomics/Intermediate_Files/$i.masked.cnv.div.ae.bed
done
```

``` bash
ls ./Population_Genomics/Intermediate_Files/*.masked.cnv.div.ae.bed| sed 's/.masked.cnv.div.ae.bed//g'  |sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "bedtools intersect -b ./Population_Genomics/Intermediate_Files/{}.masked.cnv.div.ae.bed -a ./other_files/10kb.bed -wa -wb | bedtools merge -d -1 -c 8,9,10,11 -o mean,mean,first,first > ./Population_Genomics/Intermediate_Files/{}.cnv.div.ae.per10kb.bed"

ls ./Population_Genomics/Intermediate_Files/*.masked.cnv.div.bed| sed 's/.masked.cnv.div.bed//g'  |sed 's:\./Population_Genomics/Intermediate_Files/::g' | parallel "bedtools intersect -b ./Population_Genomics/Intermediate_Files/{}.masked.cnv.div.bed -a ./other_files/10kb.bed -wa -wb | bedtools merge -d -1 -c 8,9,10,11 -o mean,mean,first,first > ./Population_Genomics/Intermediate_Files/{}.cnv.div.per10kb.bed"

cd ./Population_Genomics/Intermediate_Files/

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") *.cnv.div.per10kb.bed > ../pop.windowed.cnv
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") *.cnv.div.ae.per10kb.bed > ../pop.windowed.ae.cnv

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") <(grep WILD ../pop.windowed.cnv)  > pop.wild.windowed.cnv
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") <(grep WILD ../pop.windowed.ae.cnv)  > pop.wild.windowed.ae.cnv

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") <(grep SELECTED ../pop.windowed.cnv)  > pop.selected.windowed.cnv
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tPOP\tGROUP") <(grep SELECTED ../pop.windowed.ae.cnv)  > pop.selected.windowed.ae.cnv

bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.cnv| sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.wild.windowed.cnv
bedtools merge -i <(mawk '!/CHR/' pop.wild.windowed.ae.cnv| sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.wild.windowed.ae.cnv

bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.cnv | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.selected.windowed.cnv
bedtools merge -i <(mawk '!/CHR/' pop.selected.windowed.ae.cnv | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.selected.windowed.ae.cnv
bedtools merge -i <(mawk '!/CHR/ && !/LOLA/ && !/DEBY/' pop.selected.windowed.cnv | sort -k1,1 -k2,2n ) -c 4,5,7 -d -1 -o mean,mean,first > ../mean.pop.selected.windowed.nom.cnv

cd .. 

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.wild.windowed.cnv mean.pop.selected.windowed.cnv > mean.pop.windowed.cnv
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.wild.windowed.ae.cnv mean.pop.selected.windowed.ae.cnv > mean.pop.windowed.ae.cnv
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.wild.windowed.cnv mean.pop.selected.windowed.nom.cnv > mean.pop.windowed.nom.cnv

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.wild.windowed.cnv > mean.pop.wild.windowed.cnv.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.selected.windowed.cnv > mean.pop.selected.windowed.cnv.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.wild.windowed.ae.cnv > mean.pop.wild.windowed.ae.cnv.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tDIVERSITY\tAR\tGROUP") mean.pop.selected.windowed.ae.cnv > mean.pop.selected.windowed.ae.cnv.txt

mkdir -p ./Output/CNV_diversity
cp pop.windowed.cnv mean.pop.windowed.cnv mean.pop.windowed.nom.cnv pop.windowed.ae.cnv mean.pop.windowed.ae.cnv ./Output/CNV_diversity
```

``` r
cnv.wvd.dataframe<-read.table(paste(prefix,"pop.windowed.cnv", sep=""), sep="\t", header=T)
cnv.wvd.dataframe$GROUP <- factor(cnv.wvd.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
cnv.wvd.dataframe$POP <- factor(cnv.wvd.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","CL","SL","UMFS","NEH","DEBY","LOLA","OBOYS"))

cnv.wvd.mean.dataframe<-read.table(paste(prefix,"mean.pop.windowed.cnv", sep=""), sep="\t", header=T)
cnv.wvd.mean.dataframe$GROUP <- factor(cnv.wvd.mean.dataframe$GROUP, levels=c("WILD", "SELECTED" ))

cnv.wvd.mean.nom.dataframe<-read.table(paste(prefix,"mean.pop.windowed.nom.cnv", sep=""), sep="\t", header=T)
cnv.wvd.mean.nom.dataframe$GROUP <- factor(cnv.wvd.mean.nom.dataframe$GROUP, levels=c("WILD", "SELECTED" ))

filename <- paste(prefix,"pop.windowed.ae.cnv", sep="")
cnv.wvd.ae.dataframe<-read.table(filename, sep="\t", header=T)
cnv.wvd.ae.dataframe$GROUP <- factor(cnv.wvd.ae.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
cnv.wvd.ae.dataframe$POP <- factor(cnv.wvd.ae.dataframe$POP, levels=c("HI","SM","HC","CS","CLP","HC-VA","UMFS","NEH"))

filename <- paste(prefix,"mean.pop.windowed.ae.cnv", sep="")
cnv.wvd.mean.ae.dataframe<-read.table(filename, sep="\t", header=T)
cnv.wvd.mean.ae.dataframe$GROUP <- factor(cnv.wvd.mean.ae.dataframe$GROUP, levels=c("WILD", "SELECTED" ))
```

``` r
bd.cnv.ar <-ggplot(cnv.wvd.dataframe[which(cnv.wvd.dataframe$AR > 0),], aes(x=AR,y = POP))+
 geom_density_ridges(aes(fill=GROUP), rel_min_height = 0.001,quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), alpha=0.85)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.dataframe$POP)))+
  ggtitle("Wild vs. Selected CNV Allelic Richness")+
  theme_classic() +
  labs(x="CNV Richness") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank()) +
  coord_cartesian(xlim=c(0.5,4.5))
  
bd.cnv.ae.ar <-ggplot(cnv.wvd.ae.dataframe[which(cnv.wvd.ae.dataframe$AR > 0),], aes(x=AR,y = POP))+
  geom_density_ridges(aes(fill=GROUP), rel_min_height = 0.001,quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), alpha=0.85)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.ae.dataframe$POP)))+
  ggtitle("Wild vs. Selected CNV Allelic Richness (Atlantic Samples Only)")+
  theme_classic() +
  labs(x="CNV Richness") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank()) +
  coord_cartesian(xlim=c(0.5,4.5))

bd.cnv.gd <-ggplot(cnv.wvd.dataframe[which(cnv.wvd.dataframe$DIVERSITY > -1),], aes(x=DIVERSITY,y = POP))+
  geom_density_ridges(aes(fill=GROUP), rel_min_height = 0.001,quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), alpha=0.85)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.dataframe$POP)))+
  ggtitle("Wild vs. Selected CNV Diversity")+
  theme_classic() +
  labs(x="CNV Diversity") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank()) 

bd.cnv.ae.gd <-ggplot(cnv.wvd.ae.dataframe[which(cnv.wvd.ae.dataframe$DIVERSITY > -1),], aes(x=DIVERSITY,y = POP))+
  geom_density_ridges(aes(fill=GROUP), rel_min_height = 0.001,quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), alpha=0.85)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.ae.dataframe$POP)))+
  ggtitle("Wild vs. Selected CNV Diversity (Atlantic Samples Only)")+
  theme_classic() +
  labs(x="CNV Diversity") +theme(axis.title.x = element_markdown(),  axis.title.y = element_blank()) 
```

``` r
bc.cnv.ar <-ggplot(cnv.wvd.mean.dataframe, aes(x=AR,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.mean.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. Selected CNV Allelic Richness")+
  theme_classic() +
  #scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  #coord_cartesian(xlim=c(0,0.5))+
  labs(x="CNV Richness") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())

bc.cnv.ae.ar <-ggplot(cnv.wvd.mean.ae.dataframe, aes(x=AR,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.mean.ae.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. Selected CNV Allelic Richness (Atlantic Samples Only)")+
  theme_classic() +
  #scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  #coord_cartesian(xlim=c(0,0.5))+
  labs(x="CNV Richness") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())
```

``` r
bc.cnv.gd <-ggplot(cnv.wvd.mean.dataframe, aes(x=DIVERSITY,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.mean.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. CNV Diversity")+
  theme_classic() +
  #scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  #coord_cartesian(xlim=c(0,0.5))+
  labs(x="Observed CNV Diversity") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())

bc.cnv.ae.gd <-ggplot(cnv.wvd.mean.ae.dataframe, aes(x=DIVERSITY,y = GROUP))+
  geom_violin(aes(color=GROUP,fill= GROUP),trim=TRUE) +
  geom_boxplot(aes(fill=GROUP), width=0.1,outlier.shape = NA, outlier.color = "black")+
  stat_summary(fun=mean, geom="point", shape=23, size=2)+
  scale_fill_manual(values=c("#D55E00", "#009E73"))+
  scale_color_manual(values=c("#D55E00", "#009E73"))+
  scale_y_discrete(limits = rev(levels(cnv.wvd.mean.ae.dataframe$GROUP)))+
  #ylab("Origin")+
  ggtitle("Wild vs. CNV Diversity (Atlantic Samples Only)")+
  theme_classic() +
  #scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=c("0.0","0.1","0.2","0.3","0.4","0.5*"))+
  #coord_cartesian(xlim=c(0,0.5))+
  labs(x="Observed CNV Diversity") +
  theme(legend.position = c(0.85, 0.5), axis.title.y = element_blank(),axis.text.y=element_blank(), plot.title = element_markdown(),axis.title.x = element_markdown())
```

``` r
cnv.plot <- function(wild.filename,selected.filename){
wcnv <- read.table(wild.filename, header = TRUE)
scnv <- read.table(selected.filename, header = TRUE)

wcnv <- wcnv %>% dplyr::rename(WILD_DIVERSITY= DIVERSITY)
scnv <- scnv %>% dplyr::rename(SELECTED_DIVERSITY= DIVERSITY)
wcnv <- wcnv %>% dplyr::rename(WILD_AR= AR)
scnv <- scnv %>% dplyr::rename(SELECTED_AR= AR)

w.s.cnv <- join(wcnv,scnv, by =c("CHROM", "BIN_START"),type="full")

w.s.cnv$DIVERSITY_DIFF <- wcnv$WILD_DIVERSITY - scnv$SELECTED_DIVERSITY
w.s.cnv$AR_DIFF <- wcnv$WILD_AR - scnv$SELECTED_AR

SNP<-c(1:(nrow(w.s.cnv)))
mydf<-data.frame(SNP,w.s.cnv)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf
mydf$CHR <- as.numeric(mydf$CHR)  


  plot.val <- function(stat, subset2){
    if(subset2 == "WILD"){
      pre <- "m2"
      title <-  "Wild"
      colpal <- c("#D55E00", "#8f3e00")
      }
    else if (subset2 =="SELECTED"){
      pre <- "m3"
      title <-  "Selected"
      colpal <- c("#009E73", "#00543d")
    }
    else {
      pre <- "md4"
      title <-  "Wild - Domestic"
      colpal <- c("#7d5700", "#634500")
    }
  
    
    if(stat=="AR") {
        ylabel <- paste("CNV","Richness", sep=" ")
        ymax <- 6
    } 
    else{
        ylabel <- paste("CNV","Diversity", sep=" ")
        ymax <- 1
    }
    if(subset2 == "DIFF"){
      ylabel <- paste("Difference in", ylabel, sep =" ")
       p <- paste(stat,subset2,sep="_")
    }
    else{ p <- paste(subset2,stat,sep="_")}
    
    m2 <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue=p,snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)

    dfm <- m2[[1]]
    dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

    dfmsplit <- split(dfm, dfm$chrom)

      xbreaks <- sapply(dfmsplit,function(x) {
      midpoint <- length(x$index)/2
      if(midpoint <1) midpoint <- 1
      return(x$index[midpoint])
      })

      m2.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
        geom_point(size=0.5, aes(alpha=0.33))+  
        scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = 
                             c(0,max(dfm$index)+10)) +
        guides(colour = "none",alpha= "none") +
        labs(x = "Chromosome") +
        ggtitle(title) +theme_classic()+
        labs(y=ylabel) +theme(axis.title.y = element_markdown())+ 
        scale_color_manual(values=colpal) 
      
      if(stat == "AR"){
        m2.total <- m2.total + scale_y_sqrt(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),ymax))
        post <- "ar"
        }
        else{m2.total <- m2.total +scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker-0.001,na.rm=TRUE),ymax))
        post="div"
        }
      if( p == "AR_DIFF" ){
        m2.total <- m2.total + scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker)-0.1,max(dfm$marker)+0.1))
        post <- "ar"
      }
      if( p == "DIVERSITY_DIFF" ){
        m2.total <- m2.total + scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker)-0.1,max(dfm$marker)+0.1))
        post <- "div"
      }
      
      fig.title <- paste(pre,post,sep = ".")
      assign(fig.title, m2.total, envir = .GlobalEnv)

  }
  
  plot.val("AR","WILD")
  plot.val("AR","SELECTED")
  plot.val("AR","DIFF")
  plot.val("DIVERSITY","WILD")
  plot.val("DIVERSITY","SELECTED")
  plot.val("DIVERSITY","DIFF")

}
```

``` r
cnv.plot(paste(prefix,"mean.pop.wild.windowed.cnv.txt", sep=""),paste(prefix,"mean.pop.selected.windowed.cnv.txt", sep=""))
```

``` r
png(filename=paste(prefix.OFS,"Figure.S6.CNV.Richness.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AADD
BBDD
CCDD
EEEE
FFFF
"
m2.ar+m3.ar+md4.ar+bd.cnv.ar+bc.cnv.ar+bc.cnv.ae.ar +plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')
```

    ## Picking joint bandwidth of 0.0877

``` r
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OFS,"Figure.S7.CNV.Diversity.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AADD
BBDD
CCDD
EEEE
FFFF
"
m2.div+m3.div+md4.div+bd.cnv.gd+bc.cnv.gd+bc.cnv.ae.gd +plot_layout(design=layout, guides = "collect")+ plot_annotation(tag_levels = 'A')
```

    ## Picking joint bandwidth of 0.0345

``` r
dev.off()
```

    ## png 
    ##   2

``` r
filenamew <- paste(prefix,"mean.pop.wild.windowed.ae.cnv.txt", sep="")
filenames <- paste(prefix,"mean.pop.selected.windowed.ae.cnv.txt", sep="")

cnv.plot(filenamew,filenames)
```

``` r
png(filename=paste(prefix.OFS,"Figure.S8.CNV.ae.Richness.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AADD
BBDD
CCDD
EEDD
EEDD
EEDD
"
m2.ar+m3.ar+md4.ar+bd.cnv.ae.ar+bc.cnv.ae.ar +plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')
```

    ## Picking joint bandwidth of 0.0641

``` r
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OFS,"Figure.S9.CNV.ae.Diversity.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "
AADD
BBDD
CCDD
EEDD
EEDD
EEDD
"
m2.div+m3.div+md4.div+bd.cnv.ae.gd+bc.cnv.ae.gd +plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')
```

    ## Picking joint bandwidth of 0.0344

``` r
dev.off()
```

    ## png 
    ##   2

``` r
cnv.wvd.mean.dataframe %>% t_test(formula = AR ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     0.389 19919.   0.697 two.sided    0.00229 -0.00924   0.0138

``` r
cnv.wvd.mean.dataframe %>% t_test(formula = AR ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     0.389 19919.   0.651 less         0.00229     -Inf   0.0120

``` r
group_by(cnv.wvd.mean.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(AR, na.rm = TRUE),
    sd = sd(AR, na.rm = TRUE),
    se=std.error(AR, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD     10001  1.85 0.402 0.00402
    ## 2 SELECTED 10001  1.84 0.429 0.00429

``` r
cnv.wvd.df <- cnv.wvd.dataframe[which(cnv.wvd.dataframe$POP != "LOLA" & cnv.wvd.dataframe$POP != "DEBY"),]


group_by(cnv.wvd.mean.nom.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(AR, na.rm = TRUE),
    sd = sd(AR, na.rm = TRUE),
    se=std.error(AR, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD     10001  1.85 0.402 0.00402
    ## 2 SELECTED 10001  1.82 0.469 0.00469

``` r
cnv.wvd.mean.nom.dataframe %>% t_test(formula = AR ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df     p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>       <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      5.01 19549. 0.000000544 two.sided     0.0310   0.0188   0.0431

``` r
group_by(cnv.wvd.mean.ae.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(AR, na.rm = TRUE),
    sd = sd(AR, na.rm = TRUE),
    se=std.error(AR, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD      9218  1.91 0.404 0.00421
    ## 2 SELECTED  9218  1.88 0.518 0.00539

``` r
cnv.wvd.mean.ae.dataframe %>% t_test(formula = AR ~ GROUP, alternative = "greater",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df   p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>     <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      4.07 17414. 0.0000232 greater       0.0279   0.0166      Inf

``` r
cnv.wvd.mean.dataframe %>% t_test(formula = DIVERSITY ~ GROUP, alternative = "two.sided",order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -2.45 19941.  0.0145 two.sided   -0.00528 -0.00950 -0.00105

``` r
cnv.wvd.mean.dataframe %>% t_test(formula = DIVERSITY ~ GROUP, alternative = "less", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1     -2.45 19941. 0.00723 less        -0.00528     -Inf -0.00173

``` r
group_by(cnv.wvd.mean.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(DIVERSITY, na.rm = TRUE),
    sd = sd(DIVERSITY, na.rm = TRUE),
    se=std.error(DIVERSITY, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD     10001 0.324 0.148 0.00148
    ## 2 SELECTED 10001 0.329 0.157 0.00157

``` r
cnv.wvd.mean.nom.dataframe %>% t_test(formula = DIVERSITY ~ GROUP, alternative = "greater", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      2.04 19580.  0.0208 greater      0.00462 0.000890      Inf

``` r
group_by(cnv.wvd.mean.nom.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(DIVERSITY, na.rm = TRUE),
    sd = sd(DIVERSITY, na.rm = TRUE),
    se=std.error(DIVERSITY, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD     10001 0.324 0.148 0.00148
    ## 2 SELECTED 10001 0.319 0.172 0.00172

``` r
cnv.wvd.mean.ae.dataframe %>% t_test(formula = DIVERSITY ~ GROUP, alternative = "two.sided", order = c("WILD", "SELECTED"))
```

    ## # A tibble: 1 × 7
    ##   statistic   t_df p_value alternative estimate lower_ci upper_ci
    ##       <dbl>  <dbl>   <dbl> <chr>          <dbl>    <dbl>    <dbl>
    ## 1      1.35 17396.   0.176 two.sided    0.00340 -0.00153  0.00832

``` r
group_by(cnv.wvd.mean.ae.dataframe, GROUP) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(DIVERSITY, na.rm = TRUE),
    sd = sd(DIVERSITY, na.rm = TRUE),
    se=std.error(DIVERSITY, na.rm = TRUE) )
```

    ## # A tibble: 2 × 5
    ##   GROUP    count  mean    sd      se
    ##   <fct>    <int> <dbl> <dbl>   <dbl>
    ## 1 WILD      9218 0.347 0.148 0.00154
    ## 2 SELECTED  9218 0.344 0.190 0.00198

## Detecting Large Inversions

### Run PCAdapt on full dataset

``` r
path_to_file <- "./Population_Genomics/Masked.NoInbred.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
```

``` r
x <- pcadapt(input = filename, K = 10)
```

``` r
plot(x , option = "manhattan")
```

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->

``` r
filename <- paste(prefix.I,"NoInbredLociPCAadpt.pvalues.txt", sep="")
write(paste(CHR.m.ni,POS.m.ni, x$pvalues, sep='\t'),filename)
```

``` r
filename <- paste(prefix.I,"NoInbredLociPCAadpt.qvalues.txt", sep="")
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.0001
outliers <- which(qval < alpha)
length(outliers)
```

    ## [1] 111367

``` r
write(paste(CHR.m.ni,POS.m.ni, qval, sep='\t'),filename)
```

``` bash
mawk '{print $1 "\t" $2-1 "\t" $2 "\t" $3 }' ./Population_Genomics/Intermediate_Files/NoInbredLociPCAadpt.pvalues.txt | sort -k1,1 -k2,2n --parallel=40 > ./Population_Genomics/Intermediate_Files/sorted.NoInbredPCAdapt.bed

bedtools intersect -b ./Population_Genomics/Intermediate_Files/sorted.NoInbredPCAdapt.bed -a <(sort -k 1,1 -k2,2n ./other_files/10kb.bed) -wa -wb -sorted | bedtools merge -i - -d -1 -c 7 -o mean > ./Population_Genomics/Intermediate_Files/10kb.mean.pcadapt.bed

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tP_VALUE") ./Population_Genomics/Intermediate_Files/10kb.mean.pcadapt.bed  > ./Population_Genomics/Intermediate_Files/10kb.mean.pcadapt.bed.txt

mkdir -p ./Population_Genomics/Output/pcadapt
cp ./Population_Genomics/Intermediate_Files/sorted.NoInbredPCAdapt.bed ./Population_Genomics/Intermediate_Files/10kb.mean.pcadapt.bed.txt ./Population_Genomics/Output/pcadapt
```

### For finding bounds of inversion

``` r
filename <- paste(prefix.I,"10kb.mean.pcadapt.bed.txt", sep="")
pv <- read.table(filename, header=TRUE)

w.d.pi$PI_DIFF <- w.d.pi$DIFF

w.d.het$HET_DIFF <- w.d.het$DIFF

join1<- join(w.d.tajd,w.d.het, by=c("CHROM", "BIN_START"),type="full")
join2 <- join(join1,w.d.pi, by=c("CHROM", "BIN_START"),type="full")
join3 <- join(join2,pv, by=c("CHROM", "BIN_START"),type="full")



pal1 <- wes_palette("IsleofDogs1")
pal2 <- wes_palette("Zissou1")
pal <- c(pal1[1],pal2[1])

SNP <-c(1:(nrow(join3)))
mydf<-data.frame(SNP,join3)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf

mydf$CHR <- as.numeric(mydf$CHR)  

md4.pi <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="PI_DIFF",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "PI_DIFF", :
    ## NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.pi[[1]]
dfm$PI_DIFF[is.na(dfm$PI_DIFF)] <- 0
dfm$HET_DIFF[is.na(dfm$HET_DIFF)] <- 0
#dfm$P_VALUE[is.na(dfm$P_VALUE)] <- 1
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))


Chrom_Calc <- function(dfm,chrom) {
  dfm.sub <- dfm[which(dfm$CHR == chrom),]
  dfm.sub <- dfm.sub %>% dplyr::mutate(mean_pi100 = zoo::rollmean(PI_DIFF, k = 100,fill = c(0,0,0), align="center"), 
                                       mean_pi500 =zoo::rollmean(PI_DIFF, k = 500,fill = c(0,0,0), align="center"),
                                       mean_het100 = zoo::rollmean(HET_DIFF, k = 100, fill = c(0,0,0),alight="center"),
                                       mean_het500 = zoo::rollmean(HET_DIFF, k = 500,fill = c(0,0,0), align="center"),
                                       mean_pv100 = zoo::rollmean(-log10(P_VALUE), k = 100,fill = NA, align="center"),  
                                       mean_pv500 = zoo::rollmean(-log10(P_VALUE), k = 500,fill = NA, align="center"), 
                                       mean_pv50 = zoo::rollmean(-log10(P_VALUE), k = 50,fill = NA, align="center"))
  #assign(paste("dfm.chr.",chrom,sep="")) <- dfm.sub
  return(dfm.sub)
}
dfm.subs <- list()

for(i in 1:10){
  dfm.subs[[i]] <- Chrom_Calc(dfm, chroms[i])
  
}

dfm <- bind_rows(dfm.subs)


### For finding bounds of inversion ###
quant=0.1
quant2=0.1
quant3=0.01
qpi100 <- quantile(dfm$mean_pi100, quant2)
qpi500 <-quantile(dfm$mean_pi500 , quant)
qhet100 <-quantile(dfm$mean_het100, quant2)
qhet500 <-quantile(dfm$mean_het500, quant)
qpv100 <-quantile(dfm$mean_pv100, 1-quant2, na.rm = TRUE)
qpv500 <-quantile(dfm$mean_pv500, 1-quant, na.rm = TRUE)
qpv50 <-quantile(dfm$mean_pv50, 1-quant2, na.rm = TRUE)


dfm.sub1a <- dfm[
  which( dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100 & dfm$mean_het500 < qhet500
       | dfm$mean_pi500 < qpi500 & dfm$mean_pi100  < qpi100  & dfm$mean_het500 < qhet500 
       | dfm$mean_pi100 < qpi100 & dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100  & dfm$mean_het500 < qhet500 
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv50 > qpv50
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv500 > qpv500
       | dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi100 < qpi100 & dfm$mean_het500 < qhet500 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi100 < qpi100 & dfm$mean_het100 < qhet100 &  dfm$mean_pv100 > qpv100  
       | dfm$mean_pv50 > qpv50 &  dfm$mean_pv100 > qpv100 & dfm$mean_pi100 < qpi100 & dfm$mean_het100 < qhet100),]

summary(dfm.sub1a$CHROM)    
```

    ## NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 
    ##         251           0           0           0        1464         732 
    ## NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1 
    ##           0           0           0           0

``` r
dfm.sub1 <- dfm[
  which( dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 
       | dfm$mean_pi100 < qpi100 & dfm$mean_het100 < qhet100 & dfm$CHR == 6    
       | dfm$mean_pv50 > qpv50 &  dfm$mean_pv100 > qpv100 &  dfm$mean_pv500 > qpv500   
       | dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100 & dfm$mean_het500 < qhet500
       | dfm$mean_pi500 < qpi500 & dfm$mean_pi100  < qpi100  & dfm$mean_het500 < qhet500 
       | dfm$mean_pi100 < qpi100 & dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100  & dfm$mean_het500 < qhet500 
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv50 > qpv50
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi500 < qpi500 & dfm$mean_het500 < qhet500 &  dfm$mean_pv500 > qpv500
       | dfm$mean_pi500 < qpi500 & dfm$mean_het100 < qhet100 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi100 < qpi100 & dfm$mean_het500 < qhet500 &  dfm$mean_pv100 > qpv100
       | dfm$mean_pi100 < qpi100 & dfm$mean_het100 < qhet100 &  dfm$mean_pv100 > qpv100  
       | dfm$mean_pv50 > qpv50 &  dfm$mean_pv100 > qpv100 & dfm$mean_pi100 < qpi100 & dfm$mean_het100 < qhet100),]


dfm.sub1 <- dfm.sub1[which(dfm.sub1$CHR == 1 | dfm.sub1$CHR == 6 | dfm.sub1$CHR == 5),]

summary(dfm.sub1$CHROM)    
```

    ## NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 
    ##         380           0           0           0        1498         944 
    ## NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1 
    ##           0           0           0           0

``` r
dfm.plot <- anti_join(dfm, dfm.sub1, by=c("CHROM", "BIN_START"))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

md4.pi <- ggplot(dfm.plot, aes(x= index, y=marker, colour = as.factor(chrom_alt), fill = as.factor(chrom_alt)))+
  geom_point(alpha=0.35,shape=19, size=0.5)+  
  #geom_point(aes(y=mean_pi100,x=index), color = "green",size=0.75)+
  #geom_point(aes(y=mean_pi500,x=index), color = "yellow",size=0.75)+
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.006,0.007))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") + ggtitle("Wild - Domestic") +
  theme_classic()+
  labs(y="&Delta; &pi;") +theme(axis.title.y = element_markdown(), axis.title.x = element_blank())+
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=c(pal[1],pal[2]))


md4.pi <-md4.pi+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black", aes( fill = as.factor(chrom_alt)))


md4.het <- ggplot(dfm.plot, aes(x= index, y=HET_DIFF, colour = as.factor(chrom_alt)))+
                    #size=(0.1+ abs(marker))^2)
  geom_point(alpha=0.25, size =0.5)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.4,0.4))+
  labs(x = "Chromosome") + 
  theme_classic()+
  labs(y="&Delta; Heterozygosity") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=c(pal[1],pal[2]))

md4.het <-md4.het+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black", aes( fill = as.factor(chrom_alt)))



mp.total <- ggplot(dfm.plot, aes(x= index, y=-log10(P_VALUE), colour = as.factor(chrom_alt)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker),3))+
  theme_classic()+
  labs(y="-log(*P-value*)") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=c(pal[1],pal[2]))


mp.total <-mp.total+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black", aes( fill = as.factor(chrom_alt)))
```

``` r
mp.total / md4.pi / md4.het 
```

    ## Warning: Removed 910 rows containing missing values (geom_point).

    ## Warning: Removed 38 rows containing missing values (geom_point).

    ## Warning: Removed 117 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_point).

![](Oyster_Genome_Population_Genomic_Analysis_files/figure-gfm/unnamed-chunk-105-1.png)<!-- -->

``` r
## Upon visual inspection the boundary on Chromosome six was extended to 44000000
        
dfm.sub1.chr1 <- dfm.sub1[which(dfm.sub1$CHR == 1),]
dfm.sub1.chr5 <- dfm.sub1[which(dfm.sub1$CHR == 5),]
dfm.sub1.chr6 <- dfm.sub1[which(dfm.sub1$CHR == 6),]
#dfm.sub1.chr6 <- rbind(dfm.sub1.chr6,c(0,"NC_035785.1",43870000,44000000,0,0,0))


bed.df <- rbind(c(as.character(dfm.sub1.chr1$CHROM[1]),min(dfm.sub1.chr1$BIN_START),max(dfm.sub1.chr1$BIN_END)),c(as.character(dfm.sub1.chr5$CHROM[1]),min(dfm.sub1.chr5$BIN_START),max(dfm.sub1.chr5$BIN_END)),c(as.character(dfm.sub1.chr6$CHROM[1]),min(dfm.sub1.chr6$BIN_START),max(dfm.sub1.chr6$BIN_END,"44000000")))

write(paste(bed.df[,1],bed.df[,2], bed.df[,3], sep='\t'),"./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed")
```

``` bash
mkdir -p ./Population_Genomics/Output/structural_variants
cp ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed ./Population_Genomics/Output/structural_variants
```

# Exploring Genomic Outliers within the context of genomic architecture

## Signal of differentiation across large inversions

``` bash
#bedtools intersect -a mean.wild.dome.diff.sites.pi -b Detected_large_inversions.bed | mawk '$4 < 0' > #big.inv.pi.site.diff.bed
#
#bedtools intersect -a mean.wild.dome.diff.sites.het -b Detected_large_inversions.bed | mawk '$4 < 0' > #big.inv.het.site.diff.bed
#
#
#cat big.inv.het.site.diff.bed big.inv.pi.site.diff.bed | sort -k1,1 -k2,2n | bedtools merge -i - >  big.inv.pi.het.site.bed
#
#bedtools intersect -a Masked_NoInbred_FST_Outflank.bed -b big.inv.pi.het.site.bed | mawk '$5 ~/TR/' > #inversion.fst.pi.het.outliers.bed

bedtools intersect -a ./Population_Genomics/Masked_NoInbred_FST_Outflank.bed -b ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed | mawk '$5 ~/TR/' > ./Population_Genomics/Intermediate_Files/inversion.fst.large_inversion.outliers.bed

bcftools view --threads 40 ./Population_Genomics/SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz -R ./Population_Genomics/Intermediate_Files/inversion.fst.large_inversion.outliers.bed -O z -o ./Population_Genomics/Large.inversion.outliers.vcf.gz

bcftools index ./Population_Genomics/Large.inversion.outliers.vcf.gz

mkdir -p ./Population_Genomics/Output/Outliers
cp ./Population_Genomics/Large.inversion.outliers.vcf.gz ./Population_Genomics/Intermediate_Files/inversion.fst.large_inversion.outliers.bed ./Population_Genomics/Output/Outliers
```

``` r
source("./scripts/genotype_plot.R")
vcf.inv <- read.vcfR( "./Population_Genomics/Large.inversion.outliers.vcf.gz", verbose = FALSE )

gid.inv <- vcfR2genind(vcf.inv, return.alleles = TRUE) 

pop(gid.inv) <- popmap.ni$POP
```

``` r
vcf_file = "./Population_Genomics/Large.inversion.outliers.vcf.gz"

inv.gp <- genotype_plot(vcf=vcf_file,chr= "NC_035780.1",popmap = gp.popmap, start=1, end=94778955,cluster=FALSE,colour_scheme=c(col_pal[1],wes_palette("IsleofDogs1"))) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 138
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 138
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 138
    ##   row_num: 0
    ## Processed variant: 138
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
dfx <- min(ggplot_build(inv.gp$genotypes)$data[[1]]$x) 
dfmax <- max(ggplot_build(inv.gp$genotypes)$data[[1]]$x)

ddelta <- (dfmax-dfx)/20
dfx <- dfx -ddelta

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
```

    ## Joining by: IND

``` r
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)

gp.plot.inv.1.out <- inv.gp$genotypes + geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X) -ddelta,max(ggplot_build(inv.gp$genotypes)$data[[1]]$x) +1 )) + scale_color_manual(values=p_noinbred , name="Population/Line") +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_text(size=12))+
scale_shape_manual(values=c(18,16),name="Origin") +
scale_size_manual(values=c(3.5,3),name="Origin", guide="none") +
guides(fill = guide_legend(override.aes = list(size=6)))+ 
ggtitle("Chr. 1") +
theme(plot.title = element_text(margin=margin(0,0,-40,20))) +
#guides(color=guide_legend(override.aes = list(size=6, shape=c(16,16,18,18,16,16,18,16,16,18,16,16,18))))+
#guides(shape = guide_legend(override.aes = list(size=6)))
guides(color="none")+
guides(shape = "none")  
```

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
inv.gp <- genotype_plot(vcf=vcf_file,chr= "NC_035784.1",popmap = gp.popmap, start=1, end=64020000,cluster=FALSE,colour_scheme=c(col_pal[5],wes_palette("IsleofDogs1"))) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 522
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 522
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 522
    ##   row_num: 0
    ## Processed variant: 522
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
dfx <- min(ggplot_build(inv.gp$genotypes)$data[[1]]$x) 
dfmax <- max(ggplot_build(inv.gp$genotypes)$data[[1]]$x)

ddelta <- (dfmax-dfx)/20
dfx <- dfx -ddelta

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
```

    ## Joining by: IND

``` r
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)


gp.plot.inv.5.1.out <- inv.gp$genotypes + geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X) -ddelta,max(ggplot_build(inv.gp$genotypes)$data[[1]]$x) +1 )) + scale_color_manual(values=p_noinbred , name="Population/Line") +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y=element_blank())+
scale_shape_manual(values=c(18,16),name="Origin") +
  ggtitle("Chr. 5") +
scale_size_manual(values=c(3.5,3),name="Origin", guide="none") +
  theme(plot.title = element_text(margin=margin(0,0,-40,20))) +
guides(fill = "none")+ 
guides(color="none")+
guides(shape = "none")
```

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
inv.gp <- genotype_plot(vcf=vcf_file,chr= "NC_035784.1",popmap = gp.popmap, start=64030000, end=81090000,cluster=FALSE,colour_scheme=c(col_pal[5],wes_palette("IsleofDogs1"))) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 627
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 627
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 627
    ##   row_num: 0
    ## Processed variant: 627
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
dfx <- min(ggplot_build(inv.gp$genotypes)$data[[1]]$x) 
dfmax <- max(ggplot_build(inv.gp$genotypes)$data[[1]]$x)

ddelta <- (dfmax-dfx)/20
dfx <- dfx -ddelta

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
```

    ## Joining by: IND

``` r
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)

gp.plot.inv.5.2.out <- inv.gp$genotypes + geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X) -ddelta,max(ggplot_build(inv.gp$genotypes)$data[[1]]$x) +1 )) + scale_color_manual(values=p_noinbred , name="Population/Line") +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y=element_blank())+
scale_shape_manual(values=c(18,16),name="Origin") +
ggtitle("Chr. 5 (6)") +
scale_size_manual(values=c(3.5,3),name="Origin", guide="none") +
  theme(plot.title = element_text(margin=margin(0,0,-40,20))) +
guides(fill = "none")+ 
guides(color="none")+
guides(shape = "none")
```

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
inv.gp <- genotype_plot(vcf=vcf_file,chr= "NC_035785.1",popmap = gp.popmap, start=1, end=100000000,cluster=FALSE,colour_scheme=c(col_pal[6],wes_palette("IsleofDogs1"))) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 227
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 227
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 227
    ##   row_num: 0
    ## Processed variant: 227
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
dfx <- min(ggplot_build(inv.gp$genotypes)$data[[1]]$x) 
dfmax <- max(ggplot_build(inv.gp$genotypes)$data[[1]]$x)

ddelta <- (dfmax-dfx)/20
dfx <- dfx -ddelta

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
```

    ## Joining by: IND

``` r
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)

gp.plot.inv.6.1.out <- inv.gp$genotypes + geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X) -ddelta,max(ggplot_build(inv.gp$genotypes)$data[[1]]$x) +1 )) + scale_color_manual(values=p_noinbred , name="Population/Line") +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y=element_blank())+
scale_shape_manual(values=c(18,16),name="Origin") +
ggtitle("Chr. 6") +
scale_size_manual(values=c(3.5,3),name="Origin", guide="none") +
theme(plot.title = element_text(margin=margin(0,0,-40,20))) +
guides(fill = "none")+ 
guides(color="none")+
guides(shape = "none")
```

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

``` r
pca.inv <- dudi.pca(gid.inv ,scannf=FALSE,scale=FALSE,nf=2)
pca.inv.p.df <- pca.inv$li
pca.inv.p.df$POP <- popmap.ni$POP
pca.inv.p.df$TYPE <- popmap.ni$Type


pca.inv.plot <- ggplot(pca.inv.p.df, aes(x=Axis1,y=Axis2,color=POP,fill=POP,shape=TYPE, size=3, alpha=0.95)) + geom_point(alpha=0.95, size = 4)+
scale_shape_manual(values=c(23,21),name="Origin") +
scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") +
xlab("PC1")+ ylab("PC2") +
guides(fill = guide_legend(override.aes = list(size =4, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =4)))+ 
#xlim(-110,65)+ ylim(-40,90) +

geom_dl(alpha=1,method = list(cex = 1.25, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=12)) + theme(legend.position = "right")

png(filename=paste(prefix.OE,"PCA.inv.1.2.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca.inv.plot
dev.off()
```

    ## png 
    ##   2

``` r
fst<-read.table("./Population_Genomics/masked.OF_fst.10kb.txt", header=TRUE)

join4 <- join(join3, fst,by=c("CHROM", "BIN_START"),type="full")


pal <- c(pal3,pal1[1],pal2[1],pal3)
SNP <-c(1:(nrow(join4)))
mydf<-data.frame(SNP,join4)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf
mydf$CHR <- as.numeric(mydf$CHR)  



md4.pi <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="PI_DIFF",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue = "PI_DIFF", :
    ## NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.pi[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))


dfm.sub1 <- dfm[which(dfm$CHR == 1 & dfm$BIN_START >= 40520000 & dfm$BIN_START <= 41910000 |dfm$CHR == 5 & dfm$BIN_START >= 61760000 & dfm$BIN_START <= 80090000 | dfm$CHR == 6 & dfm$BIN_START >= 30520000 & dfm$BIN_START <= 44080000),]

#dfm.sub1 <- dfm[which(dfm$CHR == 5 & dfm$BIN_START > 61450000 & dfm$BIN_START < 80090000 | dfm$CHR == 6 & dfm$BIN_START > 30440000 & dfm$BIN_START < 44080000),]

dfm.sub1 <- dfm.sub1[which(dfm.sub1$PI_DIFF < 0 & dfm.sub1$HET_DIFF < 0), ]

dfm.plot <- anti_join(dfm, dfm.sub1, by=c("CHROM", "BIN_START"))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

md4.pi <- ggplot(dfm.plot, aes(x= index, y=marker, colour = as.factor(chrom_alt), fill = as.factor(CHR)))+
  geom_point(alpha=0.35,shape=19, size=0.5)+  
  #geom_point(aes(y=mean_pi100,x=index), color = "green",size=0.75)+
  #geom_point(aes(y=mean_pi500,x=index), color = "yellow",size=0.75)+
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.006,0.007))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") + ggtitle("Wild - Domestic") +
  theme_classic()+
  labs(y=expression(Delta*" "*pi)) +theme(axis.title.x = element_blank())+
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal)


md4.pi <-md4.pi+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black", aes( fill = as.factor(CHR)))


md4.tajd <- ggplot(dfm.plot, aes(x= index, y=TAJD_DIFF, colour = as.factor(chrom_alt),size=(0.1+ abs(marker))^2))+
  geom_point(alpha=0.85)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE,size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(-3.3,2.5))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") + ggtitle("Wild - Domestic") +
  theme_classic()+
  labs(y="Difference in Tajima's *D*") +theme(axis.title.y = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#7d5700", "#634500"))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
md4.het <- ggplot(dfm.plot, aes(x= index, y=HET_DIFF, colour = as.factor(chrom_alt), fill = as.factor(CHR)))+
                    #size=(0.1+ abs(marker))^2)
  geom_point(alpha=0.25, size =0.5)+  
  scale_size_continuous(range=c(0.01,1.5))+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.4,0.4))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  labs(x = "Chromosome") + 
  #ggtitle("Wild - Domestic") +
  theme_classic()+
  labs(y="&Delta; Heterozygosity") +theme(axis.title.y = element_markdown())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal)

#dfm.sub1 <- dfm[which(dfm$CHR == 5 & dfm$BIN_START > 60600000 & dfm$BIN_START < 80200000 | dfm$CHR == 6 & dfm$BIN_START > 29900000 & dfm$BIN_START < 44500000),]


md4.het <-md4.het+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black")




m2.m <- ggplot(dfm.plot, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt), fill = as.factor(CHR)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") +ggtitle("*F<sub>ST</sub>*") +
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal)


m2o.m <-m2.m+geom_point(data=dfm.sub1[which(dfm.sub1$OUTLIER == TRUE),],alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black")


mp.total <- ggplot(dfm.plot, aes(x= index, y=-log10(P_VALUE), colour = as.factor(chrom_alt), fill = as.factor(CHR)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker),2))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") + ggtitle("*V<sub>ST</sub>*") +
  theme_classic()+
  labs(y="-log(*P-value*)") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  #ylab(*bquote(V[ST])*)+
  #theme(axis.title.y = element_text(face = "italic"))+
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal)



mp.total <-mp.total+geom_point(data=dfm.sub1,alpha=0.8,size=1.5,shape=21, stroke = 0.01, color="black")


layout <- "
AAAAAAAAFFGGHHIIKK
AAAAAAAAFFGGHHIIKK
BBBBBBBBFFGGHHIIKK
BBBBBBBBFFGGHHIIKK
CCCCCCCCFFGGHHIIKK
CCCCCCCCFFGGHHIIKK
DDDDDDDDFFGGHHIIKK
DDDDDDDDFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK
JJJJEEEEFFGGHHIIKK"

layout <- "
AAAAAAAAGGHHIIJJKK
AAAAAAAAGGHHIIJJKK
BBBBBBBBGGHHIIJJKK
BBBBBBBBGGHHIIJJKK
CCCCCCCCGGHHIIJJKK
CCCCCCCCGGHHIIJJKK
DDDDDDDDGGHHIIJJKK
DDDDDDDDGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK
EEEEFFFFGGHHIIJJKK"

filename <- paste(prefix.OF,"Figure.3.Large.Inversions.png", sep="")

png(filename=filename, type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

m2o.m+(mp.total+ plot_layout(tag_level = "new"))+(md4.pi+ plot_layout(tag_level = "new"))+ (md4.het + plot_layout(tag_level = "new")) + finished_mapnolm_small + pca.inv.plot +gp.plot.inv.1.out + (gp.plot.inv.5.1.out+ plot_layout(tag_level = "new")) + (gp.plot.inv.5.2.out+ plot_layout(tag_level = "new")) + (gp.plot.inv.6.1.out+ plot_layout(tag_level = "new"))   + guide_area() + plot_layout(design=layout, guides = "collect") + theme(legend.text = element_text(size=10)) + plot_annotation(tag_levels = 'A') 
```

    ## Warning: Removed 1063 rows containing missing values (geom_point).

    ## Warning: Removed 976 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 122 rows containing missing values (geom_point).

    ## Warning: Removed 921 rows containing missing values (geom_point).

``` r
dev.off()
```

    ## png 
    ##   2

## Outliers outside of large Inversions

``` r
fst <- read.table("./Population_Genomics/masked.OF_fst.10kb.txt", header=TRUE)
fst <- fst %>% dplyr::mutate(fst50 = zoo::rollmean(WEIGHTED_FST, k = 5,fill = c(0,0,0)), fst100 = zoo::rollmean(WEIGHTED_FST, k = 10,fill = c(0,0,0)), fst500 = zoo::rollmean(WEIGHTED_FST, k = 50,fill = c(0,0,0)), fst1000 = zoo::rollmean(WEIGHTED_FST, k = 100,fill = c(0,0,0)))

SNP <-c(1:(nrow(fst)))
mydf<-data.frame(SNP,fst)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf

mydf$CHR <- as.numeric(mydf$CHR)  

md4.fst <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WEIGHTED_FST",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue =
    ## "WEIGHTED_FST", : NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.fst[[1]]

dfmsplit <- split(dfm, dfm$chrom)

xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})  


top95.dfm <- dfm[which(dfm$WEIGHTED_FST > quantile(dfm$WEIGHTED_FST,0.95) & dfm$fst50 > quantile(dfm$fst50,0.95) & dfm$fst100 > quantile(dfm$fst100,0.95)),]

write(paste(top95.dfm$CHROM, top95.dfm$BIN_START,top95.dfm$BIN_END,top95.dfm$WEIGHTED_FST, sep='\t'),
      "./Population_Genomics/top95.rolling.FST.bed")
```

``` bash
bedtools merge -i ./Population_Genomics/top95.rolling.FST.bed -c 4 -o max -d 10000 | mawk '$3 - $2 > 100000' > ./Population_Genomics/95%.outlier.merged.bed

bedtools intersect -a ./Population_Genomics/top95.rolling.FST.bed -b ./Population_Genomics/95%.outlier.merged.bed > ./Population_Genomics/95%.outlier.10kb.bed


bedtools intersect -v -a ./Population_Genomics/95%.outlier.10kb.bed -b ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed > ./Population_Genomics/95%.outlier.10kb.no_large_inversion.bed


bedtools intersect -v -a ./Population_Genomics/95%.outlier.10kb.no_large_inversion.bed -b ./Population_Genomics/Intermediate_Files/delly.merged.inversions.masked.bed > ./Population_Genomics/95%.outlier.10kb.no_inversion.bed 


bedtools intersect -a ./Population_Genomics/95%.outlier.10kb.no_large_inversion.bed -b ./Population_Genomics/Intermediate_Files/delly.merged.inversions.masked.bed  > ./Population_Genomics/95%.outlier.10kb.delly.inversion.bed



bedtools intersect -a <( grep TRUE ./Population_Genomics/Masked_NoInbred_FST_Outflank.bed) -b ./Population_Genomics/95%.outlier.10kb.no_inversion.bed | grep TRUE > ./Population_Genomics/95%.outlier.SNPs.no_inversion.bed

bedtools intersect -a <( grep TRUE ./Population_Genomics/Masked_NoInbred_FST_Outflank.bed) -b ./Population_Genomics/95%.outlier.10kb.delly.inversion.bed | grep TRUE > ./Population_Genomics/95%.outlier.SNPs.delly.inversion.bed


bedtools intersect -a ./Population_Genomics/masked.OF_fst.10kb.bed -b ./Population_Genomics/95%.outlier.10kb.no_inversion.bed | grep TRUE | cut -f1-4,6 > ./Population_Genomics/95%.10kb.OF.FST.bed

bedtools intersect -a ./Population_Genomics/masked.OF_fst.10kb.bed -b ./Population_Genomics/95%.outlier.10kb.delly.inversion.bed | grep TRUE | cut -f1-4,6 > ./Population_Genomics/95%.10kb.OF.delly.inversion.FST.bed

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\t95")  <(mawk '{print $0 "\tTRUE"}' ./Population_Genomics/95%.10kb.OF.FST.bed) > ./Population_Genomics/95%.outlier.10kb.ni.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\tDI95")  <(mawk '{print $0 "\tTRUE"}' ./Population_Genomics/95%.10kb.OF.delly.inversion.FST.bed) > ./Population_Genomics/95%.outlier.10kb.delly.inversion.txt
```

``` r
fst<-read.table("./Population_Genomics/masked.OF_fst.10kb.txt", header=TRUE)
Ffst <- read.table("./Population_Genomics/95%.outlier.10kb.ni.txt", header=TRUE)
DIfst <- read.table("./Population_Genomics/95%.outlier.10kb.delly.inversion.txt", header=TRUE)
big.fst <- join(fst,Ffst, by =c("CHROM", "BIN_START"),type="full")
big.fst <- join(big.fst,DIfst, by =c("CHROM", "BIN_START"),type="full")

SNP <-c(1:(nrow(big.fst)))
mydf<-data.frame(SNP,big.fst)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf

mydf$CHR <- as.numeric(mydf$CHR)  

md4.fst <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WEIGHTED_FST",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue =
    ## "WEIGHTED_FST", : NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.fst[[1]]


dfm.sub1 <- dfm[which(dfm$X95 == TRUE),]
dfm.sub2 <- dfm[which(dfm$DI95 == TRUE),]

#dfm.sub1 <- top95.dfm
  
dfm.plot <- anti_join(dfm, dfm.sub1, by=c("CHROM", "BIN_START"))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})  

m2.m <- ggplot(dfm.plot, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt), fill = as.factor(CHR)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") +ggtitle("*F<sub>ST</sub>*") +
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal)
  #scale_fill_manual(values=c(pal[1],pal[2]))

m2o.m <-m2.m+geom_point(data=dfm.sub1,alpha=0.8,size=2.5,shape=21, stroke = 0.01, color="black" ) + ggtitle("")
m2oDI.m <-m2.m+geom_point(data=dfm.sub2,alpha=0.8,size=2.5,shape=21, stroke = 0.01, color="black" )
```

``` bash
bcftools view --threads 40 -R ./Population_Genomics/95%.outlier.SNPs.no_inversion.bed ./Population_Genomics/SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz -O z -o ./Population_Genomics/Top.outliers.noinversions.vcf.gz

bcftools index ./Population_Genomics/Top.outliers.noinversions.vcf.gz

bcftools view --threads 40 -R ./Population_Genomics/95%.outlier.SNPs.delly.inversion.bed ./Population_Genomics/SNP.MASKED.noinbred.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.vcf.gz -O z -o ./Population_Genomics/Top.delly.inversion.outliers.vcf.gz

bcftools index ./Population_Genomics/Top.delly.inversion.outliers.vcf.gz

cp ./Population_Genomics/Top.delly.inversion.outliers.vcf.gz ./Population_Genomics/Top.outliers.noinversions.vcf.gz ./Population_Genomics/95%.outlier.SNPs.no_inversion.bed ./Population_Genomics/95%.outlier.SNPs.delly.inversion.bed ./Population_Genomics/Output/Outliers
```

``` r
vcf.out <- read.vcfR( "./Population_Genomics/Top.outliers.noinversions.vcf.gz", verbose = FALSE )
vcf.out.inv <- read.vcfR( "./Population_Genomics/Top.delly.inversion.outliers.vcf.gz", verbose = FALSE )

gid.out <- vcfR2genind(vcf.out, return.alleles = TRUE)
gid.out.inv <- vcfR2genind(vcf.out.inv, return.alleles = TRUE)

pop(gid.out) <- popmap.ni$POP
pop(gid.out.inv) <- popmap.ni$POP
```

``` r
gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.ni$IND
gp.popmap$pop <- popmap.ni$POP

target <- c("CL", "SL", "OBOYS","LOLA","DEBY","HC-VA","CLP","CS","HC","NEH","UMFS","SM","HI")
target <- rev(target)
gp.popmap <-gp.popmap %>% arrange(factor(pop, levels = target))


Chrom_Outlier_GP_Plot <- function(chrom,ncbi, vcf, pre) {

plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("IsleofDogs1"))

outlier.gp <- genotype_plot(vcf=vcf,chr= ncbi,popmap = gp.popmap, start=1, end=94778955,cluster=FALSE,colour_scheme=plot_pal) 


dfx <- min(ggplot_build(outlier.gp$genotypes)$data[[1]]$x) 
dfmax <- max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x)
dflen <- length(ggplot_build(outlier.gp$genotypes)$data[[1]]$x)/78
#dfgt <-length(ggplot_build(outlier.gp$genotypes))

dloc <- (dfmax-dfx)/dflen
print(dloc)
print(dflen)

if(dloc > 50000 & dflen > 50){
  dloc <- dloc/3
}

if(dloc > 50000 & dflen < 50){
  dloc <- dloc/2
}

#ddelta <- dloc*dflen/20

if(dflen > 500){
  #ddelta <- dloc*dflen/20
  ddelta <- dloc*50
  dddelta <- ddelta/10000
  }else if(dflen < 500 & dflen > 200){
  ddelta <- dloc*15
  dddelta <- ddelta/1000
  }else if(dflen < 200 & dflen > 100){
  ddelta <- dloc*10
  dddelta <- ddelta/100
  }else if(dflen < 100 & dflen > 20){
  ddelta <- dloc*5
  dddelta <- ddelta/10
  } else{
  ddelta <- dloc*5
  dddelta <- ddelta/1
  }


  
dfx <- dfx -ddelta/2
#ddelta <- ddelta/2

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)


fig <- outlier.gp$genotypes + geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X) - ddelta/2 , max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x) +dddelta)) + scale_color_manual(values=p_noinbred , name="Population/Line") +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_blank())+
scale_shape_manual(values=c(18,16),name="Origin") +
scale_size_manual(values=c(3.0,2.5),name="Origin", guide="none") +
ggtitle(paste("Chr. ", chrom))+
#guides(fill = guide_legend(override.aes = list(size=6)))+ 
#guides(color=guide_legend(override.aes = list(size=6, shape=c(16,16,18,18,16,16,18,16,16,18,16,16,18))))+
#guides(shape = guide_legend(override.aes = list(size=6)))
guides(fill = "none")+ 
guides(color="none")+
theme(plot.title = element_text(margin=margin(0,0,-36,20))) +
guides(shape = "none") 

assign(paste("gp.plot",chrom,pre,"out", sep="."),fig, envir = globalenv())


}
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 159
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 159
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 159
    ##   row_num: 0
    ## Processed variant: 159
    ## All variants processed
    ## [1] 780.7816
    ## [1] 159
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 871
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 871
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 871
    ##   row_num: 0
    ## Processed variant: 871
    ## All variants processed
    ## [1] 1945.132
    ## [1] 871
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 68
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 68
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 68
    ##   row_num: 0
    ## Processed variant: 68
    ## All variants processed
    ## [1] 1705.327
    ## [1] 68
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 2065
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2065
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 2065
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2065
    ## All variants processed
    ## [1] 6442.976
    ## [1] 2065
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 31
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 31
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 31
    ##   row_num: 0
    ## Processed variant: 31
    ## All variants processed
    ## [1] 64217.98
    ## [1] 31

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 111
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 111
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 111
    ##   row_num: 0
    ## Processed variant: 111
    ## All variants processed
    ## [1] 2531.589
    ## [1] 111
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 877
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 877
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 877
    ##   row_num: 0
    ## Processed variant: 877
    ## All variants processed
    ## [1] 1444.147
    ## [1] 877
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 207
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 207
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 207
    ##   row_num: 0
    ## Processed variant: 207
    ## All variants processed
    ## [1] 532.3866
    ## [1] 207
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 510
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 510
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 510
    ##   row_num: 0
    ## Processed variant: 510
    ## All variants processed
    ## [1] 84953.94
    ## [1] 510
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 295
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 295
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 295
    ##   row_num: 0
    ## Processed variant: 295
    ## All variants processed
    ## [1] 2683.175
    ## [1] 295
    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 104
    ##   header_line: 105
    ##   variant count: 8
    ##   column count: 87
    ## Meta line 104 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 8
    ##   Character matrix gt cols: 87
    ##   skip: 0
    ##   nrows: 8
    ##   row_num: 0
    ## Processed variant: 8
    ## All variants processed
    ## [1] 369.5781
    ## [1] 8

``` r
png(filename=paste(prefix.OF,"Figure.4.Large.FST.Outlier.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

layout <-" 
AAAAAAAAAAAAA
BBBCCCDDEEEEF
BBBCCCDDEEEEF
BBBCCCDDEEEEF
BBBCCCDDEEEEF
BBBCCCDDEEEEF
"

fst.patch <- m2o.m + gp.plot.2.SNP.out  + (gp.plot.3.SNP.out + plot_layout(tag_level = "new")) +  (gp.plot.5.SNP.out+ plot_layout(tag_level = "new")) + (gp.plot.7.SNP.out+ plot_layout(tag_level = "new")) + (gp.plot.8.SNP.out+ plot_layout(tag_level = "new"))  + plot_layout(design=layout, guides="collect") 

#fst.patch & theme(plot.title = element_text(margin=margin(0,0,-10,0)))


 fst.patch  + plot_annotation(tag_levels = 'A') 
dev.off()
```

    ## png 
    ##   2

``` r
filename <- paste(prefix.OFS,"Figure.S10.Large.FST.DI.Outlier.png", sep="")
png(filename=filename, type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

layout <-" 
AAAAAAAAAAAAAAAAAAAAAA
BBCCCCCCCCDDDEEEEFFFFG
BBCCCCCCCCDDDEEEEFFFFG
BBCCCCCCCCDDDEEEEFFFFG
BBCCCCCCCCDDDEEEEFFFFG
BBCCCCCCCCDDDEEEEFFFFG
"



#m2oDI.m / (gp.plot.2.SNPDI.out | gp.plot.3.SNPDI.out |  gp.plot.4.SNPDI.out |  gp.plot.5.SNPDI.out | gp.plot.7.SNPDI.out | gp.plot.8.SNPDI.out ) + plot_layout(design=layout, guides="collect")
m2oDI.m + gp.plot.2.SNPDI.out  + (gp.plot.3.SNPDI.out + plot_layout(tag_level = "new")) +  (gp.plot.4.SNPDI.out+ plot_layout(tag_level = "new")) +  (gp.plot.5.SNPDI.out+ plot_layout(tag_level = "new")) + (gp.plot.7.SNPDI.out+ plot_layout(tag_level = "new")) + (gp.plot.8.SNPDI.out+ plot_layout(tag_level = "new")) + plot_layout(design=layout, guides="collect") + plot_annotation(tag_levels = 'A') 
dev.off()
```

    ## png 
    ##   2

## Structural Variant Divergence

``` bash
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tNAME\tTYPE\tLENGTH") ./Population_Genomics/sv.full.bed > ./Population_Genomics/sv.full.bed.txt

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tNAME\tTYPE\tLENGTH") ./Population_Genomics/sv.full.nomissing.bed > ./Population_Genomics/sv.full.nomissing.bed.txt
```

``` r
sv_data <- read.table("./Population_Genomics/sv.full.bed.txt",header = TRUE)
sv_data_nm <- read.table("./Population_Genomics/sv.full.nomissing.bed.txt",header = TRUE)




group_by(sv_data, TYPE) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(LENGTH, na.rm = TRUE),
    sd = sd(LENGTH, na.rm = TRUE),
    se=std.error(LENGTH, na.rm = TRUE) )%>% mutate_if(is.numeric, round, 0)
```

    ## # A tibble: 5 × 5
    ##   TYPE   count   mean     sd    se
    ##   <fct>  <dbl>  <dbl>  <dbl> <dbl>
    ## 1 BND    33012      1      0     0
    ## 2 DEL   216912    823  10412    22
    ## 3 DUP    15347  14784  57501   464
    ## 4 INS     7310     28      4     0
    ## 5 INV     6808 118478 235738  2857

``` r
group_by(sv_data_nm, TYPE) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(LENGTH, na.rm = TRUE),
    sd = sd(LENGTH, na.rm = TRUE),
    se=std.error(LENGTH, na.rm = TRUE) )%>% mutate_if(is.numeric, round, 0)
```

    ## # A tibble: 5 × 5
    ##   TYPE  count  mean     sd    se
    ##   <fct> <dbl> <dbl>  <dbl> <dbl>
    ## 1 BND    1153     1      0     0
    ## 2 DEL   21158   354   5508    38
    ## 3 DUP    1611  5996  40036   997
    ## 4 INS     413    29      4     0
    ## 5 INV     519 52839 155339  6819

``` r
mask_cn <- read.table("./Population_Genomics/masked.cnv.tab",header = FALSE)

mask_cn <- mask_cn[which(m.freq.ni > 0.05),]

getVst <- function(dat, group1, group2) {
  groupLevels <- levels(droplevels(group2))
  j <- 0
  total_dat <-0
  for (i in levels(group2)){
  j <- j+1
  assign(paste("dat",j, sep=""), na.omit(dat[group1==groupLevels[groupLevels==i]]))
  total_dat <- cbind(total_dat,eval(as.name(paste("dat",j,sep=""))))
  
  }
  
  total_dat <- total_dat[,-1]

  Vtotal <- rowVars(as.matrix(total_dat))
  
  j <- 0
  Vgroup <- 0
  for (i in levels(droplevels(group2))){
  j <- j+1
  dat <- eval(as.name(paste("dat",j,sep="")))
  
  Vgroup <- Vgroup + rowVars(as.matrix(dat))*length(dat[1,])
  
  }
  Vgroup <- Vgroup / length(total_dat[1,])
 
  
  Vst <- (Vtotal-Vgroup) / Vtotal
  

  
  return(Vst)
}  
mask_cn$VST <- getVst(mask_cn[,-c(1,2,3)],popmap$POP,popmap.ni$POP)

mask_cn$VSTW <- getVst(mask_cn[,-c(1,2,3,94)],popmap$POP,popmap.wild$POP)

mask_cn$VSTWAE <- getVst(mask_cn[,-c(1,2,3,94,95)],popmap$POP,popmap.wildAE$POP)
mask_cn$VSTWvD <- getVst(mask_cn[,-c(1,2,3,94,95,96)],popmap$Type,popmap.ni$Type)
```

``` r
mask_cnVST <- mask_cn[complete.cases(mask_cn$VST),]

mask_cnWVST <- mask_cn[complete.cases(mask_cn$VSTW),]

mask_cnWAEVST <- mask_cn[complete.cases(mask_cn$VSTWAE),]

mask_cnWvDVST <- mask_cn[complete.cases(mask_cn$VSTWvD),]

write(paste(mask_cnVST$V1,mask_cnVST$V2,mask_cnVST$VST, sep='\t'),"./Population_Genomics/Intermediate_Files/masked_vst.tab")


write(paste(mask_cnWVST$V1,mask_cnWVST$V2,mask_cnWVST$VSTW, sep='\t'),"./Population_Genomics/Intermediate_Files/masked_wild_vst.tab")


write(paste(mask_cnWAEVST$V1,mask_cnWAEVST$V2,mask_cnWAEVST$VSTWAE, sep='\t'),"./Population_Genomics/Intermediate_Files/masked_wildAE_vst.tab")

summary(mask_cnVST$VST)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -0.16434 -0.01449  0.04545  0.07086  0.12000  0.91913

``` r
std.error(mask_cnVST$VST)
```

    ## [1] 0.001086641

``` bash
mawk '{print $1 "\t" $2-1 "\t" $2 "\t" $3}' ./Population_Genomics/Intermediate_Files/masked_vst.tab > ./Population_Genomics/Intermediate_Files/masked_vst.bed

bedtools intersect -b <(sort -k 1,1 -k2,2n ./Population_Genomics/Intermediate_Files/masked_vst.bed) -a <(sort -k 1,1 -k2,2n ./other_files/10kb.bed) -wa -wb -sorted |bedtools merge -d -1 -c 7 -o mean > ./Population_Genomics/Intermediate_Files/mask.vst.10kb.bed

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_VST") ./Population_Genomics/Intermediate_Files/mask.vst.10kb.bed > ./Population_Genomics/Intermediate_Files/masked.all.windowed.weir.vst

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_VST") ./Population_Genomics/Intermediate_Files/masked_vst.bed > ./Population_Genomics/Intermediate_Files/masked_vst.bed.txt
```

``` r
#mvst <- read.table("./Population_Genomics/Intermediate_Files/masked.all.windowed.weir.vst", header = TRUE)
#mvst <- read.table("masked.wild.AE.windowed.weir.vst", header = TRUE)

mvst <- read.table("./Population_Genomics/Intermediate_Files/masked_vst.tab", header = FALSE)
colnames(mvst) <- c("CHROM","BIN_START", "WEIGHTED_VST")

mvst <- mvst %>% dplyr::rename(MASK_VST= WEIGHTED_VST)
```

``` r
SNP<-c(1:(nrow(mvst)))
mydf<-data.frame(SNP,mvst)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf
mydf$CHR <- as.numeric(mydf$CHR)  

m2.cnv <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="MASK_VST",snp="SNP",logTransform=FALSE,sigLine=NA, ymax=1, relative.positions = TRUE)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- m2.cnv[[1]]
dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

pal = wes_palette("Zissou1", 10, type = "continuous")

m2.cnv.total <- ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5, aes(alpha=0.33))+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(min(dfm$marker),1))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") +
  #ggtitle("Masked") +
  theme_classic()+
  labs(y="*V<sub>ST</sub>*") +theme(axis.title.y = element_markdown())+ 
  #ylab(*bquote(V[ST])*)+
  #theme(axis.title.y = element_text(face = "italic"))+
  scale_color_manual(values=c(pal[1], pal[3])) 
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
#ylab(expression(paste("C", H[4], " (", mu,"mol ", L^-1,")"))) +

mydf.sig <- mydf[mydf$MASK_VST>quantile(mydf$MASK_VST, 0.995,na.rm=TRUE) ,]
mydf.sig$snp <- mydf.sig$SNP
dfm.sub <- merge(dfm,mydf.sig, by = "snp")

nrow(mydf.sig)
```

    ## [1] 70

``` r
m2.cnv <-m2.cnv.total+geom_point(data=dfm.sub,shape = 17,alpha=1,size=1.5)
```

``` r
gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.ni$IND
gp.popmap$pop <- popmap.ni$POP

target <- c("CL", "SL", "OBOYS","LOLA","DEBY","HC-VA","CLP","CS","HC","NEH","UMFS","SM","HI")
target <- rev(target)
gp.popmap <-gp.popmap %>% arrange(factor(pop, levels = target))

max_cn <-mask_cn[which(mask_cn$VST > quantile(mask_cn$VST,c(0.995),na.rm=TRUE)),]
max_cn %>% 
  mutate(V1 = str_replace(V1, "NC_035780.1", "1")) %>% 
  mutate(V1 = str_replace(V1, "NC_035781.1", "2")) %>% 
  mutate(V1 = str_replace(V1, "NC_035782.1", "3")) %>% 
  mutate(V1 = str_replace(V1, "NC_035783.1", "4")) %>% 
  mutate(V1 = str_replace(V1, "NC_035784.1", "5")) %>% 
  mutate(V1 = str_replace(V1, "NC_035785.1", "6")) %>% 
  mutate(V1 = str_replace(V1, "NC_035786.1", "7")) %>% 
  mutate(V1 = str_replace(V1, "NC_035787.1", "8")) %>%
  mutate(V1 = str_replace(V1, "NC_035788.1", "9")) %>% 
  mutate(V1 = str_replace(V1, "NC_035789.1", "10")) -> max_cn
max_cnvs <- paste(max_cn$V1,max_cn$V2,sep="_")
max_cnv_df <- t(max_cn[,4:93])
rownames(max_cnv_df) <- popmap$IND
colnames(max_cnv_df) <- max_cnvs
max_cnv_df <- max_cnv_df[c(1:36,40:45,51:56,57:62,67:90),]
max_cnv_df <- as.data.frame(max_cnv_df)


max_cnv_df <-max_cnv_df %>% dplyr::arrange(factor(row.names(max_cnv_df), levels = gp.popmap$ind))

melted_mask_cnv_df <- reshape2::melt(as.matrix(max_cnv_df))
colnames(melted_mask_cnv_df)<- c("IND","LOC","CN")
cnv.loc <-str_split(melted_mask_cnv_df$LOC,"_", simplify = TRUE)
melted_mask_cnv_df$CHR <- cnv.loc[,1]
melted_mask_cnv_df$BP <- cnv.loc[,2]
melted_mask_cnv_df$CHR <- factor(melted_mask_cnv_df$CHR, levels= c(1,2,3,4,5,6,7,8,9,10))

cnv.plot <- ggplot(melted_mask_cnv_df) + geom_tile(aes(x=LOC,y=IND,fill=CN)) + 
  scale_y_discrete(limits = rev(levels(as.factor(melted_mask_cnv_df$IND)))) + 
  scale_fill_gradientn(colours=pal)+
   facet_grid(~ CHR, switch = "x", scales = "free_x", space = "free_x") +
  scale_x_discrete(expand= c(0,0))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) 
  #guides(fill= guide_legend(order=1))


dfx <- min(ggplot_build(cnv.plot)$data[[1]]$x) 
dfmax <- max(ggplot_build(cnv.plot)$data[[1]]$x)

ddelta <- (dfmax-dfx)/50
dfx <- dfx -ddelta

new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- gp.popmap$ind
new_df <- join(new_df, popmap.ni)
```

    ## Joining by: IND

``` r
new_df$Y <- 78:1
new_df$X <- rep(c(dfx-0,dfx-ddelta/2,dfx-ddelta,dfx-ddelta,dfx-ddelta/2,dfx-0),13)

lab.plot <- ggplot()+ geom_point(data=new_df,aes(x=X,y=Y,shape=Type, color=POP, size=Type)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(new_df$X -0.05), max(new_df$X +.05)) )+ 
  scale_y_continuous(expand = c(0, 0), limits = c(min(new_df$Y -0.5), max(new_df$Y +0.5)) )+ 
  scale_color_manual(values=p_noinbred , name="Population/Line", guide="none") +
  theme_void() +
theme(legend.text = element_text(size=12),legend.position = "right", legend.title=element_text(size=14),axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y=element_blank())+
scale_shape_manual(values=c(18,16),name="Origin", guide="none") +
scale_size_manual(values=c(3.5,3),name="Origin", guide="none") 
```

``` r
max_cn %>% 
  mutate(V1 = str_replace(V1,  "^10$", "NC_035789.1")) %>% 
  mutate(V1 = str_replace(V1,  "^2$", "NC_035781.1")) %>% 
  mutate(V1 = str_replace(V1,  "^3$", "NC_035782.1")) %>% 
  mutate(V1 = str_replace(V1,  "^4$", "NC_035783.1")) %>% 
  mutate(V1 = str_replace(V1,  "^5$", "NC_035784.1")) %>% 
  mutate(V1 = str_replace(V1,  "^6$", "NC_035785.1")) %>% 
  mutate(V1 = str_replace(V1,  "^7$", "NC_035786.1")) %>% 
  mutate(V1 = str_replace(V1,  "^8$", "NC_035787.1")) %>%
  mutate(V1 = str_replace(V1,  "^9$", "NC_035788.1")) %>% 
  mutate(V1 = str_replace(V1,  "^1$", "NC_035780.1"))  -> max_cnvw
write(paste(max_cnvw$V1,max_cnvw$V2,max_cnvw$V3, sep= "\t"),"./Population_Genomics/Intermediate_Files/masked_full_outlier_vst.loci")

summary(max_cn$VST)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.6169  0.6438  0.7120  0.7209  0.7813  0.9191

``` r
std.error(max_cn$VST)
```

    ## [1] 0.009050884

``` bash
grep -f <(cut -f 3 ./Population_Genomics/Intermediate_Files/masked_full_outlier_vst.loci) ./Population_Genomics/delly.cnv.masked.bed > ./Population_Genomics/Output/Outliers/outlier.cnv.bed 
```

#### CNV PCA

``` r
pca2 <- dudi.pca(mask_gen.niw,scannf=FALSE,scale=FALSE,nf=4)
#pca3 <- dudi.pca(mask_gen.wild,scannf=FALSE,scale=FALSE,nf=4)
pca4 <- dudi.pca(mask_gen.wildAE,scannf=FALSE,scale=FALSE,nf=4)

out.cnvs <- apply(max_cn[,1:3],1,paste,collapse="_")
out.cnvs <- str_replace_all(out.cnvs, fixed(" "), "")
mask_gen.niw.outlier <- mask_gen.niw[loc=out.cnvs]

pca2_o <- dudi.pca(mask_gen.niw.outlier,scannf=FALSE,scale=FALSE,nf=4)
```

``` r
dffcnv <- as.data.frame(pca2$li[,1:4])
colnames(dffcnv) <- c("PC1","PC2","PC3","PC4")
dffcnv$POP <- popmap.ni$POP
dffcnv$Type <- popmap.ni$Type

dffcnv.o <- as.data.frame(pca2_o$li[,1:4])
colnames(dffcnv.o) <- c("PC1","PC2","PC3","PC4")
dffcnv.o$POP <- popmap.ni$POP
dffcnv.o$Type <- popmap.ni$Type
```

``` r
png(filename=paste(prefix.OE,"PCA.CNV.allpops.1.2.png", sep=""), type="cairo",units="px", width=3500, height=3500, res=200, bg="transparent")
pca_cnv <- ggplot(data = dffcnv,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 6)+
scale_shape_manual(values=c(23,21),name="Origin") +
scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") + 
  #ggtitle("CNV PCA of all non-inbred Populaitons/Lines") + 
guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha=FALSE,size=FALSE , shape= guide_legend(override.aes = list(size =6)))+ 
geom_dl(alpha=1,method = list(cex = 2, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=20))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
pca_cnv

dev.off()
```

    ## png 
    ##   2

``` r
pca_cnv <- ggplot(data = dffcnv,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 5)+
scale_shape_manual(values=c(23,21),name="Origin") +
scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") + 
  #ggtitle("CNV PCA of all non-inbred Populaitons/Lines") + 
guides(fill = guide_legend(override.aes = list(size =4, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha="none",size="none" , shape= guide_legend(override.aes = list(size =4), order=2))+ 
geom_dl(alpha=1,method = list(cex = 1, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=12))

pca_cnv_o <- ggplot(data = dffcnv.o,aes( x=PC1, y=PC2,color=POP,fill=POP,shape=Type, size=3, alpha=0.95)) +
geom_point(alpha=0.95, size = 5)+
scale_shape_manual(values=c(23,21),name="Origin", ) +
scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") + 
  #ggtitle("CNV PCA of all non-inbred Populaitons/Lines") + 
guides(fill = "none", alpha="none",size="none" , shape= "none", color="none")+ 
geom_dl(alpha=1,method = list(cex = 1, font=2,"smart.grid"), aes(label=POP))+
theme_bw() + theme(text = element_text(size=12))






pcag7 <- ggplot(pca4$li, aes(x=Axis1,y=Axis2,color=popmap.wildAE$POP,fill=popmap.wildAE$POP,shape=popmap.wildAE$Type, size=3, alpha=0.95)) +   geom_point()+
scale_fill_manual(values=p_wildAE, name="Population/Line") + scale_color_manual(values=p_wildAE, name="Population/Line") +guides(size=FALSE, alpha=FALSE)+ggtitle("PCA of Atlantic wild Populaitons")+ labs(x="PC1", y="PC2")+theme_bw() + theme(text = element_text(size=20))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
pcag8 <- ggplot(pca4$li, aes(x=Axis3,y=Axis4,color=popmap.wildAE$POP,fill=popmap.wildAE$POP,shape=popmap.wildAE$Type, size=3, alpha=0.95)) +  geom_point()+
scale_fill_manual(values=p_wildAE, name="Population/Line") + scale_color_manual(values=p_wildAE, name="Population/Line") +guides(size=FALSE, alpha=FALSE)+ggtitle("PCA of Atlantic wild Populaitons")+
labs(x="PC3", y="PC4")+theme_bw() + theme(text = element_text(size=20))
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
png(filename=paste(prefix.OE,"CNV.PCA.wildAE.1.2.png", sep=""), type="cairo",units="px", width=2500, height=2500, res=200, bg="transparent")
pcag7
dev.off()
```

    ## png 
    ##   2

``` r
png(filename=paste(prefix.OE,"CNV.PCA.wildAE.3.4.png", sep=""), type="cairo",units="px", width=2500, height=2500, res=200, bg="transparent")
pcag8
dev.off()
```

    ## png 
    ##   2

``` r
layout <- "ABBBBBBBBBBBBBBBBBBBBBB"
s_p <- lab.plot + cnv.plot +plot_layout(design=layout, guides = "collect")


png(filename=paste(prefix.OFS,"Figure.S11.CNV_Divergence.png", sep=""), type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")
layout <- "AAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDFFF
           AAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDFFF
           BCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEFFF
           BCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEFFF
           BCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEFFF
           BCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGFFF
           BCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGFFF
           BCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGFFF"

layout <- "AAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCEEE
           AAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCEEE
           BBBBBBBBBBBBBBBBBBBBBBBDDDDDDDDDEEE
           BBBBBBBBBBBBBBBBBBBBBBBDDDDDDDDDEEE
           BBBBBBBBBBBBBBBBBBBBBBBDDDDDDDDDEEE
           BBBBBBBBBBBBBBBBBBBBBBBFFFFFFFFFEEE
           BBBBBBBBBBBBBBBBBBBBBBBFFFFFFFFFEEE
           BBBBBBBBBBBBBBBBBBBBBBBFFFFFFFFFEEE"


p3 <- m2.cnv + wrap_elements(full = s_p) +  finished_mapnolm_small + pca_cnv +guide_area() + pca_cnv_o +plot_layout(design=layout, guides = "collect") 
p3 + plot_annotation(tag_levels = 'A')
dev.off()
```

    ## png 
    ##   2

## Big Outlier Signal in Wild Atlantic Samples

#### Wild

``` bash
bedtools intersect -b <(sort -k 1,1 -k2,2n ./Population_Genomics/Masked_WildAE_FST_Outflank.bed) -a <(sort -k 1,1 -k2,2n ./other_files/10kb.bed) -wa -wb -sorted  |bedtools merge -d -1 -c 7,8 -o mean,distinct | mawk '{ if ($0 ~ /TRUE/) {print $1 "\t" $2 "\t" $3 "\t" $4 "\tTRUE"} else {print $0}}' > ./Population_Genomics/masked.wildae.OF_fst.10kb.bed

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\tOUTLIER") ./Population_Genomics/masked.wildae.OF_fst.10kb.bed > ./Population_Genomics/masked.wildae.OF_fst.10kb.bed.txt

mkdir -p ./Population_Genomics/Output/Outliers/Wild_Atlantic
cp ./Population_Genomics/masked.wildae.OF_fst.10kb.bed ./Population_Genomics/masked.wildae.OF_fst.10kb.bed.txt ./Population_Genomics/Output/Outliers/Wild_Atlantic
```

``` r
fst <- read.table("./Population_Genomics/masked.wildae.OF_fst.10kb.bed.txt", header=TRUE)


fst <- fst %>% dplyr::mutate(fst50 = zoo::rollmean(WEIGHTED_FST, k = 5,fill = c(0,0,0)), fst100 = zoo::rollmean(WEIGHTED_FST, k = 10,fill = c(0,0,0)))

SNP <-c(1:(nrow(fst)))
mydf<-data.frame(SNP,fst)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf

mydf$CHR <- as.numeric(mydf$CHR)  

md4.fst <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WEIGHTED_FST",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue =
    ## "WEIGHTED_FST", : NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.fst[[1]]

top95.dfm <- dfm[which(dfm$WEIGHTED_FST > quantile(dfm$WEIGHTED_FST,0.95)) ,]

summary(top95.dfm)
```

    ##       SNP                CHROM       BIN_START            BIN_END         
    ##  Min.   :   40   NC_035784.1:662   Min.   :   200000   Min.   :   210000  
    ##  1st Qu.:14341   NC_035785.1:299   1st Qu.: 29907500   1st Qu.: 29917500  
    ##  Median :23008   NC_035783.1:213   Median : 42900000   Median : 42910000  
    ##  Mean   :19909   NC_035780.1:199   Mean   : 43933059   Mean   : 43943059  
    ##  3rd Qu.:26310   NC_035782.1:161   3rd Qu.: 64565000   3rd Qu.: 64575000  
    ##  Max.   :36798   NC_035788.1: 99   Max.   :102790000   Max.   :102800000  
    ##                  (Other)    :211                                          
    ##   WEIGHTED_FST      OUTLIER            fst50                fst100         
    ##  Min.   :0.09891   Mode :logical   Min.   :-0.0008561   Min.   :-0.002521  
    ##  1st Qu.:0.11111   FALSE:1657      1st Qu.: 0.0688297   1st Qu.: 0.054622  
    ##  Median :0.13102   TRUE :187       Median : 0.0986512   Median : 0.085135  
    ##  Mean   :0.14643                   Mean   : 0.1019942   Mean   : 0.087299  
    ##  3rd Qu.:0.16092                   3rd Qu.: 0.1291862   3rd Qu.: 0.115184  
    ##  Max.   :0.73485                   Max.   : 0.5213716   Max.   : 0.362241  
    ##                                                                            
    ##       CHR             chrom           bp                pvalue       
    ##  Min.   : 1.000   5      :662   Min.   :   200000   Min.   :0.09891  
    ##  1st Qu.: 4.000   6      :299   1st Qu.: 29907500   1st Qu.:0.11111  
    ##  Median : 5.000   4      :213   Median : 42900000   Median :0.13102  
    ##  Mean   : 4.773   1      :199   Mean   : 43933059   Mean   :0.14643  
    ##  3rd Qu.: 6.000   3      :161   3rd Qu.: 64565000   3rd Qu.:0.16092  
    ##  Max.   :10.000   9      : 99   Max.   :102790000   Max.   :0.73485  
    ##                   (Other):211                                        
    ##       snp            index           chrom_alt             marker       
    ##  Min.   :   40   Min.   :   29.56   Length:1844        Min.   :0.09891  
    ##  1st Qu.:14341   1st Qu.:14351.44   Class :character   1st Qu.:0.11111  
    ##  Median :23008   Median :22667.67   Mode  :character   Median :0.13102  
    ##  Mean   :19909   Mean   :19757.97                      Mean   :0.14643  
    ##  3rd Qu.:26310   3rd Qu.:26322.91                      3rd Qu.:0.16092  
    ##  Max.   :36798   Max.   :36646.87                      Max.   :0.73485  
    ## 

``` r
write(paste(top95.dfm$CHROM, top95.dfm$BIN_START,top95.dfm$BIN_END,top95.dfm$WEIGHTED_FST, sep='\t'),
      "./Population_Genomics/top95.wildAE.FST.bed")
```

``` bash
bedtools intersect -a ./Population_Genomics/top95.wildAE.FST.bed -b ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed  > ./Population_Genomics/95%.outlier.wildae.10kb.large_inversions.bed

bedtools subtract -a ./Population_Genomics/top95.wildAE.FST.bed -b ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed | bedtools subtract -a - -b ./Population_Genomics/sv.full.bed > ./Population_Genomics/95%.outlier.wildae.10kb.no_inversions.bed

bedtools intersect -a <( grep TRUE ./Population_Genomics/Masked_WildAE_FST_Outflank.bed) -b ./Population_Genomics/95%.outlier.wildae.10kb.large_inversions.bed | grep TRUE > ./Population_Genomics/95%.outlier.SNPs.wildae.large_inversions.bed


bedtools intersect -a <( grep TRUE ./Population_Genomics/Masked_WildAE_FST_Outflank.bed) -b ./Population_Genomics/95%.outlier.wildae.10kb.no_inversions.bed | grep TRUE > ./Population_Genomics/95%.outlier.SNPs.wildae.no_inversions.bed

bedtools subtract -a ./Population_Genomics/top95.wildAE.FST.bed -b ./Population_Genomics/Intermediate_Files/Detected_large_inversions.bed | bedtools intersect -a - -b ./Population_Genomics/delly.merged.inversions.masked.bed > ./Population_Genomics/95%.outlier.SNPs.wildae.10kb.small_inversions.bed

bedtools intersect -a <( grep TRUE ./Population_Genomics/Masked_WildAE_FST_Outflank.bed) -b ./Population_Genomics/95%.outlier.SNPs.wildae.10kb.small_inversions.bed | grep TRUE > ./Population_Genomics/95%.outlier.SNPs.wildae.small_inversions.bed


bcftools view --threads 40 -R ./Population_Genomics/95%.outlier.SNPs.wildae.large_inversions.bed ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz -O z -o ./Population_Genomics/Top.wildae.outliers.large_inversions.vcf.gz

bcftools view --threads 40 -R ./Population_Genomics/95%.outlier.SNPs.wildae.no_inversions.bed ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz -O z -o ./Population_Genomics/Top.wildae.outliers.no_inversions.vcf.gz

bcftools view --threads 40 -R ./Population_Genomics/95%.outlier.SNPs.wildae.small_inversions.bed ./Population_Genomics/SNP.MASKED.wildAE.g1.maf05.max2alleles.FIL.vcf.gz -O z -o ./Population_Genomics/Top.wildae.outliers.small_inversions.vcf.gz

bcftools index ./Population_Genomics/Top.wildae.outliers.small_inversions.vcf.gz
bcftools index ./Population_Genomics/Top.wildae.outliers.no_inversions.vcf.gz
bcftools index ./Population_Genomics/Top.wildae.outliers.large_inversions.vcf.gz

bedtools intersect -a ./Population_Genomics/masked.wildae.OF_fst.10kb.bed -b ./Population_Genomics/95%.outlier.SNPs.wildae.large_inversions.bed -wa | grep TRUE |bedtools merge -i - -d -1 -c 4,5 -o first| cut -f1-4,6 > ./Population_Genomics/95%.10kb.WILDAE.OF.FST.large_inversions.bed

bedtools intersect -a ./Population_Genomics/masked.wildae.OF_fst.10kb.bed -b ./Population_Genomics/95%.outlier.SNPs.wildae.no_inversions.bed -wa | grep TRUE |bedtools merge -i - -d -1 -c 4,5 -o first| cut -f1-4,6 > ./Population_Genomics/95%.10kb.WILDAE.OF.FST.no_inversions.bed

bedtools intersect -a ./Population_Genomics/masked.wildae.OF_fst.10kb.bed -b ./Population_Genomics/95%.outlier.SNPs.wildae.small_inversions.bed -wa | grep TRUE |bedtools merge -i - -d -1 -c 4,5 -o first| cut -f1-4,6 > ./Population_Genomics/95%.10kb.WILDAE.OF.FST.small_inversions.bed

cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\t95")  <(mawk '{print $0 "\tTRUE"}' ./Population_Genomics/95%.10kb.WILDAE.OF.FST.small_inversions.bed) > ./Population_Genomics/95%.outlier.wildae.10kb.small_inversions.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\t95NI")  <(mawk '{print $0 "\tTRUE"}' ./Population_Genomics/95%.10kb.WILDAE.OF.FST.no_inversions.bed) > ./Population_Genomics/95%.outlier.wildae.10kb.no_inversions.txt
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tWEIGHTED_FST\t95NLI")  <(mawk '{print $0 "\tTRUE"}' ./Population_Genomics/95%.10kb.WILDAE.OF.FST.large_inversions.bed) > ./Population_Genomics/95%.outlier.wildae.10kb.large_inversions.txt


cp ./Population_Genomics/Top.wildae.outliers.small_inversions.vcf.gz ./Population_Genomics/Top.wildae.outliers.no_inversions.vcf.gz ./Population_Genomics/Top.wildae.outliers.large_inversions.vcf.gz  ./Population_Genomics/95%.outlier.SNPs.wildae.large_inversions.bed ./Population_Genomics/95%.outlier.SNPs.wildae.no_inversions.bed ./Population_Genomics/95%.outlier.SNPs.wildae.small_inversions.bed ./Population_Genomics/Output/Outliers/Wild_Atlantic
```

``` r
fst<-read.table("./Population_Genomics/masked.wildae.OF_fst.10kb.bed.txt", header=TRUE)
Ffst <- read.table("./Population_Genomics/95%.outlier.wildae.10kb.small_inversions.txt", header=TRUE)
FfstNI <- read.table("./Population_Genomics/95%.outlier.wildae.10kb.no_inversions.txt", header=TRUE)
FfstNLI <- read.table("./Population_Genomics/95%.outlier.wildae.10kb.large_inversions.txt", header=TRUE)

big.fst <- join(fst,Ffst, by =c("CHROM", "BIN_START"),type="full")
big.fst <- join(big.fst,FfstNI, by =c("CHROM", "BIN_START"),type="full")
big.fst <- join(big.fst,FfstNLI, by =c("CHROM", "BIN_START"),type="full")

SNP <-c(1:(nrow(big.fst)))
mydf<-data.frame(SNP,big.fst)
mydf$CHR <- mydf$CHROM
mydf %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> mydf

mydf$CHR <- as.numeric(mydf$CHR)  

md4.fst <-ggman(mydf,chrom="CHR",bp="BIN_START",pvalue="WEIGHTED_FST",snp="SNP",logTransform=FALSE,sigLine=NA, relative.positions = TRUE)
```

    ## Warning in ggman(mydf, chrom = "CHR", bp = "BIN_START", pvalue =
    ## "WEIGHTED_FST", : NaNs produced

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

``` r
dfm <- md4.fst[[1]]


dfm.sub2 <- dfm[which(dfm$X95 == TRUE),]
dfm.sub1 <- dfm[which(dfm$X95NLI == TRUE),]
dfm.sub3 <- dfm[which(dfm$X95NI == TRUE),]

dfm.plot1 <- anti_join(dfm, dfm.sub1, by=c("CHROM", "BIN_START"))
dfm.plot2 <- anti_join(dfm, dfm.sub2, by=c("CHROM", "BIN_START"))
dfm.plot3 <- anti_join(dfm, dfm.sub3, by=c("CHROM", "BIN_START"))

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})  

m2.m1 <- ggplot(dfm.plot1, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.13,1))+
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") +ggtitle("*F<sub>ST</sub>*") +
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  #scale_fill_manual(values=col_pal[c(1,2,3,4,5,6,7,9)])
  scale_fill_manual(values=col_pal[c(5,6)])
  #scale_fill_manual(values=c(pal[1],pal[2]))

m2o.m1 <-m2.m1+geom_point(data=dfm.sub1,alpha=0.8,size=3.0,shape=25, stroke = 0.01, color="black", aes( fill = as.factor(CHR)))

m2.m2 <- ggplot(dfm.plot2, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5,alpha=0.25)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.13,1))+
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal[c(1,2,3,4,5,9)])

m2o.m2 <-m2.m2+geom_point(data=dfm.sub2,alpha=0.8,size=3.0,shape=24, stroke = 0.01, color="black", aes( fill = as.factor(CHR)))


m2.m3 <- ggplot(dfm.plot3, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5,alpha=0.2)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.13,0.8))+
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal[c(1,2,3,4,5,7,9)])


m2o.m3 <-m2.m3+geom_point(data=dfm.sub3,alpha=0.8,size=3.0,shape=23, stroke = 0.01, color="black", aes( fill = as.factor(CHR)))

m2.m3 <- ggplot(dfm.plot3, aes(x= index, y=WEIGHTED_FST, colour = as.factor(chrom_alt)))+
  geom_point(size=0.5,alpha=0.2)+  
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0)) +
  guides(colour = "none",alpha="none",size="none", fill="none") +
  scale_y_continuous(expand = c(0,0),limits=c(-0.13,0.8))+
  
  #scale_alpha_continuous(range=c(0.1,0.95))+
  #labs(x = "Chromosome") +ggtitle("*F<sub>ST</sub>*") +
  theme_classic() +
  labs(y="*F<sub>ST</sub>*") +theme(axis.title.y = element_markdown(),plot.title = element_markdown(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("#5e5e5e","#242424"))+
  scale_fill_manual(values=col_pal[c(1,2,3,4,5,6,7,9)])
  #scale_fill_manual(values=c(pal[1],pal[2]))

m2o.mALL <- m2.m3+geom_point(data=dfm.sub3,alpha=0.8,size=3.0,shape=23, stroke = 0.1, color="black", aes( fill = as.factor(CHR))) + geom_point(data=dfm.sub1,alpha=0.8,size=3.0,shape=25, stroke = 0.1, color="black", aes( fill = as.factor(CHR))) + geom_point(data=dfm.sub2,alpha=0.8,size=3.0,shape=24, stroke = 0.1, color="black", aes( fill = as.factor(CHR)))
```

``` r
vcf.out1 <- read.vcfR( "./Population_Genomics/Top.wildae.outliers.large_inversions.vcf.gz", verbose = FALSE )
vcf.out3 <- read.vcfR( "./Population_Genomics/Top.wildae.outliers.no_inversions.vcf.gz", verbose = FALSE )
vcf.out2 <- read.vcfR( "./Population_Genomics/Top.wildae.outliers.small_inversions.vcf.gz", verbose = FALSE )

gid.out1 <- vcfR2genind(vcf.out1, return.alleles = TRUE) 
gid.out2 <- vcfR2genind(vcf.out2, return.alleles = TRUE) 
gid.out3 <- vcfR2genind(vcf.out3, return.alleles = TRUE) 

pop(gid.out1) <- popmap.wildAE$POP
pop(gid.out2) <- popmap.wildAE$POP
pop(gid.out3) <- popmap.wildAE$POP
```

``` r
shape_n <- c( "Within Large Inversions" , "Within Small Inversions" , "Outside Structural Variants"  )
shapes <- c( "Within Large Inversions" = 25, "Within Small Inversions" =24, "Outside Structural Variants" =23 )

pca.out1 <- dudi.pca(gid.out1 ,scannf=FALSE,scale=FALSE,nf=2)
pca.out1.p.df <- pca.out1$li
pca.out1.p.df$POP <- popmap.wildAE$POP
pca.out1.p.df$TYPE <- factor(shape_n)[1]

PCA1PC1 <- round(pca.out1$eig[1]/sum(pca.out1$eig)*100,1)
PCA1PC2 <- round(pca.out1$eig[2]/sum(pca.out1$eig)*100,1)


pca.out2 <- dudi.pca(gid.out2 ,scannf=FALSE,scale=FALSE,nf=2)
pca.out2.p.df <- pca.out2$li
pca.out2.p.df$POP <- popmap.wildAE$POP
pca.out2.p.df$TYPE <- factor(shape_n)[2]

PCA2PC1 <- round(pca.out2$eig[1]/sum(pca.out2$eig)*100,1)
PCA2PC2 <- round(pca.out2$eig[2]/sum(pca.out2$eig)*100,1)


pca.out3 <- dudi.pca(gid.out3 ,scannf=FALSE,scale=FALSE,nf=2)
pca.out3.p.df <- pca.out3$li
pca.out3.p.df$POP <- popmap.wildAE$POP
pca.out3.p.df$TYPE <- factor(shape_n)[3]

PCA3PC1 <- round(pca.out3$eig[1]/sum(pca.out3$eig)*100,1)
PCA3PC2 <- round(pca.out3$eig[2]/sum(pca.out3$eig)*100,1)

pca.out.plot1 <- ggplot(pca.out1.p.df, aes(x=Axis1,y=Axis2,color=POP,fill=POP, shape=TYPE,size=3, alpha=0.9)) + geom_point(alpha=0.9, size = 4, color="black")+
scale_shape_manual(values=shapes,name="Genomic Architecture\nof Outlier SNPs") +
scale_fill_manual(values=p_wildAE, name="Population") + scale_color_manual(values=p_wildAE , name="Population") +
xlab(paste("PC1 (",PCA1PC1,"%)", sep = ""))+ ylab(paste("PC2 (",PCA1PC2,"%)", sep = "")) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
guides(fill = "none", alpha="none",size="none", shape=guide_legend(order=2) ,color="none")+ 
#xlim(-110,65)+ ylim(-40,90) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP, color=POP))+
theme_bw() + theme(text = element_text(size=12)) + theme(legend.position = "right", legend.title.align=0.5, legend.text.align = 0.5 )

pca.out.plot2 <- ggplot(pca.out2.p.df, aes(x=Axis1,y=Axis2,color=POP,fill=POP, shape=TYPE,size=3, alpha=0.9)) + geom_point(alpha=0.9, size = 4, color="black")+
scale_shape_manual(values=shapes,name="Genomic Architecture\nof Outlier SNPs") +
scale_fill_manual(values=p_wildAE, name="Population") + scale_color_manual(values=p_wildAE , name="Population") +
xlab(paste("PC1 (",PCA2PC1,"%)", sep = ""))+ ylab(paste("PC2 (",PCA2PC2,"%)", sep = "")) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
guides(fill = "none", alpha="none",size="none", shape=guide_legend(order=2) ,color="none")+ 
#xlim(-110,65)+ ylim(-40,90) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP, color=POP))+
theme_bw() + theme(text = element_text(size=12)) + theme(legend.position = "right", legend.title.align=0.5, legend.text.align = 0.5 )

pca.out.plot3 <- ggplot(pca.out3.p.df, aes(x=Axis1,y=Axis2,color=POP,fill=POP, shape=TYPE,size=3, alpha=0.9)) + geom_point(alpha=0.9, size = 4, color="black")+
scale_shape_manual(values=shapes,name="Genomic Architecture\nof Outlier SNPs") +
scale_fill_manual(values=p_wildAE, name="Population") + scale_color_manual(values=p_wildAE , name="Population") +
xlab(paste("PC1 (",PCA3PC1,"%)", sep = ""))+ ylab(paste("PC2 (",PCA3PC2,"%)", sep = "")) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP))+
guides(fill = "none", alpha="none",size="none", shape=guide_legend(order=2) ,color="none")+ 
#xlim(-110,65)+ ylim(-40,90) +
geom_dl(alpha=1,method="smart.grid", aes(label=POP, color=POP))+
theme_bw() + theme(text = element_text(size=12)) + theme(legend.position = "right", legend.title.align=0.5, legend.text.align = 0.5 )
```

``` r
gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.wildAE$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.wildAE$IND
gp.popmap$pop <- popmap.wildAE$POP



Chrom_Outlier_GP_Plot <- function(chrom,ncbi,vcf,prefix) {
plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("IsleofDogs1"))


outlier.gp <- genotype_plot(vcf=vcf,chr= ncbi,popmap = gp.popmap, start=1, end=104168038,cluster=TRUE,colour_scheme=plot_pal) 


outlier.gp$genotypes <- outlier.gp$genotypes  + guides(fill = "none") + theme(plot.title = element_text(margin=margin(0,0,-36,0))) 
outlier.gp$genotypes <- outlier.gp$genotypes + ggtitle(paste("Chr. ", chrom))
new_df <- data.frame(matrix(ncol=1,nrow=36, dimnames=list(NULL, c("IND"))))
new_df$IND <- outlier.gp$dendro_labels
new_df <- join(new_df, popmap.ni)

dendro_with_tips.outlier <- outlier.gp$dendrogram + 
geom_jitter(data=new_df,aes(x=1:length(IND),y=-2.5,color=POP,fill=POP,shape=Type),size=2, alpha=0.95, height=1.5, width = 0)+
scale_shape_manual(values=c(21),name="Origin") +
scale_fill_manual(values=p_wildAE, name="Population") + scale_color_manual(values=p_wildAE , name="Population") +
guides(alpha="none",size="none" , shape="none",fill="none",color="none")+
theme(legend.position = "none", legend.title=element_text(size=12))


layout <-"
AAABB#
"

fig <- dendro_with_tips.outlier + outlier.gp$genotypes+ plot_layout(design = layout) +  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

fig <- fig +  plot_layout(tag_level = "new")

assign(paste("gp.plot",chrom,prefix,"out", sep="."),fig, envir = globalenv())



}

Chrom_Outlier_GP_Plot_legend <- function(chrom,ncbi,vcf,prefix, start=1, end=104168038, bool=TRUE) {
plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("IsleofDogs1"))

outlier.gp <- genotype_plot(vcf=vcf,chr= ncbi,popmap = gp.popmap, start=start, end=end,cluster=TRUE,colour_scheme=plot_pal) 

outlier.gp$genotypes <- outlier.gp$genotypes + theme(legend.text = element_text(size=6),legend.position = "bottom", legend.title=element_text(size=8), plot.title = element_text(margin=margin(0,0,-36,0))) + guides(fill = guide_legend(override.aes = list(size=2))) +ggtitle(paste("Chr. ", chrom))
new_df <- data.frame(matrix(ncol=1,nrow=36, dimnames=list(NULL, c("IND"))))
new_df$IND <- outlier.gp$dendro_labels
new_df <- join(new_df, popmap.ni)

dendro_with_tips.outlier <- outlier.gp$dendrogram + 
  geom_jitter(data=new_df,aes(x=1:length(IND),y=-2.5,color=POP,fill=POP,shape=Type),size=2, alpha=0.95, height=1.5, width = 0)+
  scale_shape_manual(values=c(21),name="Origin") +
scale_fill_manual(values=p_wildAE, name="Population") + scale_color_manual(values=p_wildAE , name="Population") +
  guides(alpha="none",size="none" , shape="none",fill="none",color="none")+
  theme(legend.position = "none", plot.title = element_blank(), legend.title=element_text(size=12))


if (bool == FALSE){

layout <-"
AAABB#
CCCCCC
"
fig <- dendro_with_tips.outlier + outlier.gp$genotypes + guide_area() + plot_layout(design = layout, guides="collect", heights=c(10,1)) &  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

fig <- fig +  plot_layout(tag_level = "new")

}else{
layout <-"
AAABB#
"

fig <- dendro_with_tips.outlier + outlier.gp$genotypes+ plot_layout(design = layout) +  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

fig <- fig +  plot_layout(tag_level = "new")  
  
}


assign(paste("gp.plot",chrom,prefix,"out", sep="."),fig, envir = globalenv())


}
```

``` r
vcf1 <- "./Population_Genomics/Top.wildae.outliers.large_inversions.vcf.gz"
vcf2 <- "./Population_Genomics/Top.wildae.outliers.small_inversions.vcf.gz"
vcf3 <- "./Population_Genomics/Top.wildae.outliers.no_inversions.vcf.gz"


Chrom_Outlier_GP_Plot_legend(chroms[5],ncbis[5],vcf1,"vcf1", 1, 64020000, FALSE) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 603
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 603
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 603
    ##   row_num: 0
    ## Processed variant: 603
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot_legend(chroms[5],ncbis[5],vcf1,"vcf1a",64030000, 81090000, FALSE ) 
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 152
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 152
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 152
    ##   row_num: 0
    ## Processed variant: 152
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
gp.plot.5.vcf1a.out[[2]] <- gp.plot.5.vcf1a.out[[2]] + ggtitle("Chr. 5 (6)")
Chrom_Outlier_GP_Plot_legend(chroms[6],ncbis[6], vcf1, "vcf1",1,104168038, FALSE)
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 142
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 142
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 142
    ##   row_num: 0
    ## Processed variant: 142
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot_legend(chroms[1],ncbis[1], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 13
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 13
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 13
    ##   row_num: 0
    ## Processed variant: 13
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[2],ncbis[2], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 39
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 39
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 39
    ##   row_num: 0
    ## Processed variant: 39
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[3],ncbis[3], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 48
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 48
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 48
    ##   row_num: 0
    ## Processed variant: 48
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[4],ncbis[4], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 79
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 79
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 79
    ##   row_num: 0
    ## Processed variant: 79
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[5],ncbis[5], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 14
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 14
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 14
    ##   row_num: 0
    ## Processed variant: 14
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[9],ncbis[9], vcf2, "vcf2")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 9
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 9
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 9
    ##   row_num: 0
    ## Processed variant: 9
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot_legend(chroms[1],ncbis[1], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 27
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 27
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 27
    ##   row_num: 0
    ## Processed variant: 27
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[2],ncbis[2], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 68
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 68
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 68
    ##   row_num: 0
    ## Processed variant: 68
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[3],ncbis[3], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 137
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 137
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 137
    ##   row_num: 0
    ## Processed variant: 137
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[4],ncbis[4], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 92
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 92
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 92
    ##   row_num: 0
    ## Processed variant: 92
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[5],ncbis[5], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 18
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 18
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 18
    ##   row_num: 0
    ## Processed variant: 18
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[7],ncbis[7], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 18
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 18
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 18
    ##   row_num: 0
    ## Processed variant: 18
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
Chrom_Outlier_GP_Plot(chroms[9],ncbis[9], vcf3, "vcf3")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 103
    ##   header_line: 104
    ##   variant count: 4
    ##   column count: 45
    ## Meta line 103 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 4
    ##   Character matrix gt cols: 45
    ##   skip: 0
    ##   nrows: 4
    ##   row_num: 0
    ## Processed variant: 4
    ## All variants processed

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

    ## Joining by: IND

``` r
layout <- "
#AAABBBGG
#AAABBBGG
CCCCCCCCC
DDDEEEFFF
DDDEEEFFF
"

png(filename=paste(prefix.OF,"Figure.5.WildAe.FST.outliers.across.architecture.png", sep=""), type="cairo",units="px", width=6220, height=3500, res=300, bg="transparent")
finished_AE + pca_wildAE_12p + m2o.mALL +pca.out.plot1  +  pca.out.plot2  + pca.out.plot3   + guide_area() +plot_layout(design=layout, guides = "collect")  + plot_annotation(tag_levels = 'A')
dev.off()
```

    ## png 
    ##   2

``` r
layout<-"
ABBDDEE
#CC####
     "
sub_pp <- gp.plot.5.vcf1.out | (gp.plot.5.vcf1a.out | gp.plot.6.vcf1.out ) + plot_layout(heights = c(100,1))
sub_pp <- sub_pp + plot_layout(tag_level = "new")

layout <- "
AAAAAAAAA#
BBBBBBBCCE
BBBBBBBCCE
BBBBBBBDDE
BBBBBBBDDE
"

pca_wildAE_12pp <- pca_wildAE_12p + guides(shape="none")

plot_fin_1 <- m2o.m1 + wrap_elements(full=sub_pp) +  pca_wildAE_12pp  +   pca.out.plot1+  guide_area()+  plot_layout(design=layout, guides = "collect") + theme(legend.text =element_text(size=10)) + plot_annotation(tag_levels = 'A')


png(filename=paste(prefix.OFS,"Figure.S12.WildAe.FST.Large_inversions.outliers.png", sep=""), type="cairo",units="px", width=6800, height=3800, res=300, bg="transparent")
  plot_fin_1
  dev.off()
```

    ## png 
    ##   2

``` r
  layout <-"
  AAAAAAAAA
  CDEFGHIIJ
  CDEFGHBBJ
  "
    
plot_fin_2 <-m2o.m2 +  pca.out.plot2 + gp.plot.1.vcf2.out + gp.plot.2.vcf2.out + gp.plot.3.vcf2.out + gp.plot.4.vcf2.out + gp.plot.5.vcf2.out + gp.plot.9.vcf2.out + pca_wildAE_12pp + guide_area() +  plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A') + theme(legend.text =element_text(size=10)) 
  

plot_sub <- (gp.plot.1.vcf2.out | gp.plot.2.vcf2.out | gp.plot.3.vcf2.out | gp.plot.4.vcf2.out | gp.plot.5.vcf2.out | gp.plot.9.vcf2.out)/guide_area() + plot_layout(tag_level = 'keep',guides = "collect", heights = c(50, 1)) 

   layout <-"
    AAAAAAAAAA
    BBBBBBBCCE
    BBBBBBBCCE
    BBBBBBBDDE
    BBBBBBBDDE
    "     
        
plot_fin_2 <-m2o.m2 + wrap_elements(full=plot_sub)+  pca_wildAE_12pp +  pca.out.plot2  +guide_area()  +  plot_layout(design=layout, guides = "collect", tag_level = 'new') + plot_annotation(tag_levels = 'A') + theme(legend.text =element_text(size=10)) 
    
    
png(filename=paste(prefix.OFS,"Figure.S13.WildAe.FST.small_inversions_outliers.png", sep=""), type="cairo",units="px", width=6120, height=3420, res=300, bg="transparent")
      plot_fin_2
      dev.off()  
```

    ## png 
    ##   2

``` r
plot_sub <- (gp.plot.1.vcf3.out | gp.plot.2.vcf3.out | gp.plot.3.vcf3.out | gp.plot.4.vcf3.out | gp.plot.5.vcf3.out | gp.plot.7.vcf3.out | gp.plot.9.vcf3.out )/guide_area() + plot_layout(tag_level = 'keep',guides = "collect", heights = unit(c(500, 0.25), c('null', 'cm'))) 



   layout <-"
    AAAAAAAAAA
    BBBBBBBCCE
    BBBBBBBCCE
    BBBBBBBDDE
    BBBBBBBDDE
    "     
  
plot_fin_3 <-m2o.m3 +  wrap_elements(full=plot_sub) +  pca_wildAE_12pp + pca.out.plot3  + guide_area() +  plot_layout(design=layout, guides = "collect") + plot_annotation(tag_levels = 'A')+ theme(legend.text =element_text(size=10)) 

png(filename=paste(prefix.OFS,"Figure.S14.WildAe.FST.no_inversion_outliers.png", sep=""), type="cairo",units="px", width=6120, height=3420, res=300, bg="transparent")
  plot_fin_3
  dev.off()  
```

    ## png 
    ##   2

## FST genome scan using pairwise contrasts across salinity gradient contrasts in two Atlantic estuaries
