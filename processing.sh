######################
# ACTIVATE ENVIRONMENT
######################
conda create --name pacbio
conda activate pacbio
conda install -c bioconda raven-assembler
conda install -c bioconda prokka 
conda install -c bioconda minipolish
conda install -c bioconda flye
conda install -c bioconda krona
conda install -c bioconda seqtk
conda install -c bioconda miniasm 
conda install -c bioconda quast 
conda install -c bioconda bbmap 
conda install -c bioconda mmseqs2 

##################
# INSTALL METAMAPS
##################
git clone https://github.com/DiltheyLab/MetaMaps.git 
cd MetaMaps
./bootstrap.sh
./configure --prefix=/home/allie/src/ --with-boost=/usr/lib/x86_64-linux-gnu
sed -i 's/\/\/lib//' Makefile # edit makefile so it can find correct libraries
make metamaps
cp metamaps /home/allie/src/miniconda3/envs/pacbio/bin/
# download miniSeq+H.tar.gz database from website, upload to database folder
# https://www.dropbox.com/s/6p9o42yufx2vkae/miniSeq%2BH.tar.gz?dl=0
tar xzf miniSeq+H.tar.gz

#########
# MAPPING
#########
ls *fastq.gz | while read line; do metamaps mapDirectly --all -r ~/refdb/databases/miniSeq+H/DB.fa -t 55 -q $line -o $line.classification --maxmemory 400 1>$line.out 2>$line.err; done

################
# CLASSIFICATION
################
ls *classification | sed 's/.classification//' | while read line; do metamaps classify -t 55 --mappings $line.classification --DB ~/refdb/databases/miniSeq+H 1>$line.map.out 2>$line.map.err; done

###############
# PLOT TAXONOMY
###############
# download and update taxonomy files
ktUpdateTaxonomy.sh
bash /home/allie/src/miniconda3/envs/pacbio/opt/krona/updateTaxonomy.sh
# plot metamaps taxonomy output
ls *krona | sed 's/.fastq.gz.classification.EM.reads2Taxon.krona//' | parallel 'ktImportTaxonomy {}.fastq.gz.classification.EM.reads2Taxon.krona -o {}.krona.html'

##############
# HUMAN FILTER
##############
# get non human read ids
ls *reads2Taxon | sed 's/.fastq.gz.*//' | while read line; do grep -vE "9606|9605" $line.fastq.gz.classification.EM.reads2Taxon | awk '{print $1}' > $line.nonhum.ids; done
# pull only non human reads (fastq)
seqtk subseq demultiplex.bc1009_BAK8A_OA--bc1009_BAK8A_OA.hifi_reads.fastq.gz metamaps/demultiplex.bc1009_BAK8A_OA--bc1009_BAK8A_OA.hifi_reads.nonhum.ids > bc1009.nonhum.fq
seqtk subseq demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz metamaps/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.nonhum.ids > bc1011.nonhum.fq
# cat together
cat bc1009.nonhum.fq bc1011.nonhum.fq > all.fq
# convert to fasta
seqtk seq -a all.fq > all.fa
ls bc* | parallel 'gzip {}' &

##########
# ASSEMBLY
##########
# first just trying a full assembly (with all non-human data)
# raven
mkdir raven_assembled && cd raven_assembled
raven -t 55 ../all.fa 1>raven.out 2>raven.err
# minipolish
mkdir minipolish_assembled && cd minipolish_assembled
minimap2 -t 55 -x ava-ont ../all.fq ../all.fq > overlaps.paf-
miniasm -f ../all.fq overlaps.paf > assembly.gfa
minipolish -t 55 --pacbio ../all.fq assembly.gfa > polished.gfa
# convert from gfa to fasta
awk '/^S/{print ">"$2"\n"$3}' polished.gfa | fold > polished.fa
# flye
mkdir flye_assembled && cd flye_assembled
flye --pacbio-raw ../all.fq -t 55 -o ~/pacbio-6aug21/flye_assembled/ --meta --keep-haplotypes

###################
# EVALUATE ASSEMBLY
###################
quast-download-silva
metaquast.py -o metaquast -k --circos --threads 55 --fragmented flye_assembled/assembly.fasta minipolish_assembled/polished.fa raven_assembled/raven.out

#############
# BIN CONTIGS
#############
# install
conda install -c bioconda maxbin2
# use script from graphbin2 to get edge fasta sequences
wget https://raw.githubusercontent.com/Vini2/GraphBin2/master/support/gfa2fasta.py
chmod +x gfa2fasta.py
# add import subprocess
python3 ./gfa2fasta.py --graph assembly_graph.gfa --assembler flye --output gfa2fasta
# remap raw reads to contigs
mkdir bbmap && cd bbamp
mapPacBio.sh ref=~/pacbio-6aug21/flye_assembled/gfa2fasta/edges.fasta in=../all.fa out=bbmap.sam covstats=coverage.stats nodisk ignorebadquality
# get coverage statistics
pileup.sh in=bbmap.sam out=stats.txt
# bin with maxbin2
run_MaxBin.pl -contig gfa2fasta/edges.fasta -abund ~/pacbio-6aug21/bbmap/stats.txt -thread 55 -out ../maxbin

###################################
# ASSIGN TAXONOMY TO BINNED CONTIGS
###################################
conda install -c bioconda cat
wget https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz
tar -xvzf CAT_prepare_20210107.tar.gz
CAT contigs -c flye_assembled/assembly.fasta -d ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_CAT_database/ -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/
CAT add_names --force -i out.CAT.contig2classification.txt -o out.CAT.contig2classification.names.txt -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ --only_official
# bin contigs (not edges)
mapPacBio.sh ref=~/pacbio-6aug21/flye_assembled/assembly.fasta in=../all.fa out=bbmap.sam covstats=coverage.stats nodisk ignorebadquality
pileup.sh in=bbmap.sam out=stats.txt overwrite=true
# bin with maxbin2
run_MaxBin.pl -contig ~/pacbio-6aug21/flye_assembled/assembly.fasta -abund ~/pacbio-6aug21/bbmap/stats.txt -thread 55 -out ../maxbin/maxbin 1>../maxbin/maxbin.out
# assign taxonomy to bins
CAT bins --force --compress --sensitive -f 0.5 -b maxbin/ -d ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_CAT_database/ -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ -s .fasta
CAT add_names --force -i out.BAT.bin2classification.txt -o out.BAT.bin2classification.names.txt -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ --only_official
CAT summarise -i out.BAT.bin2classification.names.txt -o out.BAT.summary.txt

##############################
# ANNOTATE GENOMES WITH PROKKA
##############################
prokka --setupdb
prokka --listdb
mkdir prokka_annotated
ls *fasta | sed 's/.fasta//' | while read line; do prokka --outdir ../prokka_annotated/$line/ --prefix $line.annotated --addgenes --cpus 0  $line.fasta; done

##############################################
# MAP SAMPLES SEPARATELY TO EACH ANNOTATED BIN
##############################################
# bc1009 : CF-PF
ls | while read line; do mapPacBio.sh ref=$line\/$line.annotated.fna in=../bc1009.nonhum.fa.gz out=$line.bc1009.sam covstats=$line.bc1009.cov.stats nodisk ignorebadquality; done
ls *sam | sed 's/.bc1009.sam//' | while read line; do mapPacBio.sh ref=$line\/$line.annotated.fna in= ../bc1011.nonhum.fa.gz out=$line.bc1011.sam covstats=$line.bc1011.cov.stats nodisk ignorebadquality; done

############
# SAM TO BAM
############
# problem with conda version of samtools?
conda install -c bioconda samtools=1.9 --force-reinstall
# sam to bam
ls *sam | sed 's/.sam//' | parallel 'samtools view -bS {}.sam > {}.bam'
rm maxbin*sam 
# sort bam files
ls *bam | sed 's/.bam//' | parallel 'samtools view {}.bam -o {}.sort.bam'
# clean up
ls maxbin*bam | grep -v "sort" | parallel 'rm {}'

#################################
# PREDICTED TAXONOMY FOR EACH BIN
#################################
# pull rpoc sequence ids
ls -d */ | sed 's/\///' | while read line; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.ffn/,x)}1' $line\/$line.annotated.ffn | grep "DNA-directed RNA polymerase subunit beta'" | sed 's/\// /' | sed 's/>//' | sed 's/ maxbin.*annotated_/\t/' | sed 's/ /\t/'; done > rpoc_bin.list
cat maxbin.*/*ffn > all_bins.ffn
awk '{print $2}' rpoc_bin.list > rpoc.ids
seqtk subseq all_bins.ffn rpoc.ids > rpoc.ffn
# assign taxonomy to each read with kraken2
conda install -c bioconda kraken2 
# install standard kraken db
kraken2-build --standard --db ~/refdb/kraken2_standard_db --threads 60
# run kraken2
kraken2 --db ~/refdb/kraken2_standard_db/ --threads 40 --use-names --output rpoc.kraken.out rpoc.ffn --unclassified-out rpoc.unclassified.kraken.out --confidence 0.01
# merge taxonomic information with bin information to get annotations file
awk -F"\t" '{print $2 "\t" $3}' rpoc.kraken.out > taxonomy.txt
paste rpoc_bin.list taxonomy.txt | sed '1 i\bin_number\tlocus_tag\tfunction\tlocus_tag_2\tpredicted_taxonomy'> annotations.txt

##########################
# BUILD RPOC TAXONOMY TREE
########################## 
mafft rpoc.ffn > rpoc.align.fa
raxmlHPC-PTHREADS-SSE3 -T 40 -c 25 -m GTRCAT -e 0.001 -p 31415 -f a -N 100 -x 02938 -n rpoc.tre -s rpoc.align.fa


#################
# SAMPLE ASSEMBLY
#################
# want to assemble data from each sample separately, then taxonomically assign/map to reference genome -- different strains?
cd sample_assembly
mkdir bc1009
mkdir bc1011
flye --pacbio-raw ../demultiplex.bc1009_BAK8A_OA--bc1009_BAK8A_OA.hifi_reads.fastq.gz -t 28 -o bc1009/ --meta --keep-haplotypes
flye --pacbio-raw ../demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz -t 28 -o bc1011/ --meta --keep-haplotypes

###################
# EVALUATE ASSEMBLY
###################
metaquast -o metaquast -k --circos --threads 55 --fragmented bc1009/assembly.fasta bc1011/assembly.fasta

############################
# ASSIGN TAXONOMY TO CONTIGS
############################
CAT contigs -c bc1009/assembly.fasta -d ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_CAT_database/ -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/
mkdir cat_bc1009
mv out.CAT.* cat_bc1009/
cd cat_bc1009
CAT add_names --force -i out.CAT.contig2classification.txt -o out.CAT.contig2classification.names.txt -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ --only_official
cd ..
# do bc1011
CAT contigs -c bc1011/assembly.fasta -d ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_CAT_database/ -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/
mkdir cat_bc1011
mv out.CAT.* cat_bc1011
cd cat_bc1011/
CAT add_names --force -i out.CAT.contig2classification.txt -o out.CAT.contig2classification.names.txt -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ --only_official
# sumarize results
cd cat_bc1009
CAT summarise -i out.CAT.contig2classification.names.txt -o out.CAT.summary.txt -c ../bc1009/assembly.fasta
cd ../cat_bc1011
CAT summarise -i out.CAT.contig2classification.names.txt -o out.CAT.summary.txt -c ../bc1011/assembly.fasta

#############
# BIN CONTIGS
#############
# remap raw reads to contigs
mkdir bbmap && cd bbamp
mapPacBio.sh ref=bc1009/assembly.fasta in=../bc1009.nonhum.fa.gz out=bbmap/bbmap_1009.sam covstats=bbmap/coverage_1009.stats nodisk ignorebadquality
mapPacBio.sh ref=bc1011/assembly.fasta in=../bc1011.nonhum.fa.gz out=bbmap/bbmap_1011.sam covstats=bbmap/coverage_1011.stats nodisk ignorebadquality
# get coverage statistics
pileup.sh in=bbmap_1009.sam out=bbmap_1009.stats
pileup.sh in=bbmap_1011.sam out=bbmap_1011.stats
cd ..
# bin with maxbin2
mkdir maxbin
run_MaxBin.pl -contig bc1009/assembly.fasta -abund bbmap/bbmap_1009.stats -thread 60 -out maxbin/bc1009
run_MaxBin.pl -contig bc1011/assembly.fasta -abund bbmap/bbmap_1011.stats -thread 55 -out maxbin/bc1011

###################################
# ASSIGN TAXONOMY TO BINNED CONTIGS
###################################
CAT bins --force --compress --sensitive -f 0.5 -b maxbin/ -d ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_CAT_database/ -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ -s .fasta
CAT add_names --force -i out.BAT.bin2classification.txt -o out.BAT.bin2classification.names.txt -t ~/refdb/cat_ref/CAT_prepare_20210107/2021-01-07_taxonomy/ --only_official
CAT summarise -i out.BAT.bin2classification.names.txt -o out.BAT.summary.txt

##############################
# ANNOTATE GENOMES WITH PROKKA
##############################
prokka --setupdb
prokka --listdb
mkdir prokka_annotated
ls *fasta | sed 's/.fasta//' | while read line; do prokka --outdir ../prokka_annotated/$line/ --prefix $line.annotated --addgenes --cpus 0  $line.fasta; done





