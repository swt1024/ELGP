# ELGP: A computational framework for predicting essential lncRNA genes with lncRNA-protein-protein interaction network and multi-omics data.

In this work, we proposed a framework to identify essential lncRNAs by taking advantage of the topological feature of the lncRNA-protein-protein heterogeneous network and muti-omics features. We used multi-omics features to annotate nodes at lncRNA-protein-protein interaction network and introduced the HinSAGE algorithm to learn node representation for lncRNAs in the lncRNA-protein-protein heterogeneous network. We named this method as ELGP.

## 1. File descriptions

### 1.1 Data

Some of the data in the `data` folder has been preprocessed. For more details, please refer to the *Preprocessing* section. For some datasets that require substantial memory, only download links are provided.

#### 1.1.1 Benchmarking

- `esslnc.csv`: Essential lncRNA gene data downloaded from the dbEssLnc2.0 database.

**human**
- `esslnc.csv`: Human essential lncRNA gene data from the dbEssLnc2.0 database.
- `ess_lpi.csv`: Human essential lncRNA genes identified from the LPI network.
- `noness_lpi.csv`: Selected human non-essential lncRNA genes from the LPI network.

**mouse**
- `esslnc.csv`: Mouse essential lncRNA gene data from the dbEssLnc2.0 database.
- `ess_lpi.csv`: Mouse essential lncRNA genes identified from the LPI network.
- `noness_lpi.csv`: Selected mouse non-essential lncRNA genes from the LPI network.

#### 1.1.2 LncBook_LPI

- [lncrna_rbp_LncBook2.0.csv.gz](https://ngdc.cncb.ac.cn/lncbook/files/lncrna_rbp_LncBook2.0.csv.gz): Raw human LPI data downloaded from the LncBook 2.0 database.  
- `LncBook_LPI.csv`: Processed human LPI data derived from the LncBook 2.0 database.

#### 1.1.3 NPInter_LPI

- [lncRNA_interaction.txt.gz](http://bigdata.ibp.ac.cn/npinter5/download/file/lncRNA_interaction.txt.gz): Raw LPI dat download from the NPInterv5.0 database. 

**human**
- `NPInter_LPI.csv`: Filtered human LPI data.
- `correct_NPInter_LPI.csv`: Corrected human LPI data.

**mouse**
- `NPInter_LPI.csv`: Filtered mouse LPI data.
- `correct_NPInter_LPI.csv`: Corrected mouse LPI data.

#### 1.1.4 LPI

**human**
- `LPI.csv`: Merged and deduplicated human LPI data from the LncBook and NPInter databases.
- `lncRNA.csv`: lncRNA nodes in the human LPI network.
- `protein.csv`: Protein nodes in the human LPI network.

**mouse**
- `LPI.csv`: Mouse LPI data from the NPInter database.
- `lncRNA.csv`: lncRNA nodes in the mouse LPI network.
- `protein.csv`: Protein nodes in the mouse LPI network.

#### 1.1.5 PPI

**human**
- `BIOGRID-Homo_sapiens-4.4.241.tab3.txt`: Human PPI data from the BIOGRID database, which  were obtained from the [BIOGRID-ORGANISM-4.4.241.tab3.zip](https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.241/BIOGRID-ORGANISM-4.4.241.tab3.zip).
- `PPI.csv`: Filtered human PPI data.
- `protein_in_ppi.csv`: Proteins involved in the human PPI data.

**mouse**
- `BIOGRID-Mus_musculus-4.4.241.tab3.txt`: Mouse PPI data from the BIOGRID database, which  were obtained from the [BIOGRID-ORGANISM-4.4.241.tab3.zip](https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.241/BIOGRID-ORGANISM-4.4.241.tab3.zip).
- `PPI.csv`: Filtered mouse PPI data.
- `protein_in_ppi.csv`: Proteins involved in the mouse PPI data.

#### 1.1.6 LPPI

**human**
- `LPPI.csv`: Human LPPI data obtained by merging LPI and PPI networks.
- `protein.csv`: Protein nodes in the human LPPI network.

**mouse**
- `LPPI.csv`: Mouse LPPI data obtained by merging LPI and PPI networks.
- `protein.csv`: Protein nodes in the mouse LPPI network.
- `LPPI_updated.csv`: Mouse LPPI data after correcting invalid protein names.
- `protein_updated.csv`: Protein nodes in the corrected mouse LPPI network.

#### 1.1.7 reference_lncRNA

**human**
**`gtf`**
- [lncRNA_LncBookv2.0_GRCh38.gtf](https://ngdc.cncb.ac.cn/lncbook/files/lncRNA_LncBookv2.0_GRCh38.gtf.gz): Human lncRNA annotations from LncBook v2.0, GRCh38 assembly.  
- [NONCODEv5_human_hg38_lncRNA.gtf](http://v5.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz): Human lncRNA annotations from NONCODE v5, hg38 assembly.  
- [NONCODEv6_human_hg38_lncRNA.gtf](http://www.noncode.org/datadownload/NONCODEv6_human_hg38_lncRNA.gtf.gz): Human lncRNA annotations from NONCODE v6, hg38 assembly.  
- `ensembl`: Contains multiple versions of GTF files downloaded from Ensembl. The download links for each version are as follows:
  - [Homo_sapiens.GRCh38.113.gtf.gz](https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz)
  - [Homo_sapiens.GRCh38.112.gtf.gz](https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz)
  - [Homo_sapiens.GRCh38.111.gtf.gz](https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz)
  - [Homo_sapiens.GRCh38.110.gtf.gz](https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz)
  - [Homo_sapiens.GRCh38.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz)
  - [Homo_sapiens.GRCh38.108.gtf.gz](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz)
  - [Homo_sapiens.GRCh38.107.gtf.gz](https://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz)
  - [Homo_sapiens.GRCh38.106.gtf.gz](https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz)
  - [Homo_sapiens.GRCh38.104.gtf.gz](https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz)
  - [Homo_sapiens.GRCh38.97.gtf.gz](https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz)
  - [Homo_sapiens.GRCh38.93.gtf.gz](https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz)
  - [Homo_sapiens.GRCh38.87.gtf.gz](https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz)
  - [Homo_sapiens.GRCh38.84.gtf.gz](https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz)
  - [Homo_sapiens.GRCh38.80.gtf.gz](https://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz)
  - [Homo_sapiens.GRCh38.78.gtf.gz](https://ftp.ensembl.org/pub/release-78/gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz)
  - [Homo_sapiens.GRCh38.76.gtf.gz](https://ftp.ensembl.org/pub/release-76/gtf/homo_sapiens/Homo_sapiens.GRCh38.76.gtf.gz)

**`bed`**
- `pro_bed.sh`: Script to process `lncRNA_LncBookv2.0_GRCh38.gtf`, [NONCODEv6_hg38.lncAndGene.bed.gz](http://www.noncode.org/datadownload/NONCODEv6_hg38.lncAndGene.bed.gz), and [NONCODEv5_hg38.lncAndGene.bed.gz](http://v5.noncode.org/datadownload/NONCODEv5_hg38.lncAndGene.bed.gz) to generate the NONCODE BED files (`NONCODEv6_hg38.lncRNAGene.bed`, `NONCODEv5_hg38.lncRNAGene.bed`) and the LncBook BED file (`lncRNA_gene_LncBookv2.0_GRCh38.bed`) used in this study. 
- `gtf2bed.py`: Script to process multiple versions of GTF files downloaded from the Ensembl database to generate corresponding BED files.
- `ensembl`: Folder containing the generated Ensembl BED files.

**`transcript`**
- `get_trans.py`: Script to process GTF files from NONCODE and LncBook, generating mappings between lncRNA genes and transcripts. This results in the creation of `lncRNA_LncBookv2.0_GRCh38_trans.txt`, `NONCODEv6_human_hg38_lncRNA_trans.txt`, and `NONCODEv5_human_hg38_lncRNA_trans.txt`.
- `get_trans_ensembl.py`: Script to process GTF files from Ensembl, generating mappings between lncRNA genes and transcripts for various versions.
- `ensembl`: Folder that stores mappings between genes and transcripts for various versions.

**mouse**
**`gtf`**
- [NONCODEv5_mouse_mm10_lncRNA.gtf](http://v5.noncode.org/datadownload/NONCODEv5_mouse_mm10_lncRNA.gtf.gz): Mouse lncRNA annotations from NONCODE v5, mm10 assembly.  
- `ensembl`: Contains multiple versions of GTF files downloaded from Ensembl. The download links for each version are as follows:
  - [Mus_musculus.GRCm38.100.gtf](https://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz)
  - [Mus_musculus.GRCm38.97.gtf](https://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz)
  - [Mus_musculus.GRCm38.96.gtf](https://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz)
  - [Mus_musculus.GRCm38.93.gtf](https://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz)
  - [Mus_musculus.GRCm38.87.gtf](https://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz)
  - [Mus_musculus.GRCm38.86.gtf](https://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz)
  - [Mus_musculus.GRCm38.85.gtf](https://ftp.ensembl.org/pub/release-85/gtf/mus_musculus/Mus_musculus.GRCm38.85.gtf.gz)
  - [Mus_musculus.GRCm38.84.gtf](https://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz)

**`bed`**
- `pro_bed.sh`: Script to process [NONCODEv6_mm10.lncAndGene.bed](http://www.noncode.org/datadownload/NONCODEv6_mm10.lncAndGene.bed.gz), and [NONCODEv5_mm10.lncAndGene.bed](http://v5.noncode.org/datadownload/NONCODEv5_mm10.lncAndGene.bed.gz) to generate the NONCODE BED files (`NONCODEv6_mm10.lncRNAGene.bed`, `NONCODEv5_mm10.lncRNAGene.bed`) used in this study. 
- `gtf2bed.py`: Script to process multiple versions of GTF files downloaded from the Ensembl database to generate corresponding BED files.
- `ensembl`: Folder containing the generated Ensembl BED files.

**`transcript`**
- `get_trans.py`: Script to process GTF files from NONCODE and LncBook, generating mappings between lncRNA genes and transcripts. This results in the creation of `lncRNA_LncBookv2.0_GRCh38_trans.txt`, `NONCODEv6_human_hg38_lncRNA_trans.txt`, and `NONCODEv5_mouse_mm10_lncRNA_trans.txt`.
- `get_trans_ensembl.py`: Script to process GTF files from Ensembl, generating mappings between lncRNA genes and transcripts for various versions.
- `ensembl`: Folder that stores mappings between genes and transcripts for various versions.

#### 1.1.8 Sequence
- [GRCh38.p14.genome.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz): Human reference genome sequence file for the GRCh38.p14 release.  
- [GRCm38.p6.genome.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz): Mouse genome sequence file for the GRCm38.p6 release.  
- `process_fasta.py`: Script for preprocessing fasta files.

**human**
- [LncBookv2_OnlyLnc.fa](https://ngdc.cncb.ac.cn/lncbook/files/lncRNA_LncBookv2.0.fa.gz): Human lncRNA sequence file from LncBook v2.0.  
- [NONCODEv6_human.fa](http://www.noncode.org/datadownload/NONCODEv6_human.fa.gz): Human lncRNA sequence file from NONCODE v6.  
- [NONCODEv5_human.fa](http://v5.noncode.org/datadownload/NONCODEv5_human.fa.gz): Human lncRNA sequence file from NONCODE v5.  

**`ensembl`** - Non-coding gene sequence files downloaded from Ensembl:
- [Homo_sapiens.GRCh38.113.ncrna.fa](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.112.ncrna.fa](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.111.ncrna.fa](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.110.ncrna.fa](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.109.ncrna.fa](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.108.ncrna.fa](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.107.ncrna.fa](https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.106.ncrna.fa](https://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.104.ncrna.fa](https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.97.ncrna.fa](https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.93.ncrna.fa](https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.87.ncrna.fa](https://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.84.ncrna.fa](https://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.80.ncrna.fa](https://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.78.ncrna.fa](https://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Homo_sapiens.GRCh38.76.ncrna.fa](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)

**`processed_ensembl`** - Folder containing preprocessed Ensembl sequence files.

**mouse**
- [NONCODEv5_mouse.fa](http://v5.noncode.org/datadownload/NONCODEv5_mouse.fa.gz): Mouse lncRNA sequence file from NONCODE v5.  
- [NONCODEv6_mouse.fa](http://www.noncode.org/datadownload/NONCODEv6_mouse.fa.gz): Mouse lncRNA sequence file from NONCODE v6.  

**`ensembl`** - Non-coding gene sequence files downloaded from Ensembl for mouse:
- [Mus_musculus.GRCm38.100.ncrna.fa](https://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.97.ncrna.fa](https://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.96.ncrna.fa](https://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.93.ncrna.fa](https://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.87.ncrna.fa](https://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.86.ncrna.fa](https://ftp.ensembl.org/pub/release-86/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.85.ncrna.fa](https://ftp.ensembl.org/pub/release-85/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)
- [Mus_musculus.GRCm38.84.ncrna.fa](https://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz)

**`processed_ensembl`** - Folder containing preprocessed Ensembl sequence files for mouse.


### 1.2 preprocess

#### LPI_human
- `pro_lncbook.sh`: Script to filter LPI data from LncBook, retaining only interactions between lncRNA genes and proteins.
- `pro_NPInter.sh`: Script to filter LPI data from NPInter, retaining only interactions found in humans where both interacting molecules are categorized as 'lncRNA' and 'protein', and the interaction calss is 'binding'.
- `getCoordination.py`: Script to retrieve gene coordinates from BED files in the `data/reference_lncRNA/human/bed` folder. Genes without corresponding genomic coordinates are recorded in the `lnc_no_pos.csv` file.
- `match_lncRNA.sh`: Uses bedtools to identify lncRNA genes with exactly overlapping genomic coordinates. Pairs of completely overlapping lncRNA genes are stored in the `mapped_lncRNA.txt` file.
- `merge_all.py`: Merges and deduplicates LPI data from LncBook and NPInter.
- `correct_inter.py`: Corrected erroneous protein-coding gene names, incorrect tissue or cell line names, and replaced transcript IDs with their corresponding gene IDs in the LPI data.

#### LPI_Mouse

#### LPPI

#### benchMarking_human

#### benchMarking_Mouse

### 1.3 annotate 

**human**
- lncRNA_feature_annotate.ipynb: Annotating features for lncRNA nodes in the heterogeneous graph.
- protein_feature_annotate.ipynb:Annotating features for protein nodes in the heterogeneous graph.

**mouse**
- lncRNA_feature_annotate.ipynb: Annotating features for lncRNA nodes in the heterogeneous graph.
- protein_feature_annotate.ipynb:Annotating features for protein nodes in the heterogeneous graph.

### 1.4 HinSAGE

- train.py:Training the HinSAGE model and generating node representations for lncRNA gene nodes

### 1.5 classifier

**SVM**
- svm.ipynb: Code for testing HinSAGE parameters, tuning SVM parameters, performing cross-validation, and making predictions.

**MLP**
- mlp.ipynb: Code for tuning MLP parameters, performing cross-validation, and making predictions.

## 2. How to use ELGP
The following steps take human as an example.

### 2.1 Construct lncRNA-protein-protein heterogeneous network

Step1: Users can run `idmap.py`, input `LPI.csv`  , output `id_lncRNA.txt`, `id_protein.txt` and `id_lncRNA_protein.txt`.

### 2.2 Heterogeneous representation learning

Step1:  Users can run `./HinSAGE/train.py` to generate the node representation for lncRNA gene nodes , e.g. `lncRNA_embeddings_heart`. 

### 2.3 Select optimal non-essential lncRNAs

### 2.4 Supervised machine learning

#### 2.4.1 SVM model

| dataset |   C   | kernel function |
| :-----: | :---: | :-------------: |
|  Human  |  100  |     linear      |
|  Mouse  |   10  |     linear      |


#### 2.4.2 MLP model

| dataset | The number of neurons in the first hidden layer | The number of neurons in the second hidden layer |
| :-----: | :---------------------------------------------: | :----------------------------------------------: |
|  Human  |                      256                        |                        64                        |
|  Mouse  |                       64                        |                        256                       |


