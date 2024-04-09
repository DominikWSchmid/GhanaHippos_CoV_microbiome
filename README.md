# SARS-related beta-coronavirus infection linked to gut microbial dysbiosis in bats

We included here all primary 16S and meta data collected to arrive at the conclusions from the above mentioned manuscript. The fecal samples were collected from cave-roosting bats in Ghana over two years  between 2010 and 2012 (see methods for details). The fecal samples were used for viral and 16S microbiota screening. We found that infections with one SARS-related beta-CoV altered the gut microbial composition of its bat host and did so in a way that more intensely infected bats suffered a more heavily dysbiotic microbial configuration. 

# Description of the data and file structure

The RDS file includes a pre-filtered phyloseq object, which includes a tax_table (describing the taxonomy of described ASVs), otu_table (describing the reads assigned to each ASV per sample), phy_tree (describing the phylogenetic relationship between ASVs) and the sample_data (which contains all meta information about each sample).

The meta data contains the following columns in no particular order ("sample.ID"  [individualised identifier for the sample], "TissueID" [identifier for the tissue sample to identify bat species], "FaecalID" [identifier for the fecal sample for the respective sampled bats], "X229ELogical" [category to determine infection status with CoV-229E-like], "CT_value229E" [cycle threshold value for respective 229E infection, NAs imply CT-value above 44.00, i.e., negative], "X2bLogical" [category to determine infection status with CoV-2b], "CT_value2B" [cycle threshold value for respective 2b infection, NAs imply CT-value above 44.00, i.e., negative], "X2bBasalLogical" [category to determine infection status with CoV-2bBasal], "CT_value2Bbasal" [cycle threshold value for respective 2bBasal infection, NAs imply CT-value above 44.00, i.e., negative], "Lineage" [host bat species, i.e., H. caffer D], "Site" [sampling cave site], "AGE" [age category], "SEX" [sex category], "sample_period" [sampling period during which the sample was collected], "infection_status" [infection status as category, i.e., uninfected, 229E-infected, 2b-infected, 2bBasal-infected]

# Sharing/Access information

Links to other publicly accessible locations of the data:

[http://...](https://doi.org/10.5061/dryad.sxksn039t)

NCBI accession numbers of raw microbiota fasta files are:

Bioproject-PRJNA1096136

# Code/Software

As described in the methods, raw reads were processed with the DADA2 plug-in in QIIME 2 (v2021.8.0).

Subsequent quality control and filtering steps of the phylseq object were performed in R Studio. 

All statistical analysis of the cleaned phyloseq object was also performed in R Studio and are described in the code included in this submission. Versions of packages and parameterisation was described in the methods of the manuscript/supplementary material. In case of questions, do not hesitate to contact dominikwerner.schmid@uni-ulm.de
