#### Run in bash terminal
#These steps do not include all the analysis steps, only the most important ones.
#### Results obtained from PfamScan
#### Domain used 
# ASIC family	PF00858.27
# Cys-loop family	PF02931.26
# Cys-loop family	PF02932.19
# iGluR family	PF00060.29
# iGluR family	PF10613.12
# P2X family	PF00864.22
#for l in *fa ; do  pfam_scan.pl -fasta $l -dir /home/vic/bin/Database/Pfam_35/Pfam_LGIC -outfile $l.dom -cpu 16 -e_seq 1 ; done

#for l in *dom; do
  #filename=${l}
  #base="${filename%.dom}"  
  #grep -v -e "^$" -e "#" $l | awk -v FS=' ' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' | sed '1iseq id\talignment start\talignment end\tenvelope start\tenvelope end\thmm acc\thmm name\ttype\thmm start\thmm end\thmm length\tbit score\tE-value\tsignificance\tclan' | awk '$12 >= 11' > 1.$l.mod ; cut -f 1 1.$l.mod | sort | uniq > 2.codigos.$l.pfam ; fastagrep.pl -X -f 2.codigos.$l.pfam $base > 3.codigos.$l.pfam.fa ; cd-hit -i 3.codigos.$l.pfam.fa -o 4.codigos.$l.pfam_cd_96.fa  -c 0.96 ; done 

#### Results obtained in Diamond
#### Files available in "Results_Diamons_taxa_Pfam_Trembl_annotation"
#### DIAMOND
#diamond blastp -d /Database/trembl/01_Metazoa_not_Protostomia/Metazoa_not_Protostomia.fasta -q LGIC.fa -o matches_01_Metazoa.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/02_Protostomia_not_Spiralia/Protostomia_not_Spiralia.fasta.dmnd -q LGIC.fa -o matches_02_Protostomia.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/03_Spiralia_not_Lophotrochozoa/Spiralia_not_Lophotrochozoa.fasta.dmnd -q LGIC.fa -o matches_03_Spiralia.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/04_Lophotrochozoa_not_Platyhelminthes/Lophotrochozoa_not_Platyhelminthes.fasta.dmnd -q LGIC.fa -o matches_04_Lophotrochozoa.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/05_Platyhelminthes_not_Neodermata/Platyhelminthes_not_Neodermata.fasta.dmnd -q LGIC.fa -o matches_05_Platyhelminthes_not_Neodermata.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/06A_Neodermata_not_Monopisthocotylea/Neodermata_not_Monopisthocotylea.fasta.dmnd -q LGIC.fa -o matches_06A_Neodermata_not_Monopisthocotylea.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/06B_Neodermata_not_Polyopisthocotylea/Neodermata_not_Polyopisthocotylea.fasta.dmnd -q LGIC.fa -o matches_06B_Neodermata_not_Polyopisthocotylea.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/06C_Neodermata_not_Trematoda/Neodermata_not_Trematoda.fasta.dmnd -q LGIC.fa -o matches_06C_Neodermata_not_Trematoda.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/06D_Neodermata_not_Cestoda/Neodermata_not_Cestoda.fasta.dmnd -q LGIC.fa -o matches_06D_Neodermata_not_Cestoda.tsv --sensitive -p 16
#diamond blastp -d /Database/trembl/06E_Neodermata_not_Free-living/Neodermata_not_Free-living.fasta.dmnd -q LGIC.fa -o matches_06E_Neodermata_not_Free-living.tsv --sensitive -p 16

#cut -f 1,11,12 matches_02_Protostomia.tsv |awk '!seen[$1]++' > matches_02_Protostomia_unique.tsv
#cut -f 1,11,12 matches_03_Spiralia.tsv |awk '!seen[$1]++' > matches_03_Spiralia_unique.tsv
#cut -f 1,11,12 matches_04_Lophotrochozoa.tsv |awk '!seen[$1]++' > matches_04_Lophotrochozoa_unique.tsv
#cut -f 1,11,12 matches_05_Platyhelminthes_not_Neodermata.tsv |awk '!seen[$1]++' > matches_05_Platyhelminthes_not_Neodermata_unique.tsv
#cut -f 1,11,12 matches_06A_Neodermata_not_Monopisthocotylea.tsv |awk '!seen[$1]++' > matches_06A_Neodermata_not_Monopisthocotylea_unique.tsv
#cut -f 1,11,12 matches_06B_Neodermata_not_Polyopisthocotylea.tsv |awk '!seen[$1]++' > matches_06B_Neodermata_not_Polyopisthocotylea_unique.tsv
#cut -f 1,11,12 matches_06C_Neodermata_not_Trematoda.tsv |awk '!seen[$1]++' > matches_06C_Neodermata_not_Trematoda_unique.tsv
#cut -f 1,11,12 matches_06D_Neodermata_not_Cestoda.tsv |awk '!seen[$1]++' > matches_06D_Neodermata_not_Cestoda_unique.tsv
#cut -f 1,11,12 matches_06E_Neodermata_not_Free-living.tsv |awk '!seen[$1]++' > matches_06E_Neodermata_not_Free-living_unique.tsv
#cat matches_06*_unique.tsv > matches_06_Neodermata_unique.tsv

####For anotation Tremb_UniProti
#diamond blastp -d ~/bin/Database/trembl_Metazoa.fasta.dmnd -q Platyhelminthes_cdhit_LGIC_filtered.fa -o matches_Trembl_for_annotation.tsv --sensitive -p 16
#cut -f 1,2,11,12 matches_Trembl_for_annotation.tsv | sed "s/tr|// ; s/sp|// ; s/|/\t/" | awk '!seen[$1]++' > matches_Trembl_for_annotation_unique.tsv



#### R Script for merging and processing the results
#setwd("script_plots/")
getwd()

library(tidyverse)

#Upload results files
codigos.LGIC <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/ID_Proteins_LGIC",sep = "\t", header = FALSE) %>% rename(ID = V1)
LGIC.pfam <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/LGIC.pfam_unique.tsv",sep = "\t", header = TRUE) %>% rename(ID = seq.id) %>% select(ID, Class, Species,LGIC.Family, hmm.acc, hmm.name, bit.score, E.value)
LGIC.trembl <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_Trembl_for_annotation_unique.tsv",sep = "\t", header = FALSE)  %>% rename(Entry = V2, ID = V1)  %>% select(ID, Entry)
trembl <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/for_annotation_trembl.tsv", header = TRUE)
match01 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_01_Metazoa_unique.tsv", header = FALSE) %>% rename(Match_Metazoa_score = V3, Match_Metazoa_evalue = V2, ID = V1)  %>% select(ID, Match_Metazoa_evalue)
match02 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_02_Protostomia_unique.tsv", header = FALSE) %>% rename(Match_Protostomia_score = V3, Match_Protostomia_evalue = V2, ID = V1) %>% select(ID, Match_Protostomia_evalue)
match03 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_03_Spiralia_unique.tsv", header = FALSE) %>% rename(Match_Spiralia_score = V3, Match_Spiralia_evalue = V2, ID = V1) %>% select(ID, Match_Spiralia_evalue)
match04 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_04_Lophotrochozoa_unique.tsv", header = FALSE) %>% rename(Match_Lophotrochozoa_score = V3, Match_Lophotrochozoa_evalue = V2, ID = V1) %>% select(ID, Match_Lophotrochozoa_evalue)
match05 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_05_Platyhelminthes_not_Neodermata_unique.tsv", header = FALSE) %>% rename(Match_Platyhelminthes_score = V3, Match_Platyhelminthes_evalue = V2, ID = V1)  %>% select(ID, Match_Platyhelminthes_evalue)
match06 <- read.delim("Results_Diamons_taxa_Pfam_Trembl_annotation/matches_06_Neodermata_unique.tsv", header = FALSE) %>% rename(Match_Neodermata_score = V3, Match_Neodermata_evalue = V2, ID = V1) %>% select(ID, Match_Neodermata_evalue)

#Join the results
Evalue <- left_join(codigos.LGIC, match01 , "ID")
Evalue <- left_join(Evalue, match02 , "ID")
Evalue <- left_join(Evalue, match03 , "ID")
Evalue <- left_join(Evalue, match04 , "ID")
Evalue <- left_join(Evalue, match05 , "ID")
Evalue <- left_join(Evalue, match06 , "ID")
Evalue <- left_join(Evalue, LGIC.pfam , "ID")
Evalue <- left_join(Evalue, LGIC.trembl , "ID")
Evalue <- left_join(Evalue, trembl, "Entry")


# Transform the E-values to -log10 scale
Evalue2 <- Evalue %>% mutate(across(starts_with("Match") & ends_with("evalue"), ~ -log10(.)))
# Replace infinite values resulting from log transformation with 308 (as -log10(0) is infinite, and 308 is a common upper bound for -log10(E-value))
Evalue2 <- Evalue2 %>% mutate(across(starts_with("Match") & ends_with("evalue"), ~ replace(., is.infinite(.), 308)))
# Round the transformed E-values to the nearest integer
Evalue2 <- Evalue2 %>%   mutate(across(starts_with("Match") & ends_with("evalue"), round))
# Cap the transformed E-values at 308 (in case any value exceeds this after rounding)
Evalue2 <- Evalue2 %>% mutate(across(starts_with("Match") & ends_with("evalue"), ~ pmin(., 308)))
# Relocate a specific column (the 10th column) to appear after the 28th column
Evalue2 <- Evalue2 %>%  relocate(names(Evalue2)[10], .after = names(Evalue2)[28])

Evalue2 <- Evalue2 %>% mutate(Class.LGIC = paste(Evalue2[[8]], Evalue2[[28]], sep = "-"))


#Save the results to make analyzes and graphs
write.table(Evalue2, file = "LGIC_Platyheminthes_2.tsv", sep = "\t", row.names = FALSE)
