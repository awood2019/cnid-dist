#analyzing hydrozoan eDNA sequences

#load packages
library(dada2)

#--------code from GoMx phyla distribution script to make phyloseq object-------

#load packages
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 
library("ggplot2")

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-phy-dist")

#load components of phyloseq object: taxonomy table, count table, and sample data table
tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & phyla
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_phy: includes taxonomic information for each representative (ASV) sequence
phylum <- tax_tab$Phylum #pulling out phylum column from taxonomy table
tax_tab_phy <- tibble(phylum) #making phyla into a tibble containing phylum for each sequence
tax_tab_phy <- as.matrix(tax_tab_phy) #make tibble into matrix
row.names(tax_tab_phy) <- tax_tab$Sequence #make sequence column the row names
tax_tab_phy[is.na(tax_tab_phy)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_phy: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_phy)
count_tab_phy <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_phy) <- count_tab$...1 #make sequences the row names (ignore warning message)

#sample_info_tab_phy: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab_phy <- sample_info_tab
sample_info_tab_phy <- sample_info_tab_phy[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_phy) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object 
ASV_physeq <- phyloseq(otu_table(count_tab_phy, taxa_are_rows = TRUE), tax_table(tax_tab_phy), sample_data)
ASV_physeq <- prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'allphy_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe 
df_ASV_physeq <- ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe
head(df_ASV_physeq)

#----------end of the code from GoM phyla distribution script--------

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-cnid-dist")

#import files
genbank_seqs <- read_csv("28S-Hyd-Scy-Cub-Staur-Genbank-barcode-nodecimal.csv") #file of sequences identified as HSCS in Genbank and their accession ID
accid_taxid <- read.delim("accid_taxid.txt") #file of accession IDs and tax IDs
unique_taxid <- read.delim("unique_taxid.txt") #file of unique tax IDs
blast_tax <- read.delim("blastn_hits_taxonomy.txt") #file of tax ID and BLAST taxonomy

#create table of HSCS sequences and their taxonomy by combining accesssion ID, sequence, and taxonomy
accid_taxid_unique <- accid_taxid %>% distinct(TaxID, .keep_all = TRUE) #select unique TaxIDs that are the first occurrence in df along w/ accession ID
hscs_classifier <- left_join(genbank_seqs, accid_taxid_unique, by = "Name") #join unique TaxIDs with Genbank sequences by their accession ID
blast_tax_unique <- blast_tax %>% distinct(TaxID, .keep_all=TRUE) #select unique TaxIDs that are the first occurrence in df along w/ tax classification
hscs_classifier <- right_join(hscs_classifier, blast_tax_unique, by = "TaxID") #join df of sequences w/ classification
unique(hscs_classifier$Species) #how many species were identified?

#save table to make assignTaxa and assignSpecies classifiers in CLI
write.csv(hscs_classifier, "~/Desktop/AW RStudio/results/gomx-cnid-dist/HSCS_classifier.csv")

#select only the cnidarians
cnid.table <- df_ASV_physeq %>% 
  subset(phylum == "Cnidaria") %>%
  filter(Abundance >0) #only the samples that contained cnidarians

#create vector of only unique sequences from cnid.table so we can compare them to FASTA files
seqs <- unique(cnid.table$OTU)

#assignTaxonomy using FASTA file as reference
taxa <- assignTaxonomy(seqs, "28S-Hyd-Scy-Cub-Staur_assignTaxonomy.fasta", multi=TRUE, minBoot = 80,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) #classifying to genus level with assignTaxonomy
unname(taxa)
unique(taxa [, 3]) #see unique classes identified

#addSpecies using FASTA file as reference
genus.species <- addSpecies(taxa, "28S-Hyd-Scy-Cub-Staur_assignSpecies.fasta", allowMultiple = TRUE)  #finding 100% matches to our reference database of Gulf of Mexico ctenos with assignSpecies
unique(genus.species [, 6]) #see unique genera identified 
unique(genus.species [, 7]) #see unique species identified

#convert from vector to dataframe
HSCS.species.df <- as.data.frame(genus.species)
HSCS.species.df$seq <- row.names(HSCS.species.df) #make columun of sequences from the row names
HSCS.species.df$Genus <- ifelse(HSCS.species.df$Genus == "", NA, HSCS.species.df$Genus)

#add abundance of each unique ASV to dataframe
length(unique(HSCS.species.df[["seq"]])) #find how many unique ASVs are in table

#make table of abundance counts for each unique ctenophore ASV
counts_per_ASV <- cnid.table %>% 
  group_by(OTU) %>%
  summarize(totalcount = sum(Abundance)) 
names(counts_per_ASV)[names(counts_per_ASV) == "OTU"] <- "seq" #rename this column seq so we can combine it with HSCS.species.df

#combine abundance counts table of ctenophores w/ classifier table
HSCS.species.df <- left_join(HSCS.species.df, counts_per_ASV, by = "seq")

#save as a table
setwd("~/Desktop/AW RStudio/results/gomx-cnid-dist")
write.table(HSCS.species.df, file = 'HSCS_species_classifier_results.tsv', sep = "\t", row.names = FALSE, quote=FALSE) #writing the ASV counts table with the taxonomic classifications of each cteno ASV

#create table from cnid.table with only depth, method & abundance 
cnid.depth <- cnid.table %>%
  select(OTU, Abundance, depth_bin, CTD.ROV)

#join with table of HSCS species taxonomy
cnid.species.depth <- right_join(HSCS.species.df, cnid.depth, by = c("seq" = "OTU")) 
cnid.species.depth <- cnid.species.depth %>% #take out all samples that were blanks or controls
  filter(CTD.ROV != "Negative" & CTD.ROV != "NTC")

#figure out which depth ranges were sampled at so we can make depth bins for plots
unique(cnid.species.depth$depth_bin)

#plot by genus
ggplot(cnid.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Genus)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Genera in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by family
ggplot(cnid.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") + #facet by sampling method
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Families in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by order
ggplot(cnid.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Order)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Orders in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 







