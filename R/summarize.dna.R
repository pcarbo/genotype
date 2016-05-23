# TO DO: Explain here what this script does.
source("functions.R")

# SCRIPT PARAMETERS
# -----------------
# Tab-delimited "raw DNA" text file downloaded from ancestry.com.
geno.file <- "ancestrydna_peter_carbonetto.txt.gz"

# READ GENOTYPE DATA
# ------------------
cat("Reading genotype data from ",geno.file,".\n",sep="")
geno <- read.geno.ancestrydna(geno.file)
n    <- nrow(geno)
cat("Read genotypes called at",n,"SNPs on",nlevels(geno$chromosome),
    "chromosomes.\n")

# SUMMARIZE GENOTYPE DATA
# -----------------------
# Count number of homozygous and heterozygous genotypes on each chromosome.
print(with(geno,
           data.frame(chr   = levels(chromosome),
                      hom   = tapply(allele1 == allele2,chromosome,sum),
                      het   = tapply(allele1 != allele2,chromosome,sum),
                      total = as.vector(table(chromosome)))),
      row.names = FALSE)
