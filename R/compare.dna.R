# TO DO: Explain here what this script does.
source("functions.R")

# SCRIPT PARAMETERS
# -----------------
# Tab-delimited "raw DNA" text files downloaded from ancestry.com.
geno.file1 <- "ancestrydna_peter_carbonetto.txt.gz"
geno.file2 <- "ancestrydna_lenora_carbonetto.txt.gz"

# READ GENOTYPE DATA
# ------------------
cat("Reading genotype data from ",geno.file1,".\n",sep="")
geno1 <- read.geno.ancestrydna(geno.file1)
cat("Read genotypes called at",nrow(geno1),"SNPs on",
    nlevels(geno1$chromosome),"chromosomes.\n")

cat("Reading genotype data from ",geno.file2,".\n",sep="")
geno2 <- read.geno.ancestrydna(geno.file2)
cat("Read genotypes called at",nrow(geno2),"SNPs on",
    nlevels(geno2$chromosome),"chromosomes.\n")

# COMPARE GENOTYPES
# -----------------
# Count the number of SNPs on chromosomes 1-22 in which: (1) the
# genotype is exactly the same; (2) the genotype differs by one allele
# (letter); (3) the genotype differs in two alleles (letters).
markers <- which(is.element(geno1$chromosome,1:22))
n0 <- sum((geno1$allele1[markers] == geno2$allele1[markers] &
           geno1$allele2[markers] == geno2$allele2[markers]) |
          (geno1$allele1[markers] == geno2$allele2[markers] &
           geno1$allele2[markers] == geno2$allele1[markers]))
n2 <- sum(geno1$allele1[markers] != geno2$allele1[markers] &
          geno1$allele2[markers] != geno2$allele2[markers])
n1 <- length(markers) - n0 - n2
cat("Number of SNPs in which genotype is exactly the same:         ",n0,"\n")
cat("Number of SNPs in which genotype differs in 1 allele/letter:  ",n1,"\n")
cat("Number of SNPs in which genotype differs in 2 alleles/letters:",n2,"\n")
