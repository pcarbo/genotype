# This R script converts genotypes stored in a text file downloaded
# from your 23andme account (http://you.23andme.com/tools/data/download)
# into a PLINK .ped file. There is probably a combination of PLINK
# commands that could accomplish some of what this script does in R.
#
# Several variables need to be set to run this script:
#
#   geno.file   file containing 23andme genotype data.
#   out.prefix  name of .map and .ped files to generate.
#   iid         sample identifier in PLINK file.
#   fid         family identifier in PLINK file.
#   bim.file    Optionally, you can specify a PLINK .bim file to
#               match the SNPs against. The final .ped file will
#               contain genotypes for exactly the SNPs listed in this
#               .bim file; for any SNPs that are not in the 23andme
#               genotype file, the genotypes are set to missing ("0").
#

# SCRIPT PARAMETERS
# -----------------
# There are two variables you need to set to run this script: (1)
# geno.file, which contains the genotypes in a file downloaded from
# 23andme.com; (2) bim.file, the list of SNPs to keep (set to NULL to
# retain all the genotype data).
out.prefix <- "nick_carraway_23andme"
geno.file  <- "nick_carraway_23andme.txt"
bim.file   <- "../data/panel/1kg_hgdp.bim"
iid        <- "nc3"
fid        <- 0

# READ GENOTYPE DATA
# ------------------
# Any alleles that are not A, T, G or C are set to missing ("0").
cat("Reading genotypes from 23andme data file.\n")
geno        <- read.table(geno.file,header = FALSE,stringsAsFactors = FALSE)
names(geno) <- c("id","chr","pos","genotype")
geno        <- transform(geno,
                         A1 = substr(geno$genotype,1,1),
                         A2 = substr(geno$genotype,2,2),
                         stringsAsFactors = FALSE)
geno        <- cbind(geno,data.frame(dist = 0))
geno        <- geno[c("chr","id","dist","pos","A1","A2")]
rows        <- which(!(is.element(geno$A1,c("A","T","G","C")) &
                       is.element(geno$A2,c("A","T","G","C"))))
geno[rows,c("A1","A2")] <- "0"
cat("Loaded genotype data at",nrow(geno),"SNPs.\n")
rm(rows)

# OPTIONALLY, SELECT SNPs IN .bim FILE
# ------------------------------------
if (!is.null(bim.file)) {

  # Read the SNP data from a PLINK .bim file. Here, I ignore any
  # information provided about the genetic map distances.
  cat("Reading SNP data from PLINK .bim file.\n")
  map          <- read.table(bim.file,stringsAsFactors = FALSE)
  names(map)   <- c("chr","id","dist","pos","A1","A2")
  map[,"dist"] <- 0
  cat("Loaded map info for",nrow(map),"SNPs.\n")

  # Select SNPs that have the same identifer in the Ancestry file and
  # in the PLINK .bim file, and copy the genotypes at these SNPs into
  # a new data frame, "geno.new".
  geno.new           <- map
  geno.new[,"A1"]    <- "0"
  geno.new[,"A2"]    <- "0"
  rownames(geno)     <- geno$id
  rownames(geno.new) <- geno.new$id
  ids                <- intersect(geno$id,geno.new$id)
  geno.new[ids,]     <- geno[ids,]
  geno               <- geno.new
  rm(geno.new)

  # Summarize the results of these processing steps.
  cat("Aftering processing, we have genotypes at",
      sum(with(geno,A1 != "0" & A2 != "0")),"SNPs.\n")
  
  # Fix SNPs that differ in strand.
  cat("Fixing any strand differences to match .bim file.\n")
  for (i in c("A1","A2")) {
    geno[geno[[i]] == "A" & map$A1 != "A" & map$A2 != "A",i] <- "T"
    geno[geno[[i]] == "T" & map$A1 != "T" & map$A2 != "T",i] <- "A"
    geno[geno[[i]] == "G" & map$A1 != "G" & map$A2 != "G",i] <- "C"
    geno[geno[[i]] == "C" & map$A1 != "C" & map$A2 != "C",i] <- "G"
  }
}

# WRITE SNP DATA TO .map FILE
# ---------------------------
cat("Writing SNP data to PLINK .map file.\n")
write.table(geno[c("chr","id","dist","pos")],paste0(out.prefix,".map"),
            sep = " ",row.names = FALSE,col.names = FALSE,quote = FALSE)

# WRITE GENOTYPE DATA TO .ped FILE
# --------------------------------
cat("Writing genotypes to PLINK .ped file.\n")
write.table(t(as.matrix(c(fid,iid,0,0,0,-9,as.vector(t(geno[c("A1","A2")]))))),
            paste0(out.prefix,".ped"),sep = " ",row.names = FALSE,
            col.names = FALSE,quote = FALSE)
