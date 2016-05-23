# Function for reading genotype data from tab-delimited "raw DNA" text
# file downloaded from ancestry.com.
read.geno.ancestrydna <- function (file) {
  out <- read.table(file,sep = "\t",comment.char = "#",header = TRUE,
                    stringsAsFactors = FALSE)
  return(transform(out,chromosome = factor(chromosome)))
}
