library(readr)
intronfile <- "extdata/sra/sraintrons/chr22_SRA_introns.tsv.gz"
stopifnot(file.exists(intronfile))
introns <- read_tsv(intronfile, col_names = FALSE,
                    col_types = "ciiccccc")
colnames(introns) <- c("chr", "start", "end", "strand",
                       "left_motif", "right_motif", "sampleIdx", "counts")
sampleIdx.sp <- strsplit(introns$sampleIdx, ",")
count.sp <- strsplit(introns$counts, ",")
int_sampleIdx.sp <- lapply(sampleIdx.sp, as.integer)
int_count.sp <- lapply(count.sp, as.integer)

head(introns)

tab <- table(sapply(sampleIdx.sp, length))
tabCut <- c(tab[1:9], sum(tab[-(1:9)]))
names(tabCut)[10] <- "10+"
round(tabCut / length(sampleIdx.sp), 3)

library(Matrix)
sparse.ii <- rep(seq(along = sampleIdx.sp), times = sapply(sampleIdx.sp, length))
sparse.jj <- as.integer(unlist(sampleIdx.sp)) + 1L # But Abhi says non-zero index
sparse.xx <- as.integer(unlist(count.sp))
intM <- sparseMatrix(i = sparse.ii, j = sparse.jj, x = sparse.xx,
                     symmetric = FALSE)
dim(intM)
print(object.size(intM), units = "auto")

save(intM, file = "objects/intM.rda")
cpM <- crossprod(intM)
save(cpM, file = "objects/cpM.rda")
