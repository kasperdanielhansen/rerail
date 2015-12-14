library(readr)
library(Matrix)


for(chr in as.character(c(7:22, "X", "Y", "M"))) {
    print(chr)
    intronfile <- sprintf("extdata/sra/sraintrons/chr%s_SRA_introns.tsv.gz", chr)
    stopifnot(file.exists(intronfile))
    introns <- read_tsv(intronfile, col_names = FALSE,
                        col_types = "ciiccccc")
    colnames(introns) <- c("chr", "start", "end", "strand",
                           "left_motif", "right_motif", "sampleIdx", "counts")
    sampleIdx.sp <- strsplit(introns$sampleIdx, ",")
    count.sp <- strsplit(introns$counts, ",")

    tab <- table(sapply(sampleIdx.sp, length))
    tabCut <- c(tab[1:9], sum(tab[-(1:9)]))
    names(tabCut)[10] <- "10+"
    print(round(tabCut / length(sampleIdx.sp), 3))

    sparse.ii <- rep(seq(along = sampleIdx.sp), times = sapply(sampleIdx.sp, length))
    sparse.jj <- as.integer(unlist(sampleIdx.sp)) + 1L 
    sparse.xx <- as.integer(unlist(count.sp))
    intM <- sparseMatrix(i = sparse.ii, j = sparse.jj, x = sparse.xx,
                         symmetric = FALSE)
    print(dim(intM))
    print(object.size(intM), units = "auto")
    save(intM, file = sprintf("objects/chr%s_intM.rda", chr))
}


## cpM <- crossprod(intM)
## save(cpM, file = "objects/cpM.rda")

## system.time({ svdM <- svd(cpM) })
## save(svdM, file = "objects/svdM.rda")

metaData <- read_tsv("extdata/sra/sraintrons/parsed_biosample_data_with_srr.tsv")
save(metaData, file = "objects/metaData.rda")
