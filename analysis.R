library(Matrix)
list.files("objects")

for(chr in c(1:22, "X", "Y", "M")) {
    cat(chr, "\n")
    load(sprintf("objects/chr%s_intM.rda", chr))
    csums <- colSums(intM)
    save(csums, file = sprintf("objects/chr%s_csums.rda", chr))
}


chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
names(chrs) <- chrs
csums = lapply(chrs, function(chr) {
    load(sprintf("objects/%s_csums.rda", chr), envir = environment())
    csums
})
total.csums <- colSums(do.call(rbind, csums[1:22]))

for(chr in chrs) {
    cat(chr, "\n")
    load(sprintf("objects/%s_intM.rda", chr))
    tmp <- sweep(intM, FUN = "/", MARGIN = 2, total.csums/10^6)



for(chr in c(1:22, "X", "Y", "M")) {
    cat(chr, "\n")
    load(sprintf("objects/chr%s_csums.rda", chr))
    csums <- colSums(intM)
    save(csums, file = sprintf("objects/chr%s_csums.rda", chr))
}




intM <- rbind(chr1_int, chr2_int, chr3_int)

sweep(M, 1, colSums(M))
