INPUT_DIR <- './scRNAseq/'
OUTPUT_DIR <- './scRNAseq/'
mymeta <- read.delim(paste(INPUT_DIR, "metadata_cols.tsv", sep=""))

myday0 <- read.delim(paste(INPUT_DIR, "counts.day0.tsv", sep=""))
myday1 <- read.delim(paste(INPUT_DIR, "counts.day1.tsv", sep=""))
myday2 <- read.delim(paste(INPUT_DIR, "counts.day2.tsv", sep=""))
myday3 <- read.delim(paste(INPUT_DIR, "counts.day3.tsv", sep=""))

myday0 <- as.data.frame(t(myday0))
myday1 <- as.data.frame(t(myday1))
myday2 <- as.data.frame(t(myday2))
myday3 <- as.data.frame(t(myday3))

write.table(file = paste(INPUT_DIR, "counts.day0.t.tsv", sep=""), myday0, sep="\t", quote=FALSE, row.names=TRUE)
write.table(file = paste(INPUT_DIR, "counts.day1.t.tsv", sep=""), myday1, sep="\t", quote=FALSE, row.names=TRUE)
write.table(file = paste(INPUT_DIR, "counts.day2.t.tsv", sep=""), myday2, sep="\t", quote=FALSE, row.names=TRUE)
write.table(file = paste(INPUT_DIR, "counts.day3.t.tsv", sep=""), myday3, sep="\t", quote=FALSE, row.names=TRUE)

myday0_m <- merge(mymeta, myday0, by="donor_id")
myday1_m <- merge(mymeta, myday1, by="donor_id")
myday2_m <- merge(mymeta, myday2, by="donor_id")
myday3_m <- merge(mymeta, myday3, by="donor_id")

myday0 = melt(myday0_m, id.vars = names(myday1_m)[2:])
myday1 = melt(myday1_m, id.vars = names(myday1_m)[2:])
myday2 = melt(myday2_m, id.vars = names(myday1_m)[2:])
myday3 = melt(myday3_m, id.vars = names(myday1_m)[2:])


write.table(file = paste(OUTPUT_DIR, "counts.day0.input.tsv", sep=""), myday0, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file = paste(OUTPUT_DIR, "counts.day1.input.tsv", sep=""), myday1, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file = paste(OUTPUT_DIR, "counts.day2.input.tsv", sep=""), myday2, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file = paste(OUTPUT_DIR, "counts.day3.input.tsv", sep=""), myday3, sep="\t", quote=FALSE, row.names=FALSE)
