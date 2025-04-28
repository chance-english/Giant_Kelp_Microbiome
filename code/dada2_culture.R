path <- "/Users/chanceenglish/Desktop/Lab Shiz/ARPA-e/Giant_Kelp_Microbiome/Inputdata/fastq_cultures"

fnFS <- list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE)

fnRS <- list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE)

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"

all_orients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(forward = dna, complement = complement(dna), reverse = reverse(dna),
               reverse_complement = reverseComplement(dna))
  return(sapply(orients, toString))
}

FWD.orients <- all_orients(FWD)

REV.orients <- all_orients(REV)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits>0))
}


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFS[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRS[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFS[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRS[[1]]))

sample.names <- sapply(strsplit(basename(fnFS), "_L00"), '[',1)

filt_path <- file.path(path, "filtered_culture_amplicons_path")

filtFS <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRS <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(fnFS, filtFS, fnRS, filtRS, truncLen = c(240,150), maxN = 0, 
                     maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)

errF <- learnErrors(filtFS, multithread = TRUE)

errR <- learnErrors(filtRS, multithread = TRUE)

derepFS <- derepFastq(filtFS, verbose = TRUE)
derepRS <- derepFastq(filtRS, verbose = TRUE)

names(derepFS) <- sample.names
names(derepRS) <- sample.names

dadaFS <- dada(derepFS, err = errF, multithread = TRUE)
dadaRS <- dada(derepRS, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFS, derepFS, dadaRS, derepRS, verbose = TRUE, trimOverhang = TRUE)

head(mergers[[1]])

seqtable <- makeSequenceTable(mergers)


dim(seqtable)

df_seq_table <- as.data.frame(seqtable, row.names = TRUE)

write.table(t(seqtable), file="cluture_seqtable.txt",sep = ",",row.names=FALSE, col.names=TRUE)

write_csv(df_seq_table,  "/Users/chanceenglish/Desktop/df_seq_table.csv")

seqtable_nochim <- removeBimeraDenovo(seqtable, verbose = TRUE)

taxa <- assignTaxonomy(seqtable_nochim, "/Users/chanceenglish/Desktop/Lab Shiz/ARPA-e/Giant_Kelp_Microbiome/Inputdata/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)



