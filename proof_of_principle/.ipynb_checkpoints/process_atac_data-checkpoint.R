library(cisTopic)
setwd("/data/hoan/multiomics")

counts_mel <- read.table(file = 'data/snare_seq/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv', sep = '\t', header = TRUE,row.names = 1)

# t_counts <- t(counts_mel)
# rownames(t_counts) = colnames(counts_mel)

cisTopicObject <- createcisTopicObject(counts_mel, project.name='GSE126074_CellLineMixture')

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(10:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=20, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')

data <- cisTopicObject@selected.model$document_expects

write.table(data, file = "data/snare_seq/GSE126074_CellLineMixture_SNAREseq_chromatin_topics.tsv", sep = "\t", row.names = F)