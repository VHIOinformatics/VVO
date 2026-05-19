library(DeeDeeExperiment)
library(iSEE)
library(readxl)
library(S4Vectors)

dde <- readRDS("dde.rds")

rd <- as.data.frame(rowData(dde))

# Comprovem que les columnes existeixen
colnames(rd)

# Filtrem files amb NAs
rd_net <- rd[!is.na(rd$Start) & !is.na(rd$End) & !is.na(rd$Chr), ]

gr <- GenomicRanges::GRanges(
  seqnames = paste0("chr", as.character(rd_net$Chr)),
  ranges   = IRanges::IRanges(
    start = as.integer(rd_net$Start),
    end   = as.integer(rd_net$End)
  ),
  strand = ifelse(rd_net$Strand == 1, "+", "-")
)
names(gr) <- rownames(rd_net)

# Conservem el rowData original als mcols del GRanges
mcols(gr) <- rd_net

# Filtrem el dde i assignem el GRanges
dde_filtrat <- dde[rownames(rd_net), ] # eliminem gens sense coordenades
rowRanges(dde_filtrat) <- gr

saveRDS(dde_filtrat, file = "dde_ranges.rds")

