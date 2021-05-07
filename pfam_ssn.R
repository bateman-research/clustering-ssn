# Compare domain sequences in a subset of related Pfam families (e.g. same clan)
# Resolve domain overlaps, count same protein and tandem and species distribution
# Plot a sequence similarity network (SSN) colored by Pfam membership
# Input Pfam families in separate folders with their ALIGN, DESC, HMM and score files
# 
# Aleix Lafita - October 2019

suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(igraph))

theme_set(theme_bw() + theme(panel.grid = element_blank()))

set.seed(0)

blast_cols = gsub(" ", "", strsplit("query acc.ver, subject acc.ver, p.identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score", ",")[[1]])


###################### Argparse #############################

pfams = "examples/PF08428,examples/PF18938,examples/PF18957"
prefix = "examples/pfam_clan"
ecod = "examples/ecod_domains.fa"
bit_thr = 30
sample = 200

# Consider species distribution too
dosp = F

ALIGN_names = c("domain_id", "alignment")
score_names = c("bits", "domain_id", "range", "evalue")

# create parser object
parser = ArgumentParser(description = 'Compare domain sequences in a subset of related Pfam families (e.g. same clan)')

# specify our desired options 
parser$add_argument("-f", "--pfams", default=pfams,
                    help="Pfam IDs semicolon separated [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")
parser$add_argument("-b", "--bitscore", default=bit_thr,
                    help="Threshold for BLAST significant hits in bitscore [default \"%(default)s\"]")
parser$add_argument("-e", "--ecod", default=ecod,
                    help="Path to the ECOD fasta sequences file [default \"%(default)s\"]")
parser$add_argument("-s", "--sample", default=sample,
                    help="Number of domains to subsample from each family [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

pfams = args$pfams
prefix = args$prefix
bit_thr = as.integer(args$bitscore)
ecod = args$ecod
sample = as.integer(args$sample)

######################### Data Parsing ##########################################
message("# Parsing data from Pfam families...")

pfam.ids = unlist(strsplit(pfams, ","))

# Parse ALIGN, scores and Pfam names
data = NULL
data.ecod = NULL

for (p in pfam.ids){
  
  ALIGN = read.csv(
    paste(p, "ALIGN", sep="/"),
    sep = "",
    header = F,
    stringsAsFactors = F,
    col.names = ALIGN_names
  )
  
  scores = read.csv(
    paste(p, "scores", sep="/"),
    sep = "",
    comment.char = "#",
    header = F,
    stringsAsFactors = F,
    col.names = score_names
  )

  #scores = scores[1,]
  data.pfam = merge(ALIGN, scores, all.x = T) %>%
    mutate(
      pfam = p,
      uniprot = gsub("\\..*", "", domain_id),
      range = gsub(".*/", "", domain_id),
      start = gsub("-.*", "", range) %>% strtoi(),
      end = gsub(".*-", "", range) %>% strtoi()
    )
  
  if (dosp) {
    species = read.csv(
      paste(p, "ALIGN.tab", sep="/"),
      sep = "\t",
      comment.char = "#",
      header = T,
      stringsAsFactors = F
    ) %>% select(
      Entry,
      Taxonomic.lineage..SPECIES.,
      Taxonomic.lineage..PHYLUM.
    )
    
    data.pfam = merge(data.pfam, species, by.x="uniprot", by.y="Entry", all.x = T)
  }
  
  # Search the family models in ECOD domains
  system(sprintf(
    "hmmsearch -o %s/ECOD_hmmer.log --domtblout %s/ECOD.domtblout --domT %i %s/HMM %s",
    p, p, bit_thr, p, ecod
  ))
  #system(sprintf("esl-reformat pfam %s/ECOD.sto > %s/ECOD.pfam", p, p))
  system(sprintf(
    "grep -v '#' %s/ECOD.domtblout | awk '{print $1}' | sed '/^$/d'| sort | uniq > %s/ECOD.ids",
    p, p
  ))
  system(sprintf("grep -A1 -f %s/ECOD.ids %s | sed '/^--$/d' > %s/ECOD.fa", p, ecod, p))
  
  if (is.null(data)) {
    data = data.pfam
    system(sprintf("cat %s/ECOD.fa > %s_ecod.fa", p, prefix))
  } else {
    data = rbind(data, data.pfam)
    system(sprintf("cat %s/ECOD.fa >> %s_ecod.fa", p, prefix))
  }
}

####################### Overlaps ############################################
message("# Resolving domain hit overlaps...")

data.overlap = data[order(data$uniprot, data$start, data$bits),]
data.overlap$id = 1:nrow(data.overlap)

data.overlap2 = data.overlap
data.overlap2$id = c(nrow(data.overlap2),1:(nrow(data.overlap2)-1))

data.overlapall = merge(data.overlap, data.overlap2, by = "id")

data.overlapall = data.overlapall %>% mutate(
  overlap = ifelse( #check same uniprot ID
    uniprot.x == uniprot.y,
    ifelse( #check end smaller than start next
      end.x - start.y > 10, 
      ifelse( #check bits lower
        bits.x < bits.y, id, id+1
      ), NA), NA
  )
)

overlaps = data.overlapall$overlap
overlaps = unique(overlaps[!is.na(overlaps)])

data.nonoverlap = data.overlap[!is.element(data.overlap$id, overlaps),]
data.overlapping = data.overlapall[!is.na(data.overlapall$overlap),]

data.overlapping.dir = data.overlapping %>%
  group_by(pfam.x, pfam.y) %>%
  summarise(count = length(uniprot.x))

#print(data.overlapping.dir)

x = data.overlapping
y = data.nonoverlap
venn = list()
for (pf in pfam.ids){
  venn[[pf]] = c(x[x$pfam.x == pf,]$id,
                 x[x$pfam.y == pf,]$id,
                 y[y$pfam == pf,]$id) %>% unique()
}

d = euler(venn)

pdf(sprintf("%s_overlap_venn.pdf", args$prefix), 6, 4)
plot(d, #fills = c("green", "red", "blue", "orange"), alpha = 0.6,
     quantities = list(fontsize = 12))
log = dev.off()

data.domains= data.nonoverlap %>% select(domain_id, pfam)

write.table(
  data.domains,
  sprintf("%s_domains.tsv", args$prefix),
  sep = "\t",
  row.names = F,
  quote = F
)

x = data.nonoverlap
venn = list()
for (pf in pfam.ids){
  venn[[pf]] = x[x$pfam == pf,]$uniprot %>% unique()
}

d = euler(venn)

pdf(sprintf("%s_protein_venn.pdf", args$prefix), 6, 4)
plot(d, quantities = list(fontsize = 12))
log = dev.off()


####################### Tandem ############################################
message("# Analysis of tandem and same protein occurances...")

data.sameprotein = data.nonoverlap %>%
  group_by(uniprot, pfam) %>%
  summarize(num = length(range)) %>%
  group_by(uniprot) %>%
  summarize(ribs = paste(pfam, collapse = ","),
            sum = sum(num),
            num = paste(num, collapse = ",")
            )

data.tandem1 = data.nonoverlap
data.tandem1$id = 1:nrow(data.tandem1)

data.tandem2 = data.tandem1
data.tandem2$id = c(nrow(data.tandem2),1:(nrow(data.tandem2)-1))

data.tandemall = merge(data.tandem1, data.tandem2, by = "id")

data.tandem = data.tandemall %>%
  mutate(distance = ifelse(
    uniprot.x == uniprot.y,
    start.y - end.x,
    #max(start.y - end.x, end.y - start.x),
    1000)) %>%
  filter(distance < 30)

x = data.tandem
venn = list()
for (pf in pfam.ids){
  venn[[pf]] = c(x[x$pfam.x == pf,]$id, x[x$pfam.y == pf,]$id) %>% unique()
}

d = euler(venn)

pdf(sprintf("%s_tandem_venn.pdf", args$prefix), 6, 4)
plot(d, quantities = list(fontsize = 12))
log = dev.off()

data.tandem = data.tandem %>%
  mutate(
    tandem = paste(pfam.x, pfam.y)
  )

p = ggplot(data.tandem %>% filter(pfam.x == pfam.y), aes(x = distance)) + 
  geom_histogram(
    alpha = 0.3,
    #aes(color = pfam),
    color = "black",
    boundary = 0, 
    binwidth = 1
  ) + 
  facet_wrap( ~ tandem, ncol = 3) +
  xlab("Linker length")

pdf(sprintf("%s_tandem_linker.pdf", args$prefix), 8, 3)
plot(p)
log = dev.off()

####################### Length ############################################
message("# Looking at the domain length distribution...")

data.length = data.nonoverlap %>%
  mutate(
    hmmalign = gsub("\\.", "", gsub("[a-z]", "", alignment)),
    sequence = toupper(gsub("[.-]", "", alignment)),
    reglen = end - start + 1,
    seqlen = nchar(sequence),
    alnlen = nchar(gsub("-", "", hmmalign)),
    pfam = factor(pfam, levels = pfam.ids),
    gaps = str_count(hmmalign, "-"),
    gapsN = nchar(hmmalign) - nchar(gsub("^-+", "", hmmalign)),
    gapsC = nchar(hmmalign) - nchar(gsub("-+$", "", hmmalign)),
    gapsI = gaps - gapsN - gapsC,
    ins = str_count(alignment, "[a-z]"),
    insN = nchar(gsub("\\.", "", alignment)) - nchar(gsub("^[a-z]+", "", gsub("\\.", "", alignment))),
    insC = nchar(gsub("\\.", "", alignment)) - nchar(gsub("[a-z]+$", "", gsub("\\.", "", alignment))),
    insI = ins - insN - insC,
    corrlen = seqlen + gapsN - insN + gapsC - insC
  )


p = ggplot(data.length, aes(x = corrlen, fill = pfam)) + 
  geom_line(
    aes(color = pfam),
    #stat = "density"
    stat = "bin", boundary = 0, binwidth = 1
  ) + 
  geom_histogram(position = "identity", boundary = 0, binwidth = 1, alpha = 0.3) + 
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  ) + xlab("Domain length")

pdf(sprintf("%s_length.pdf", args$prefix), 5, 5)
plot(p)
log = dev.off()


####################### Subsample #########################################
message("# Subsample domain sequences for BLAST similarity network...")

# Subsample the families so they are in more or less the same proportion
set.seed(0) # reproducible

data.domains.subsample = data.length %>% 
  group_by(pfam) %>% 
  sample_n(min(length(domain_id), sample))

# Write table of subsampled IDs and Pfam annotation
write.table(
  data.domains.subsample %>% select(domain_id, pfam),
  sprintf("%s_domains_subsample.tsv", args$prefix),
  sep = "\t",
  row.names = F,
  quote = F
)

# Write FASTA file of subsampled sequences
write.fasta(
  as.list(data.domains.subsample$sequence), 
  data.domains.subsample$domain_id,
  sprintf("%s_domains_subsample.fa", prefix)
)

system(sprintf("cat %s_domains_subsample.fa %s_ecod.fa > %s_blast.fa", prefix, prefix, prefix))


######################## SSN BLAST ###########################################

message("# Running BLASTp on sequences...")
system(sprintf("makeblastdb -in %s_blast.fa -dbtype prot", prefix))
system(sprintf(
  "blastp -db %s_blast.fa -query %s_blast.fa -out %s_blast.tsv -evalue 0.01 -outfmt 7 -max_target_seqs 10000", 
  prefix, prefix, prefix
))

# Parse the BLAST results
blast = read.table(
  sprintf("%s_blast.tsv", prefix),
  sep = "\t",
  header = F,
  stringsAsFactors = F,
  comment.char = "#",
  col.names = blast_cols
)

message("# BLAST results table parsed successfully...")

# Remove duplicates and cleanup - speedup following steps
blast.all = blast %>%
  filter(
    bitscore > bit_thr,
    subjectacc.ver != queryacc.ver
  ) %>%
  rowwise() %>%
  mutate(id = paste(min(queryacc.ver, subjectacc.ver), max(queryacc.ver, subjectacc.ver))) %>%
  ungroup() %>%
  group_by(id) %>%
  slice(which.max(bitscore)) %>%
  ungroup()

# Include structures from ECOD - if file is empty seqinr returns error
data.graph = data.domains.subsample
try({
  ecod.fa = read.fasta(sprintf("%s_ecod.fa", prefix))
  ecod.ids = data.frame(domain_id=unique(names(ecod.fa)), pfam = "1ECOD")
  data.graph = merge(data.domains.subsample, ecod.ids, all = T)
})

# Construct a graph with all edges
graph = graph_from_data_frame(
  blast.all,
  vertices = data.graph,
  directed = F
)

# Delete nodes with low degree - extra cleaning stage
v = degree(graph)
graph = delete_vertices(graph, names(v[v < 2]))

# For very large networks save it to a PNG
if (vcount(graph) < 2000) {
  pdf(sprintf("%s_SSN.pdf", prefix), 100, 100)
} else {
  png(sprintf("%s_SSN.png", prefix), 200, 200, "cm", res = 100)
}
plot(
  graph,
  layout=layout_with_fr(graph, grid = "nogrid"),
  #layout=layout_with_kk,
  vertex.size=ifelse(vertex_attr(graph, "pfam") == "1ECOD", 3, 1),
  vertex.label=vertex_attr(graph, "domain_id"),
  vertex.label.color = "black",
  #vertex.label.size = 0.1,
  #label.family = "serif",
  vertex.color=adjustcolor(as.integer(factor(vertex_attr(graph, "pfam"))), alpha.f = 0.2),
  vertex.frame.color = NA,
  #edge.width = 5 * edge_attr(graph, "bitscore") / edge_attr(graph, "alignmentlength"),
  edge.width = 0.1,
  edge.arrow.size=0.1
)
log = dev.off()

message("# Sequence similarity network SSN created successfully!")

######################## Species ###########################################

if (dosp) {
data.species = data.nonoverlap %>%
  filter(Taxonomic.lineage..SPECIES. != "") %>%
  mutate(pfam = factor(pfam, levels = pfam.ids))


p = ggplot(data, aes(
  y = Taxonomic.lineage..SPECIES., 
  x = pfam, 
  fill = Taxonomic.lineage..PHYLUM.)
  ) + 
  geom_tile() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(sprintf("%s_species.pdf", args$prefix), 10, 40)
plot(p)
log = dev.off()

data.phylum = data.species %>%
  #filter(Taxonomic.lineage..PHYLUM. == "Firmicutes") %>%
  select(pfam, Taxonomic.lineage..PHYLUM., Taxonomic.lineage..SPECIES.) %>%
  unique()

data.phylum.sp = data.phylum %>%
  group_by(pfam, Taxonomic.lineage..PHYLUM.) %>%
  summarise(sp_num = length(Taxonomic.lineage..SPECIES.))

p = ggplot(
  data.phylum.sp, 
  aes(y = Taxonomic.lineage..PHYLUM., 
      x = pfam, 
      label = sp_num)
  ) + geom_tile() +
  geom_text(color = "white") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

pdf(sprintf("%s_phylum_species.pdf", args$prefix), 5, 8)
plot(p)
log = dev.off()

venn = list()
for (pf in pfam.ids){
  venn[[pf]] = data.phylum[data.phylum$pfam == pf,]$Taxonomic.lineage..SPECIES.
}

d = euler(venn)

pdf(sprintf("%s_species_venn.pdf", args$prefix), 5, 4)
plot(d, quantities = list(fontsize = 12))
log = dev.off()
}
