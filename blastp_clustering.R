# Cluster a set of protein sequences using connected components of BLAST hits
# Run all against all BLAST search, build an igraph and extract connected components
# Beware of the potential file size and runtime for files with more than 1000 proteins
#
# Aleix Lafita
# created June 2019
# modified October 2019

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(data.table))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(igraph))

theme_set(theme_bw() + theme(panel.grid = element_blank()))
set.seed(0)

###################### Argparse #############################

blast_cols = gsub(" ", "", strsplit("query acc.ver, subject acc.ver, p.identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score", ",")[[1]])

seqs = "examples/sequences.fa"
prefix = "examples/sequences"
min_bit = 30
bit_thr = 100
pid = 50
cov = 50

# create parser objects
parser = ArgumentParser(description = 'Cluster a set of protein sequences using connected components of BLAST hits')

# specify our desired options 
parser$add_argument("-s", "--seqs", default=seqs,
                    help="Protein sequences in FASTA format [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")
parser$add_argument("-c", "--clustscore", default=bit_thr,
                    help="Threshold for connected components clustering in bitscore [default \"%(default)s\"]")
parser$add_argument("-t", "--hitscore", default=min_bit,
                    help="Threshold for BLAST significant hits in bitscore [default \"%(default)s\"]")
parser$add_argument("-i", "--pid", default=pid,
                    help="Minimum percentage of identity for connected components clustering [default \"%(default)s\"]")
parser$add_argument("-v", "--cov", default=cov,
                    help="Minimum percentage of identity for connected components clustering [default \"%(default)s\"]")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

seqs = args$seqs
prefix = args$prefix
bit_thr = as.numeric(args$clustscore)
min_bit = as.numeric(args$hitscore)
pid = as.numeric(args$pid)
cov = as.numeric(args$cov)

################################## BLAST #####################################

if (file.exists(sprintf("%s_blast.tsv", prefix))) {
  message(sprintf("# File '%s_blast.tsv' already exists, skipping BLASTp computation...", prefix))

} else {
  message("# Running BLASTp on sequences...")
  system(sprintf("makeblastdb -in %s -dbtype prot", seqs))
  system(sprintf(
    "blastp -db %s -query %s -out %s_blast.tsv -evalue 0.01 -outfmt 7 -max_target_seqs 10000", 
    seqs, seqs, prefix))
  system(sprintf("grep -v '#' %s_blast.tsv > %s_blast_dt.tsv", prefix, prefix))
}


# Parse the BLAST results
blast = fread( # read.table(
  sprintf("%s_blast_dt.tsv", prefix),
  sep = "\t",
  header = F,
  stringsAsFactors = F,
  #comment.char = "#",
  col.names = blast_cols
)

message("# BLAST results table parsed successfully...")

# Parse the sequences and convert to DF
seqs.fa = read.fasta(seqs)
seqs.df = data.frame(
  seqid=names(seqs.fa), 
  seq=toupper(unlist(getSequence(seqs.fa, as.string=T))),
  stringsAsFactors = F
) %>% mutate(
  seqlen = nchar(seq)
)
seqs.seqlen = seqs.df %>% select(seqid, seqlen)

# Remove duplicates and cleanup - speedup following steps
blast.all = blast %>%
  filter(
    bitscore > min_bit,
    subjectacc.ver != queryacc.ver
  ) %>%
  rowwise() %>%
  mutate(id = paste(min(queryacc.ver, subjectacc.ver), max(queryacc.ver, subjectacc.ver))) %>%
  ungroup() %>%
  group_by(id) %>%
  slice(which.max(bitscore)) %>%
  ungroup()

# Free up memory
blast = NULL

# Calculate the coverage of the alignment on query and subject
blast.cov = merge(blast.all, seqs.seqlen, by.x = "queryacc.ver", by.y = "seqid")
blast.cov = merge(blast.cov, seqs.seqlen, by.x = "subjectacc.ver", by.y = "seqid")

blast.cov = blast.cov %>%
  rowwise() %>%
  mutate(
    cov.q = abs(q.end - q.start) / seqlen.x,
    cov.s = abs(s.end - s.start) / seqlen.y,
    seqlen.min = min(seqlen.x, seqlen.y),
    cov.max = max(cov.q, cov.s) * 100,
    cov.min = min(cov.q, cov.s) * 100,
    seqlen.diff = abs(seqlen.x - seqlen.y)
  ) %>% ungroup()

# Free up more memory
blast.all = NULL

# Filter hits before the clustering step
blast.filter = blast.cov %>% 
  filter(
    bitscore > bit_thr,
    p.identity > pid,
    cov.min > cov
  )

######################### Graph clustering #################################
message("# Clustering sequences using connected components from BLAST hits...")

# Construct a graph and do the clustering with connected components
graph.filter = graph_from_data_frame(
  blast.filter,
  directed = F
)

graph.cc = components(graph.filter, mode ="strong")

# Extract the degree of nodes
node.degree = data.frame(
  id = names(degree(graph.filter)),
  degree = degree(graph.filter)
)

# Convert clustering results into a dataframe
clust = data.frame(
  id = names(graph.cc$membership),
  clust.id = graph.cc$membership,
  stringsAsFactors = F
)

# Extract cluster representatives
clust.repr = merge(node.degree, clust) %>%
  group_by(clust.id) %>%
  slice(which.max(degree)) %>%
  ungroup()

######################### Save clusters ###################################
message("# Saving sequence cluster mapping and cluster representatives...")

clust.size = clust %>%
  group_by(clust.id) %>%
  summarise(count = length(id), id = min(id)) %>%
  ungroup()

# Look at the distribution of cluster sizes
p = ggplot(clust.size, aes(x = count)) +
  geom_histogram(binwidth = 1, boundary = 0, alpha = 0.5, color = "black") +
  #geom_line(stat = "bin", binwidth = 1, boundary = 0, alpha = 0.5) +
  #scale_y_log10() +
  xlim(0, max(clust.size$count)) +
  xlab("Cluster size") +
  ylab("Number of clusters")

pdf(sprintf("%s_cluster-size.pdf", prefix), 5, 5)
print(p)
log = dev.off()

# Write a file with the results of the clustering
write.table(
  clust,
  sprintf("%s_clusters.tsv", prefix),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

# Also write the cluster representatives
write.table(
  clust.repr %>% select(id),
  sprintf("%s_cluster-repr.id", prefix),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

######################### Plot igraph ###################################
message("# Plotting BLAST hits network colored by clusters...")

blast.clust = blast.cov %>%
  filter(
    is.element(queryacc.ver, clust$id),
    is.element(subjectacc.ver, clust$id)
  )

# Construct a graph with all edges
graph = graph_from_data_frame(
  blast.clust,
  vertices = clust,
  directed = F
)

# Delete nodes with low degree - extra cleaning stage
#v = degree(graph)
#graph = delete_vertices(graph, names(v[v < 2]))

if (vcount(graph) < 2000) {
  pdf(sprintf("%s_blast-graph.pdf", prefix), 100, 100)
} else {
  png(sprintf("%s_blast-graph.png", prefix), 100, 100, "cm", res = 100)
}
plot(
  graph,
  layout=layout_with_fr(graph, grid = "nogrid"),
  #layout=layout_with_kk,
  vertex.size=1.2,
  vertex.label=as.factor(vertex_attr(graph, "clust.id")),
  vertex.label.color = "black",
  #label.family = "serif",
  vertex.color=adjustcolor(vertex_attr(graph, "clust.id"), alpha.f = 0.2),
  vertex.frame.color = NA,
  edge.arrow.size=0.1,
  edge.width = 5 * edge_attr(graph, "bitscore") / edge_attr(graph, "seqlen.min")
)
log = dev.off()

message("# BLASTp clustering is done!")
