# Create a all against all sequence identity matrix of domains in a protein (or few proteins)
# Input is a Pfam alignment - large alignments take forever, max 1000 sequences
# For convenience, it also saves a list of high identity tandem repeat domains (HITRDs)
#
# Aleix Lafita - October 2019

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

theme_set(theme_bw() + theme(panel.grid = element_blank()))

###################### Functions #############################

pident = function(x, y) {
  equals = 0
  gaps = 0
  for (i in seq(1, nchar(x), 1)){
    a = substr(x, i, i)
    b = substr(y, i, i)
    if (a != "-" && b != "-") {
      if (a == b) {
        equals = equals + 1
      }
    } else
      gaps = gaps + 1
  }
  return (equals / max((nchar(x) - gaps), 1))
}

###################### Argparse #############################

input = "examples/PF08428/ALIGN"
hitrd_pid = 0.7
adjacent_dist = 30

# create parser object
parser = ArgumentParser(
  description='All vs all sequence identity matrix of domains in a protein')

# specify our desired options 
parser$add_argument("-a", "--align", default=input,
                    help="Input Pfam alignment in Pfam format [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=input,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")
parser$add_argument("-i", "--pid", default=hitrd_pid,
                    help="Percentage of identity threshold for HITRDs [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

align = args$align
prefix = args$prefix
hitrd_pid = as.numeric(args$pid)

######################### Data Parsing #############################
message("# Parsing alignment files...")

# Parse the Pfam protein alignment
data = read.csv(
  align,
  sep = "",
  comment.char = "#",
  header = F,
  stringsAsFactors = F,
  col.names = c("domain_id", "alignment")
) %>% filter(domain_id != "//")

# Unpack the data into different columns
data.mut = data %>%
  mutate(
    seqid = gsub("/.*", "", domain_id),
    range = gsub(".*/", "" , domain_id),
    start = as.numeric(gsub("-.*", "", range)),
    end = as.numeric(gsub(".*-", "", range)),
    hmmalign = gsub("\\.", "", gsub("[a-z]", "", alignment)),
    sequence = toupper(gsub("[.-]", "", alignment)),
    reglen = end - start + 1,
    seqlen = nchar(sequence),
    alnlen = nchar(gsub("-", "", hmmalign))
  ) #%>% filter(is.element(seqid, c("Q5FIM8.1", "Q8E1C4.1", "Q9XDB6.1")))

# Length of the protein HMM model
hmmlen = nchar(data.mut$hmmalign[1])

########################### PID Heatmap ################################

if (nrow(data.mut) > 1000) {
  message("# Too many domains, would take too long. Reducing number to same protein.")
  data.sqrd = merge(data.mut, data.mut, by = "seqid") %>%
    mutate(
      seqid.x = seqid,
      seqid.y = seqid
    )
} else {
  data.sqrd = merge(data.mut, data.mut, by = c())
}

message("# Calculating pairwise sequence identities... (might take a while)")

domain_id.labels = data.mut$domain_id[order(data.mut$seqid, data.mut$start, decreasing = T)]
data.sqrd = data.sqrd %>%
  rowwise() %>%
  mutate(
    domain_id.x = factor(domain_id.x, domain_id.labels),
    domain_id.y = factor(domain_id.y, domain_id.labels),
    pid = pident(hmmalign.x, hmmalign.y)
  )

message("# Calculation done!")

# Plot the sequence identity confusion matrix
p = ggplot(data.sqrd, aes(x = domain_id.x, y = domain_id.y, fill = pid)) +
  geom_tile() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 270, hjust = 0),
    axis.text.y = element_text(angle = 0, hjust = 0)
  ) +
  scale_fill_continuous(low = "white", high = "black") +
  coord_fixed() +
  xlab("") + ylab("")

pdf(sprintf("%s_pidmat.pdf", prefix), 20, 20)
plot(p)
log = dev.off()

message(sprintf("# Matrix plot saved to %s_pidmat.pdf", prefix))

########################### HITRDs ################################

# Calculate the percentage of identity within adjacent domains

data.filt = data.sqrd %>%
  filter(alnlen.x > 20, alnlen.y > 20)

data.sameprot = data.filt %>%
  filter(
    seqid.x == seqid.y,
    start.y > start.x
  ) %>% mutate(
    distance = start.y - end.x,
    distance = ifelse(distance < 0, 0, distance)
  )

data.adj = data.sameprot %>%
  filter(
    start.x < start.y,
    distance < adjacent_dist
  )

data.nonadj = data.sameprot %>%
  filter(distance > adjacent_dist)

data.diffprot = data.filt %>%
  filter(seqid.x > seqid.y)

# Calculate percentages of identity of different groups
pid.sameprot = mean(data.sameprot$pid)
num.sameprot = nrow(data.sameprot)
pid.adj = mean(data.adj$pid)
num.adj = nrow(data.adj)
pid.nonadj = mean(data.nonadj$pid)
num.nonadj = nrow(data.nonadj)
pid.diffprot = mean(data.diffprot$pid)
num.diffprot = nrow(data.diffprot)

message("# Domain percentages of identity:")
message(sprintf(" - Same protein: %.2f (%i)", pid.sameprot, num.sameprot))
message(sprintf("   - Adjacent : %.2f (%i)", pid.adj, num.adj))
message(sprintf("   - Non-adjacent: %.2f (%i)", pid.nonadj, num.nonadj))
message(sprintf(" - Different protein: %.2f (%i)", pid.diffprot, num.diffprot))

# Write a file with a summary of the percentages of identity
write.table(
  data.frame(
    pid.sameprot, num.sameprot,
    pid.adj, num.adj,
    pid.nonadj, num.nonadj,
    pid.diffprot, num.diffprot
  ),
  sprintf("%s_pid-summary.tsv", prefix),
  sep = "\t",
  quote = F,
  row.names = F
)


# We need to duplicate the data frame 
data.adj2 = data.adj %>% 
  mutate(
    domain_id.z = domain_id.x,
    domain_id.x = domain_id.y,
    domain_id.y = domain_id.z
  ) %>% select(-domain_id.z)

data.adjpid = rbind(data.adj, data.adj2) %>% 
  ungroup() %>%
  mutate(
    domain_id = domain_id.x,
    seqid = seqid.x
  ) %>%
  group_by(seqid, domain_id) %>% 
  summarise(adjpid = max(pid)) %>%
  mutate(hitrd = ifelse(adjpid > hitrd_pid, 1, 0)) %>%
  ungroup()

# Plot the correlation of domain number and HITRD number
data.domnum = data.adjpid %>%
  group_by(seqid) %>%
  summarise(
    domnum = length(domain_id),
    hitrds = sum(hitrd)
  ) %>% ungroup()

# Plot the correlation of HITRDs and domain number
p = ggplot(data.domnum, aes(x = domnum, y = hitrds)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  geom_abline(alpha = 0.5) +
  coord_fixed(xlim = c(0, max(data.domnum$domnum)), ylim = c(0, max(data.domnum$domnum))) +
  xlab("Number of domains") + ylab("Number of HITRDs")

pdf(sprintf("%s_domnum.pdf", prefix), 6, 6)
plot(p)
log = dev.off()

# Now save the adjacent table into a list
write.table(
  data.adjpid %>% select(domain_id, adjpid),
  sprintf("%s_adjpid.tsv", prefix),
  quote = F,
  row.names = F,
  sep = "\t"
)

# Now save the HITRDs into a list
write.table(
  data.adjpid %>% filter(hitrd == 1) %>% select(domain_id, adjpid),
  sprintf("%s_hitrds.tsv", prefix),
  quote = F,
  row.names = F,
  sep = "\t"
)

message(sprintf("# List of HITRDs saved to %s_hitrds.tsv", prefix))

message("# Done!")


