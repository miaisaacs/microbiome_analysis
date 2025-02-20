# Introduction to Microbiome Analysis 

library(phyloseq)
library(microbiome)
library(ape)
library(vegan)

### if the above packages are not installed...
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages('ape')
install.packages('vegan')

### Use this dataset in package "phyloseq"
data(GlobalPatterns)

### Extract dataset

OTU_tab = t( otu_table(mb_data) )
OTU_tab[1:5, 1:20]
dim(OTU_tab) # n = 26; p = 19216
### Taxon table
Tax_tab = tax_table(mb_data)
Tax_tab[1:5,]
### Sample info
Sample_info = sample_data(mb_data)
Sample_info[1:5, ]

### Get relative abundances
scaled_OTU_tab = OTU_tab / rowSums(OTU_tab)
scaled_OTU_tab[1:5, 1:20]
rowSums(scaled_OTU_tab)

### Aggregate OTUs
Tax_tab[1:5,]
colnames(Tax_tab)
Order_mb_data = microbiome::aggregate_taxa(mb_data, 
                                           level = 'Order', 
                                           verbose = F)
Order_OTU_tab = t( otu_table(Order_mb_data) )
Order_OTU_tab[1:5, 1:10]

### Diversity calculation
### alpha diversity
alpha_div = microbiome::alpha(mb_data, 
                              index = c('shannon', 'dominance_simpson'))
alpha_div$diversity_shannon
alpha_div$dominance_simpson

### beta diversity
UN_unifrac = as.matrix( phyloseq::distance(mb_data, method = c("unifrac")) )
W_unifrac = as.matrix( phyloseq::distance(mb_data, method = c("wunifrac")) )
BC_dist = as.matrix( phyloseq::distance(mb_data, method = c("bray")) )
UN_unifrac[1:5, 1:5]
W_unifrac[1:5, 1:5]
BC_dist[1:5, 1:5]

### Cluster plots
par(mfrow = c(1,2))
### Cluster on OTUs
P_tree = phy_tree(mb_data) # tree
clustering = cluster::pam( as.dist(cophenetic(P_tree)), k = 10, diss = T)$clustering # do PAM clustering; would take a long time
load('./Clustering.Rdata') # save('clustering', file = './Clustering.Rdata')
ape::plot.phylo(P_tree, 
                show.tip.label = F, 
                direction = 'downwards', 
                edge.color = 'black', 
                edge.width = 0.5, 
                main = 'Cluster on OTUs')
ape::nodelabels(node = 1:length(P_tree$tip.label), 
                pch = 16, cex = 0.2, col = rep(rainbow(10), table(clustering)))
### Cluster on samples
HC = hclust(as.dist(BC_dist))
plot(HC, xlab = '', main = 'Cluster on Subjects')

### MDS(PCoA) plot
par(mfrow = c(1,1))
MDS_plot = phyloseq::plot_ordination(mb_data, 
                                     ordinate(mb_data, "MDS", "bray"), 
                                     color = "SampleType")
print(MDS_plot)

### Add Group 1 and 2
table(Sample_info$SampleType)
Sample_info$Group = ifelse(Sample_info$SampleType %in% c('Feces', 'Skin', 'Tongue', 'Mock'), 
                           '1', '2')
sample_data(mb_data) = Sample_info

### Test alpha-diversities between two groups
shannon_index = alpha_div$diversity_shannon
Group1_sam = which(Sample_info$Group == 1)
Group2_sam = which(Sample_info$Group == 2)
test_results = t.test(shannon_index[Group1_sam],
                      shannon_index[Group2_sam],
                      alternative = 'two.sided')
test_results
test_results$p.value

### Testing relative abundance of a taxon between two groups
temp_tab = t( otu_table(Order_mb_data) )
temp_tab = temp_tab / rowSums(temp_tab) # relative abundance
pvalue_vec = vector()
for (j in 1:ncol(temp_tab)) {
  Y = temp_tab[,j]
  ### Extract group values
  Y1 = as.vector( Y[Sample_info$Group == '1'] )
  Y2 = as.vector( Y[Sample_info$Group == '2'] )
  ### Do test
  test_results = t.test(Y1, Y2, alternative = 'two.sided')
  # test_results = wilcox.test(Y1, Y2, alternative = "two.sided")
  pvalue_vec[j] = test_results$p.value
}
# Adjust pvalues
adj_pvalue_vec_BH = p.adjust(pvalue_vec, method = 'BH')
Order_names = colnames(temp_tab)
Order_names[ which(adj_pvalue_vec_BH <= 0.1) ]

### PERMANOVA
set.seed(12345)
permanova_results = vegan::adonis(BC_dist ~ Group,
                                  data = data.frame(Sample_info),
                                  permutations = 10000)
permanova_results$aov.tab # test results
permanova_results$aov.tab['Group', 'Pr(>F)'] # pvalue


# sw2206@cumc.columbia.edu
# tw2697@cumc.columbia.edu
