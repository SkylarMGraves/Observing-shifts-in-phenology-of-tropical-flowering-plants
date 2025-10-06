#install.packages("fields")
library(rentrez)

species_list <- c("Ceiba erianthos", "Ceiba jasminodora", "Ceiba schottii", 
                  "Ceiba trischistandra", "Porcelia ponderosa", 
                  "Bougainvillea modesta", "Bougainvillea stipitata",  
                  "Barnebya harleyi", "Mimosa acutistipula",  
                  "Terminalia fagifolia", "Dioscorea bulbifera", 
                  "Combretum acutum", "Dracaena phrynioides", "Stylosanthes erecta", 
                  "Terminalia laxiflora", "Tetrorchidium didymostemon", 
                  "Vangueriella nigerica", "Diospyros lotus")

# Create a function to search and fetch GenBank accession numbers
fetch_ITS_accessions <- function(species) {
  search_term <- paste0(species, " AND internal transcribed spacer[Title]")
  result <- entrez_search(db = "nucleotide", term = search_term, retmax = 5)
  result$ids
}

# Run on species list
ITS_ids <- lapply(species_list, fetch_ITS_accessions)
names(ITS_ids) <- species_list


# Load required library
library(ape)

# Flatten accession numbers and keep track of species
accession_list <- unlist(ITS_ids)
species_lookup <- rep(names(ITS_ids), sapply(ITS_ids, length))
names(species_lookup) <- accession_list


# Download sequences from GenBank
ITS_seqs <- read.GenBank(accession_list, species.names = TRUE)

# Optional: Rename sequences to species name for clarity
names(ITS_seqs) <- gsub(" ", "_", species_lookup[names(ITS_seqs)])


# Save the sequences to a FASTA file
write.dna(ITS_seqs, file = "ITS_sequences.fasta", format = "fasta", nbcol = -1, colsep = "")
########################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")


#######################################################################################
library(msa)
library(Biostrings)

# Load sequences
seqs <- readDNAStringSet("ITS_sequences.fasta")

# Align
alignment <- msa(seqs, method = "ClustalW")  # or "Muscle", "ClustalOmega"



# Convert alignment to DNAbin (compatible with ape)
aligned_dnabin <- msaConvert(alignment, type = "ape::DNAbin")

# Save as FASTA
write.dna(aligned_dnabin, file = "aligned_ITS.fasta", format = "fasta", nbcol = -1, colsep = "")


#####################################################################################################


library(phangorn)
library(ape)

# Read alignment
aligned_seqs <- read.dna("aligned_ITS.fasta", format = "fasta")

# Convert to phyDat object (required by phangorn)
phy_data <- phyDat(aligned_seqs, type = "DNA")

# Compute distance matrix
dm <- dist.ml(phy_data)

# Build NJ tree
nj_tree <- NJ(dm)

# Optional: Optimize initial tree
nj_tree <- ladderize(nj_tree)
plot(nj_tree, main = "NJ Tree (starting tree for ML)", cex = 0.6)



# Fit a maximum likelihood tree using GTR model
fit_start <- pml(nj_tree, data = phy_data)

# Optimize model: edge lengths, substitution model, etc.
fit_ml <- optim.pml(fit_start, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 1))
################################################

any(duplicated(nj_tree$tip.label))

# Add suffixes to duplicate tip labels to ensure uniqueness
nj_tree$tip.label <- make.unique(nj_tree$tip.label)

# Also update the alignment data to match new names
names(phy_data) <- nj_tree$tip.label


fit_start <- pml(nj_tree, data = phy_data)

fit_ml <- optim.pml(fit_start, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 1))




set.seed(123)
bs <- bootstrap.pml(fit_ml, bs = 100, optNni = TRUE, multicore = FALSE, control = pml.control(trace = 0))

# Plot with bootstrap support
plotBS(fit_ml$tree, bs, p = 50, type = "phylogram", main = "ML Tree with Bootstrap", cex = 0.6)


plotBS(fit_ml$tree, bs, type = "phylogram", p = 0, main = "ML Tree with Bootstrap Support", cex = 0.6)



# Get original names from your FASTA (assuming they matched species names)
original_names <- names(read.dna("ITS_sequences.fasta", format = "fasta"))
cleaned_names <- make.unique(original_names)

# Map cleaned names back to original names (remove .1, .2 etc.)
name_map <- setNames(original_names, cleaned_names)

# Replace labels in the final tree
clean_labels <- name_map[fit_ml$tree$tip.label]
fit_ml$tree$tip.label <- clean_labels

# Also reapply for plotting bootstrap tree if using plotBS return value
# (e.g., if you did: tree_with_bs <- plotBS(...))
# tree_with_bs$tip.label <- clean_labels


plotBS(fit_ml$tree, bs, type = "phylogram", p = 0.7, main = "ML Tree with Bootstrap Support", cex = 0.6)




head(original_names, 10)

# Step 1: Original names like "Ceiba_erianthos", repeated
original_names <- names(read.dna("ITS_sequences.fasta", format = "fasta"))

# Step 2: Made unique for ML inference
cleaned_names <- make.unique(original_names)

# Step 3: Mapping from unique → species (e.g. "Ceiba_erianthos.1" → "Ceiba_erianthos")
name_map <- setNames(original_names, cleaned_names)

# Step 4: Apply clean names back to tree
fit_ml$tree$tip.label <- name_map[fit_ml$tree$tip.label]

table(fit_ml$tree$tip.label)


plotBS(fit_ml$tree, bs, type = "phylogram", p = 0,
       main = "ML Tree with Bootstrap Support", cex = 0.6)







# Add unique suffixes for plotting only
tip_labels_unique <- make.unique(fit_ml$tree$tip.label)

# Save clean → unique map for a legend if needed
label_map <- setNames(tip_labels_unique, fit_ml$tree$tip.label)

# Apply to tree
fit_ml$tree$tip.label <- tip_labels_unique

# Plot with support values
plotBS(fit_ml$tree, bs, type = "phylogram", p = 0.7,
       main = "ML Tree with Bootstrap Support", cex = 0.6)






#########################################################################

# Your original data
species <- c("Porcelia ponderosa", "Peltogyne pauciflora", "Peltogyne recifensis", "Ceiba erianthos", 
             "Ceiba jasminodora", "Ceiba schottii", "Ceiba trischistandra", "Bougainvillea modesta", 
             "Bougainvillea stipitata", "Aeschynomene martii", "Terminalia fagifolia", 
             "Caesalpinia pluviosa", "Mimosa acutistipula", "Barnebya harleyi", 
             "Pseudobombax parvifolium", "Dioscorea bulbifera", "Aeschynomene indica", 
             "Annona glauca", "Alafia barteri", "Landolphia micrantha", "Dracaena phrynioides", 
             "Combretum acutum", "Terminalia laxiflora", "Tetrorchidium didymostemon", 
             "Crotalaria mortonii", "Stylosanthes erecta", "Marantochloa leucantha", 
             "Triclisia subcordata", "Tricalysia pallens", "Vangueriella nigerica", 
             "Diospyros lotus", "Rudgea crassipetiolata", "Nymphoides aurantiaca")

DaysPerD <- c(2.534, 0.5445, 14.118, 0.8557, 2.698, 3.1352, 4.280, 1.8775, 0.6746, 1.2017,
              0.9269, 0.07914, 1.1150, 5.847, 1.7068, 0.17379, 0.7917, 0.1867, 0.7174,
              1.5413, 1.7685, 3.3746, 2.537, 10.65, 4.0833, 3.5666, 1.2530, 50.619,
              9.379, 5.273, 0.8754, 0.03695, 1.6434)

# Clean names (replace spaces with underscores)
names(DaysPerD) <- gsub(" ", "_", species)


library(dplyr)

# Extract tree tip labels (unique with .1, .2, etc)
tip_labels <- fit_ml$tree$tip.label

# Map each tip to its species base name (strip .1, .2 suffix)
species_base <- sub("\\.\\d+$", "", tip_labels)

# Create a data frame linking tip to species
tip_df <- data.frame(tip_label = tip_labels, species = species_base, stringsAsFactors = FALSE)

# Join DaysPerD values to each tip's species
tip_df$trait_value <- DaysPerD[tip_df$species]

# Aggregate trait values by tip label:
# For duplicates (e.g. Ceiba_erianthos.1, .2) we keep the same species trait value.
# So each tip gets the species trait.

# Vector for trait values matching tree tip labels exactly
trait_vec <- tip_df$trait_value
names(trait_vec) <- tip_df$tip_label




library(ape)
library(phytools)

# Keep tips with non-NA trait values only
keep_tips <- names(trait_vec)[!is.na(trait_vec)]

# Prune tree to keep only those tips
pruned_tree <- drop.tip(fit_ml$tree, setdiff(fit_ml$tree$tip.label, keep_tips))

# Prune trait vector to match pruned tree tips
trait_vec_pruned <- trait_vec[pruned_tree$tip.label]



# Make sure your tree has branch lengths
if (is.null(pruned_tree$edge.length)) {
  pruned_tree <- compute.brlen(pruned_tree, method = "Grafen")
}

# Calculate phylogenetic signal
blomberg_result <- phylosig(pruned_tree, trait_vec_pruned, method = "K", test = TRUE)

print(blomberg_result)




# Find branches with zero or near-zero length
zero_branches <- which(pruned_tree$edge.length < 1e-8)

# If any zero branches, add a tiny length to all branches to avoid zero-length issues
if(length(zero_branches) > 0){
  pruned_tree$edge.length <- pruned_tree$edge.length + 1e-6
}

blomberg_result <- phylosig(pruned_tree, trait_vec_pruned, method = "K", test = TRUE)
print(blomberg_result)






# Color tips by trait value
library(viridis)
tip_colors <- viridis(length(trait_vec_pruned))[rank(trait_vec_pruned)]

plot(pruned_tree, tip.color = tip_colors, cex = 0.8, main = "Phylogeny with Flowering Time")
tiplabels(pch = 21, bg = tip_colors, cex = 1)






library(fields)

# Use the original tree plot space
plot(pruned_tree, tip.color = tip_colors, cex = 0.8, main = "Phylogeny with Flowering Time")
tiplabels(pch = 21, bg = tip_colors, cex = 1)

# Add color key (legend)
image.plot(legend.only = TRUE, zlim = range(trait_vec_pruned), col = viridis(100),
           legend.args = list(text = 'Flowering Time', side = 4, font = 2, line = 2.5, cex=0.8))

