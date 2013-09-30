library(ProjectTemplate)
load.project()
set.seed(42)

# Cellular subpopulations of interest
subpopulations <- c("lymph", "CD3", "CD4", "CD25+CD127-", "CCR4+CD45RO+",
                    "CCR4-CD45RO+", "CCR4+CD45RO-", "CCR4-CD45RO-",
                    "CCR4+HLADR+", "CCR4-HLADR+", "CCR4+HLADR-", "CCR4-HLADR-",
                    "Memory", "Naive", "Total", "Activated")

gs_path <- "/loc/no-backup/ramey/Cytotrol/gating-sets/gs-treg"

# Loads the archived gatingSet object
gs_treg <- load_gs(gs_path)

# Creates the gating-template object from a CSV file
gt_csv <- "gating-templates/gt-treg.csv"
gating_template <- gatingTemplate(gt_csv)

# Applies OpenCyto to GatingSet
gating(gating_template, gs_treg, mc.cores = 8, parallel_type = "multicore")

# NOTE: The CCR4 gate for BSMS is poor, due to lack of separation.
# Several attempts were made to improve the gates for BSMS, but these changes
# degraded gates for other centers.
plotGate(gs_treg[5:8], c("CCR4+CD45RO+", "CCR4+CD45RO-", "CCR4-CD45RO+", "CCR4-CD45RO-"),
         xbin = 32, margin = FALSE)


# Hides intermediate helper gates
setNode(gs_treg, "boundary", FALSE)
setNode(gs_treg, "CD25CD127_transitional", FALSE)
setNode(gs_treg, "CD25", FALSE)
setNode(gs_treg, "HLADR+", FALSE)
setNode(gs_treg, "HLADR-", FALSE)
setNode(gs_treg, "CD45RO+", FALSE)
setNode(gs_treg, "CD45RO-", FALSE)
setNode(gs_treg, "CCR4+", FALSE)
setNode(gs_treg, "CCR4-", FALSE)

# Archives the GatingSet
save_gs(gs_treg, path = gs_path, overwrite = TRUE)

# Extracts the proportions of the subpopulations of interest and saves them to a CSV file
# For each file, we also add the given pData.
pop_stats <- reshape2::melt(getPopStats(gs_treg))
colnames(pop_stats) <- c("Population", "Filename", "Proportion")
pop_stats$Population <- as.character(pop_stats$Population)
pop_stats$Filename <- as.character(pop_stats$Filename)

# Splits the strings to extract the populations and their parents
population_split <- strsplit(pop_stats$Population, "/")
parent_pop <- sapply(population_split, function(x) head(tail(x, n = 2), n = 1))
parent_pop <- replace(parent_pop, parent_pop == "root", NA)
parent_pop <- replace(parent_pop, parent_pop == "", "root")
pop_stats$Population <- sapply(population_split, tail, n = 1)
pop_stats$Parent <- parent_pop

# Grabs the counts for each population and its parent
# Then appends the counts to the population statistics
pop_counts <- reshape2::melt(getPopStats(gs_treg, stat = "count"))
colnames(pop_counts) <- c("Population", "Filename", "Count")
pop_counts$Population <- as.character(pop_counts$Population)
pop_counts$Filename <- as.character(pop_counts$Filename)
pop_counts$Population <- sapply(strsplit(pop_counts$Population, "/"), tail, n = 1)
pop_stats <- merge(pop_stats, pop_counts)
pop_stats <- merge(pop_stats, pop_counts, by.x = c("Filename", "Parent"),
                   by.y = c("Filename", "Population"))
colnames(pop_stats)[5:6] <- c("Count", "Count_Parent")

# Subsets and polishes the population statistics
pop_stats <- subset(pop_stats, Population %in% subpopulations)
pop_stats <- merge(pop_stats, pData(gs_treg), by.x = "Filename", by.y = "name", sort = FALSE)
pop_stats <- subset(pop_stats, select = c(Population, Filename, Center, Sample, Replicate, Proportion, Count, Parent, Count_Parent))

pop_stats <- within(pop_stats, Population[Population == "lymph"] <- "Lymphocytes")
write.csv(pop_stats, file = "popstats/popstats-treg.csv", row.names = FALSE)
