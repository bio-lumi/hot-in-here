setwd("~/Desktop/r-summerproject/glucose-analysis")
library(qiime2R)
library(phyloseq)
library(NetCoMi)
library(igraph)
library(vegan)
library(ape)
library(dplyr)
library(viridis)


# 1. Read QIIME2 artifacts to phyloseq
gluc.physeq <- qza_to_phyloseq(
  features="glucose-clust-feature-table.qza",
  taxonomy="glucose-tax-table.qza",
  metadata = "glucose-metadata.tsv"
)

sample_data(gluc.physeq)$condition <- as.factor(sample_data(gluc.physeq)$condition)

otu_table(gluc.physeq)
tax_table(gluc.physeq)
taxa_names(gluc.physeq) <- 1:ntaxa(gluc.physeq)

# 2. Remove samples with less than 1000 reads
gluc.physeq <- prune_samples(sample_sums(gluc.physeq) > 1000, gluc.physeq)
gluc.physeq

# 3. Estimate alpha-diversity 
a.div.physeq <- estimate_richness(gluc.physeq)
write.csv(a.div.physeq, file = "~/Desktop/r-summerproject/glucose-analysis/a-div-physeq.csv")

# 4. Transform to relative abundances and plot relative abundances
gluc.physeq.family <- tax_glom(gluc.physeq, taxrank = "Family")
df.gluc.family<- psmelt(gluc.physeq.family)

family.top8 <- df.gluc.family %>%
  group_by(Family) %>%
  summarize(total.abund = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total.abund))

top.families <- head(family.top8$Family, 8)

df.gluc.family <- df.gluc.family %>%
  mutate(Family = case_when(
    Family %in% top.families ~ Family,
    Family == "uncultured" ~ "Unknown",
    TRUE ~ "Others"
  ))

df.gluc.family.merged <- df.gluc.family %>%
  group_by(condition, Family) %>%
  summarize(total.abund = sum(Abundance), .groups = "drop")

family.merged.relabund <- df.gluc.family.merged %>%
  group_by(condition) %>%
  mutate(rel.abund = total.abund / sum(total.abund)) %>%
  ungroup()

family.levels <- unique(family.merged.relabund$Family)
family.levels <- sort(family.levels[!family.levels %in% c("Others", "Unknown")])
family.levels <- c(family.levels, "Others", "Unknown")

family.merged.relabund$Family <- factor(family.merged.relabund$Family, levels = family.levels)

plot.gluc.family <- ggplot(family.merged.relabund, aes(x = condition, y = rel.abund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Temperature (°C)", y = "Relative Abundance") +
  theme_cowplot(12)

plot.gluc.family

# 5. Move phyloseq data to vegan by first converting to simpler form. In vegan, we will
# create an ordination plot to measure beta diversity. I couldn't do this in phyloseq even
# though it has a function - some technical issue that is way out of my depth to fix. 

# Convert sample data to a dataframe, and the otu table to a matrix.
str(sample_data(gluc.physeq))
str(otu_table(gluc.physeq))

otu <- otu_table(gluc.physeq)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
mat <- as(otu, "matrix")
gluc.matrix <- as.data.frame(mat)
custom_colors <- c("15" = "#f9de8a", "20" = "#ffa242", "25" = "#DC4D01", "30" =  "#792A00")

pcoa.gluc = ordinate(gluc.physeq, "PCoA", "bray") 
ordination.scores <- as.data.frame(pcoa.gluc$vectors)
gluc.ord <- data.frame(ordination.scores, sample_data(gluc.physeq))
gluc.ord.plot <- ggplot(gluc.ord, aes(x = Axis.1, y = Axis.2, color = condition)) +
  geom_point(size = 5) +
  labs(x = "PCoA.1 [47.4%]",
       y = "PCoA.2 [35.6%]",
       color = "Temperature (°C)") +
  theme_cowplot(12) +
  scale_color_manual(values = custom_colors)
gluc.ord.plot

# 6. Conduct a PERMANOVA on BC distances
gluc.bray <- phyloseq::distance(gluc.physeq, method = "bray")
gluc.meta.df <- data.frame(sample_data(gluc.physeq))
adonis2(gluc.bray ~ condition, data = gluc.meta.df)

# 7. Subset the phyloseq object according to conditions to prepare for network construction.
otu.gluc <- phyloseq::otu_table(gluc.physeq)
otu.gluc <- t(as(otu.gluc, "matrix"))
gluc.15 <- subset_samples(gluc.physeq, condition == "15")
otus.15 <- otu_table(gluc.15)
write.csv(otus.15, file = "~/Desktop/r-summerproject/glucose-analysis/cohesion/otus-15.csv")

gluc.20 <- subset_samples(gluc.physeq, condition == "20")
otus.20 <- otu_table(gluc.20)
write.csv(otus.20, file = "~/Desktop/r-summerproject/glucose-analysis/cohesion/otus-20.csv")

gluc.25 <- subset_samples(gluc.physeq, condition == "25")
otus.25 <- otu_table(gluc.25)
write.csv(otus.25, file = "~/Desktop/r-summerproject/glucose-analysis/cohesion/otus-25.csv")

gluc.30 <- subset_samples(gluc.physeq, condition == "30")
otus.30 <- otu_table(gluc.30)
write.csv(otus.30, file = "~/Desktop/r-summerproject/glucose-analysis/cohesion/otus-30.csv")

# 4. Construct the networks for each condition
########### 15C ############
gluc.15.net <- netConstruct(gluc.15,
                            measure = "pearson",
                            filtTax = "numbSamp",
                            filtTaxPar = list(numbSamp = 0.25),
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.3,
                            verbose = 2)

prop.gluc.15.net <- netAnalyze(gluc.15.net,
                               centrLCC = TRUE,
                               clustMethod = "cluster_walktrap",
                               hubPar = "eigenvector",
                               weightDeg = FALSE, normDeg = FALSE)

summary(prop.gluc.15.net)

plot.gluc.15 <- plot(prop.gluc.15.net,
                     layout = "layout_with_fr",
                     nodeColor = "cluster",
                     posCol = "#1A85FF",
                     negCol = "#D41159")
                    

############ 20C #############
gluc.20.net <- netConstruct(gluc.20,
                            measure = "pearson",
                            filtTax = "numbSamp",
                            filtTaxPar = list(numbSamp = 0.25),
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.3,
                            verbose = 2)

prop.gluc.20.net <- netAnalyze(gluc.20.net,
                               centrLCC = TRUE,
                               clustMethod = "cluster_walktrap",
                               hubPar = "eigenvector",
                               weightDeg = FALSE, normDeg = FALSE)

summary(prop.gluc.20.net)

plot.gluc.20 <- plot(prop.gluc.20.net,
                     layout = "layout_with_fr",
                     nodeColor = "cluster",
                     posCol = "#1A85FF",
                     negCol = "#D41159")

############ 25C #############
gluc.25.net <- netConstruct(gluc.25,
                            measure = "pearson",
                            filtTax = "numbSamp",
                            filtTaxPar = list(numbSamp = 0.25),
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.3,
                            verbose = 2)

prop.gluc.25.net <- netAnalyze(gluc.25.net,
                               centrLCC = TRUE,
                               clustMethod = "cluster_walktrap",
                               hubPar = "eigenvector",
                               weightDeg = FALSE, normDeg = FALSE)

summary(prop.gluc.25.net)

plot.gluc.25 <- plot(prop.gluc.25.net,
                     layout = "layout_with_fr",
                     nodeColor = "cluster",
                     posCol = "#1A85FF",
                     negCol = "#D41159")

############ 30C #############
gluc.30.net <- netConstruct(gluc.30,
                            measure = "pearson",
                            filtTax = "numbSamp",
                            filtTaxPar = list(numbSamp = 0.25),
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.3,
                            verbose = 2
)

prop.gluc.30.net <- netAnalyze(gluc.30.net,
                               centrLCC = TRUE,
                               clustMethod = "cluster_walktrap",
                               hubPar = "eigenvector",
                               weightDeg = FALSE, normDeg = FALSE)

summary(prop.gluc.30.net)

plot.gluc.30 <- plot(prop.gluc.30.net,
                     layout = "layout_with_fr",
                     nodeColor = "cluster",
                     posCol = "#1A85FF",
                     negCol = "#D41159")


# 7. Generate 100 random networks for each condition, 
# and test whether actual experimental values are non-random
# Define the conditions and corresponding OTU tables
conditions <- c("15", "20", "25", "30")
otu_tables <- list(
  "15" = gluc.15,
  "20" = gluc.20,
  "25" = gluc.25,
  "30" = gluc.30
)

# Initialize a list to store results for each condition
all_globalPropsList <- list()

# Loop through each condition
for (condition in conditions) {
  cat("Processing condition:", condition, "\n")  # Print progress
  
  # Get the corresponding OTU table
  gluc_physeq <- otu_tables[[condition]]
  
  # Get names of OTUs and other necessary info
  gluc_physeq_df <- as.data.frame(otu_table(gluc_physeq))
  num_OTUs <- sum(as(otu_table(gluc_physeq), "matrix") > 0)
  num_rows <- sum(sample_data(gluc_physeq)$condition == condition)
  num_reads <- sum(sample_sums(gluc_physeq))
  num_tables = 100
  
  # Create list to store generated OTU tables
  otu_tables_list <- list()
  
  # Generate 100 OTU tables
  set.seed(123)  # Set seed for reproducibility
  for (i in 1:num_tables) {
    selected_otus <- sample(rownames(gluc_physeq_df), size = num_OTUs, replace = TRUE)
    sampled_otu_matrix <- matrix(0, nrow = num_rows, ncol = num_OTUs)
    colnames(sampled_otu_matrix) <- selected_otus
    rownames(sampled_otu_matrix) <- paste0("Sample_", 1:num_rows)
    
    for (j in 1:num_rows) {
      proportions <- runif(num_OTUs)  # Generate random proportions
      proportions <- proportions / sum(proportions)  # Normalize to sum to 1
      sampled_otu_matrix[j, ] <- round(proportions * num_reads)  # Scale to total reads
    }
    
    otu_tables_list[[i]] <- sampled_otu_matrix  # Store the generated OTU table
  }
  
  # List to store network analysis results
  network_results <- list()
  
  # Loop over each OTU table to construct and analyze networks
  for (i in 1:num_tables) {
    cat("Processing OTU table", i, "\n")  # Print progress
    
    # Construct the network using netConstruct
    cat("Calling netConstruct for table", i, "\n")
    network <- tryCatch({
      netConstruct(
        data = otu_tables_list[[i]], 
        measure = "pearson",
        filtTax = "numbSamp", 
        filtTaxPar = list(numbSamp = 0.25),
        normMethod = "clr",           
        sparsMethod = "threshold",        
        thresh = 0.3                  
      )
    }, error = function(e) {
      cat("Error in netConstruct for table", i, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(network)) next  # Skip to the next table if network construction failed
    
    # Analyze the constructed network using netAnalyze
    cat("Calling netAnalyze for table", i, "\n")
    analysis <- tryCatch({
      netAnalyze(
        network, 
        clustMethod = "cluster_walktrap",
        hubPar = "eigenvector",
        weightDeg = FALSE, normDeg = FALSE
      )
    }, error = function(e) {
      cat("Error in netAnalyze for table", i, ":", e$message, "\n")
      return(NULL)
    })
    
    # Store the network and analysis results
    network_results[[i]] <- list(
      network = network,
      analysis = analysis
    )
  }
  
  # Initialize a list to store globalProps and globalPropsLCC
  globalPropsList <- list()
  
  # Loop through each result in network_results
  for (i in 1:length(network_results)) {
    # Extract the analysis result
    analysis <- network_results[[i]]$analysis
    
    # Initialize a list with default NA values
    globalProps <- list(
      nComp1 = NA, nComp2 = NA,
      avDiss1 = NA, avDiss2 = NA,
      avPath1 = NA, avPath2 = NA,
      clustCoef1 = NA, clustCoef2 = NA,
      modularity1 = NA, modularity2 = NA,
      vertConnect1 = NA, vertConnect2 = NA,
      edgeConnect1 = NA, edgeConnect2 = NA,
      natConnect1 = NA, natConnect2 = NA,
      density1 = NA, density2 = NA,
      pep1 = NA, pep2 = NA,
      lccSize1 = NA, lccSize2 = NA,
      lccSizeRel1 = NA, lccSizeRel2 = NA,
      avDissLCC1 = NA, avDissLCC2 = NA,
      avPathLCC1 = NA, avPathLCC2 = NA,
      clustCoefLCC1 = NA, clustCoefLCC2 = NA,
      modularityLCC1 = NA, modularityLCC2 = NA,
      vertConnectLCC1 = NA, vertConnectLCC2 = NA,
      edgeConnectLCC1 = NA, edgeConnectLCC2 = NA,
      natConnectLCC1 = NA, natConnectLCC2 = NA,
      densityLCC1 = NA, densityLCC2 = NA,
      pepLCC1 = NA, pepLCC2 = NA
    )
    
    # Check if analysis is not NULL and contains globalProps
    if (!is.null(analysis) && "globalProps" %in% names(analysis)) {
      gp <- analysis$globalProps
      
      # Replace NULL values with NA and fill globalProps
      globalProps$nComp1 <- ifelse(!is.null(gp$nComp1), gp$nComp1, NA)
      globalProps$nComp2 <- ifelse(!is.null(gp$nComp2), gp$nComp2, NA)
      globalProps$avDiss1 <- ifelse(!is.null(gp$avDiss1), gp$avDiss1, NA)
      globalProps$avDiss2 <- ifelse(!is.null(gp$avDiss2), gp$avDiss2, NA)
      globalProps$avPath1 <- ifelse(!is.null(gp$avPath1), gp$avPath1, NA)
      globalProps$avPath2 <- ifelse(!is.null(gp$avPath2), gp$avPath2, NA)
      globalProps$clustCoef1 <- ifelse(!is.null(gp$clustCoef1), gp$clustCoef1, NA)
      globalProps$clustCoef2 <- ifelse(!is.null(gp$clustCoef2), gp$clustCoef2, NA)
      globalProps$modularity1 <- ifelse(!is.null(gp$modularity1), gp$modularity1, NA)
      globalProps$modularity2 <- ifelse(!is.null(gp$modularity2), gp$modularity2, NA)
      globalProps$vertConnect1 <- ifelse(!is.null(gp$vertConnect1), gp$vertConnect1, NA)
      globalProps$vertConnect2 <- ifelse(!is.null(gp$vertConnect2), gp$vertConnect2, NA)
      globalProps$edgeConnect1 <- ifelse(!is.null(gp$edgeConnect1), gp$edgeConnect1, NA)
      globalProps$edgeConnect2 <- ifelse(!is.null(gp$edgeConnect2), gp$edgeConnect2, NA)
      globalProps$natConnect1 <- ifelse(!is.null(gp$natConnect1), gp$natConnect1, NA)
      globalProps$natConnect2 <- ifelse(!is.null(gp$natConnect2), gp$natConnect2, NA)
      globalProps$density1 <- ifelse(!is.null(gp$density1), gp$density1, NA)
      globalProps$density2 <- ifelse(!is.null(gp$density2), gp$density2, NA)
      globalProps$pep1 <- ifelse(!is.null(gp$pep1), gp$pep1, NA)
      globalProps$pep2 <- ifelse(!is.null(gp$pep2), gp$pep2, NA)
    }
    
    # Check if analysis contains globalPropsLCC
    if (!is.null(analysis) && "globalPropsLCC" %in% names(analysis)) {
      gpLCC <- analysis$globalPropsLCC
      
      # Replace NULL values with NA and fill globalPropsLCC
      globalProps$lccSize1 <- ifelse(!is.null(gpLCC$lccSize1), gpLCC$lccSize1, NA)
      globalProps$lccSize2 <- ifelse(!is.null(gpLCC$lccSize2), gpLCC$lccSize2, NA)
      globalProps$lccSizeRel1 <- ifelse(!is.null(gpLCC$lccSizeRel1), gpLCC$lccSizeRel1, NA)
      globalProps$lccSizeRel2 <- ifelse(!is.null(gpLCC$lccSizeRel2), gpLCC$lccSizeRel2, NA)
      globalProps$avDissLCC1 <- ifelse(!is.null(gpLCC$avDiss1), gpLCC$avDiss1, NA)
      globalProps$avDissLCC2 <- ifelse(!is.null(gpLCC$avDiss2), gpLCC$avDiss2, NA)
      globalProps$avPathLCC1 <- ifelse(!is.null(gpLCC$avPath1), gpLCC$avPath1, NA)
      globalProps$avPathLCC2 <- ifelse(!is.null(gpLCC$avPath2), gpLCC$avPath2, NA)
      globalProps$clustCoefLCC1 <- ifelse(!is.null(gpLCC$clustCoef1), gpLCC$clustCoef1, NA)
      globalProps$clustCoefLCC2 <- ifelse(!is.null(gpLCC$clustCoef2), gpLCC$clustCoef2, NA)
      globalProps$modularityLCC1 <- ifelse(!is.null(gpLCC$modularity1), gpLCC$modularity1, NA)
      globalProps$modularityLCC2 <- ifelse(!is.null(gpLCC$modularity2), gpLCC$modularity2, NA)
      globalProps$vertConnectLCC1 <- ifelse(!is.null(gpLCC$vertConnect1), gpLCC$vertConnect1, NA)
      globalProps$vertConnectLCC2 <- ifelse(!is.null(gpLCC$vertConnect2), gpLCC$vertConnect2, NA)
      globalProps$edgeConnectLCC1 <- ifelse(!is.null(gpLCC$edgeConnect1), gpLCC$edgeConnect1, NA)
      globalProps$edgeConnectLCC2 <- ifelse(!is.null(gpLCC$edgeConnect2), gpLCC$edgeConnect2, NA)
      globalProps$natConnectLCC1 <- ifelse(!is.null(gpLCC$natConnect1), gpLCC$natConnect1, NA)
      globalProps$natConnectLCC2 <- ifelse(!is.null(gpLCC$natConnect2), gpLCC$natConnect2, NA)
      globalProps$densityLCC1 <- ifelse(!is.null(gpLCC$density1), gpLCC$density1, NA)
      globalProps$densityLCC2 <- ifelse(!is.null(gpLCC$density2), gpLCC$density2, NA)
      globalProps$pepLCC1 <- ifelse(!is.null(gpLCC$pep1), gpLCC$pep1, NA)
      globalProps$pepLCC2 <- ifelse(!is.null(gpLCC$pep2), gpLCC$pep2, NA)
    }
    
    # Store the globalProps for the current OTU table
    globalPropsList[[i]] <- globalProps
  }
  
  # Store all globalProps for the current condition
  all_globalPropsList[[condition]] <- globalPropsList
}

# Now all_globalPropsList contains the results for each condition and each OTU table

# Initialize a list to store the rows for the dataframe
rows_list <- list()

# Loop through each condition and corresponding globalPropsList
for (condition in names(all_globalPropsList)) {
  globalPropsList <- all_globalPropsList[[condition]]
  
  # Loop through each set of global properties
  for (i in 1:length(globalPropsList)) {
    globalProps <- globalPropsList[[i]]
    
    # Create a data frame for this row
    df_row <- as.data.frame(t(unlist(globalProps)), stringsAsFactors = FALSE)
    
    # Add the condition column
    df_row$condition <- condition
    
    # Append this row to the list
    rows_list[[length(rows_list) + 1]] <- df_row
  }
}

# Combine all rows into a single dataframe
final_df <- bind_rows(rows_list)

# Rename columns to ensure they are descriptive
colnames(final_df) <- c(
  "nComp1", "nComp2", "avDiss1", "avDiss2", "avPath1", "avPath2",
  "clustCoef1", "clustCoef2", "modularity1", "modularity2",
  "vertConnect1", "vertConnect2", "edgeConnect1", "edgeConnect2",
  "natConnect1", "natConnect2", "density1", "density2",
  "pep1", "pep2", "lccSize1", "lccSize2", "lccSizeRel1", "lccSizeRel2",
  "avDissLCC1", "avDissLCC2", "avPathLCC1", "avPathLCC2",
  "clustCoefLCC1", "clustCoefLCC2", "modularityLCC1", "modularityLCC2",
  "vertConnectLCC1", "vertConnectLCC2", "edgeConnectLCC1", "edgeConnectLCC2",
  "natConnectLCC1", "natConnectLCC2", "densityLCC1", "densityLCC2",
  "pepLCC1", "pepLCC2", "condition"
)

# Print the final dataframe
print(final_df)

# Clean up dataframe 
gluc.glob.props <- final_df %>%
  select_if(~ !any(is.na(.)))

# CIs for each metric for each condition
boot.15 <- gluc.glob.props[gluc.glob.props$condition == 15, ]
write.csv(boot.15, "~/Desktop/r-summerproject/glucose-analysis/globprop-15.csv")
ci.apl.15 <- quantile(boot.15$avPath1, probs = c(0.025, 0.975))
cat("95% Confidence Interval for Average Path Length", ci.15, "\n")
ci.nc.15 <- quantile(boot.15$natConnect1, probs = c(0.025, 0.975))

# if plotting where it lies relative to bootstrapped
#ggplot(boot.15.apl, aes(x = avPath1)) +
  #geom_histogram(binwidth = 0.01, fill = "grey", color = "black", alpha = 0.7) +
  #geom_vline(xintercept = ci.15[1], linetype = "dashed", color = "orange") +
  #geom_vline(xintercept = ci.15[2], linetype = "dashed", color = "orange") +
  #geom_vline(xintercept = exp.15.apl, linetype = "solid", color = "blue") +
  #labs(x = "Average Path Length", y = "Frequency") +
  #theme_cowplot(12) +
  #scale_x_continuous(limits = c(min(boot.15.apl$avPath1) - 0.01, max(boot.15.apl$avPath1) + 0.01))

boot.20 <- gluc.glob.props[gluc.glob.props$condition == 20, ]
write.csv(boot.20, "~/Desktop/r-summerproject/glucose-analysis/globprop-20.csv")
ci.apl.20 <- quantile(boot.20$avPath1, probs = c(0.025, 0.975))
ci.nc.20 <- quantile(boot.20$natConnect1, probs = c(0.025, 0.975))
cat("95% Confidence Interval for Average Path Length", ci.ap20, "\n")

boot.25<- gluc.glob.props[gluc.glob.props$condition == 25, ]#
write.csv(boot.25, "~/Desktop/r-summerproject/glucose-analysis/globprop-25.csv")
ci.apl.25 <- quantile(boot.25$avPath1, probs = c(0.025, 0.975))
ci.nc.25 <- quantile(boot.25$natConnect1, probs = c(0.025, 0.975))
cat("95% Confidence Interval for Average Path Length", ci.25, "\n")

boot.30 <- gluc.glob.props[gluc.glob.props$condition == 30, ]
write.csv(boot.30, "~/Desktop/r-summerproject/glucose-analysis/globprop-30.csv")
ci.apl.30 <- quantile(boot.30$avPath1, probs = c(0.025, 0.975))
ci.nc.30 <- quantile(boot.30$avPath1, probs = c(0.025, 0.975))
cat("95% Confidence Interval for Average Path Length", ci.30, "\n")


