#############################
#I/O
#############################
library(tidyverse)
library(data.table)

# Read in sfari database
sfari <- read.csv("sfari_q4_2024.csv", header = TRUE)
x_linked_sfari_genes <- (sfari %>% filter(chromosome == "X"))$gene.symbol

# Read in omim data
omim <- fread("genemap_omim_030824.txt", 
               header = TRUE, fill = TRUE, sep = "\t", skip = 3) %>% as.data.frame()
colnames(omim)[1] <- "Chromosome"
omim <- omim %>% dplyr::select(Chromosome, 'Gene/Locus And Other Related Symbols', Phenotypes)
colnames(omim) <- c("Chromosome", "Gene", "Phenotype")
omim <- omim %>% filter(Chromosome == "chrX") 

#############################
#GET MODES OF INHERITANCE
#############################
# Filter OMIM dataframe to only contain entries involving one of the X-linked SFARI genes
omim_filtered <- omim %>%
  filter(sapply(Gene, function(g) any(sapply(x_linked_sfari_genes, grepl, g)))) %>%
  filter(Phenotype != "")

# Identify matched and unmatched genes
matched_genes <- x_linked_sfari_genes[sapply(x_linked_sfari_genes, function(g) any(grepl(g, omim_filtered$Gene)))]
unmatched_genes <- x_linked_sfari_genes[sapply(x_linked_sfari_genes, function(g) !any(grepl(g, omim_filtered$Gene)))]
num_matched <- length(matched_genes)
cat("Number of X-linked SFARI genes:", length(x_linked_sfari_genes), "\n")
cat("Number of X-linked SFARI genes found in OMIM:", num_matched, "\n")

# Assign explicit X-linked mode of inheritance
omim_filtered <- omim_filtered %>%
  mutate(
    X_linked_recessive = ifelse(grepl("X-linked recessive", Phenotype, ignore.case = TRUE), 1, 0),
    X_linked_dominant = ifelse(grepl("X-linked dominant", Phenotype, ignore.case = TRUE), 1, 0),
    X_linked_syndromic = ifelse(grepl("X-linked syndromic", Phenotype, ignore.case = TRUE), 1, 0)
  )

# Assign "other" X-linked category if no specific type is found
omim_filtered <- omim_filtered %>%
  mutate(
    X_linked_other = ifelse(
      grepl("X-linked", Phenotype, ignore.case = TRUE) &
      X_linked_recessive == 0 & 
      X_linked_dominant == 0 & 
      X_linked_syndromic == 0, 
      1, 0
    )
  )

# Make a row per gene 
omim_filtered <- omim_filtered %>%
  separate_rows(Gene, sep = ",\\s*")  # Splits by comma and optional spaces

# Create summary mode of inheritance per gene
omim_filtered <- omim_filtered %>%
  dplyr::select(Chromosome, Gene, X_linked_recessive, X_linked_dominant, X_linked_syndromic, X_linked_other) %>%
  unique() %>%
  mutate(
    mode_of_inheritance = case_when(
      X_linked_recessive == 1 & X_linked_dominant == 1 ~ "multiple reported modes of inheritance",  # If both XLR and XLD
      X_linked_recessive == 1 ~ "XLR",  
      X_linked_dominant == 1 ~ "XLD",  
      TRUE ~ "X-linked"                 # Otherwise, default to "X-linked"
    )
  )

# Assign final mode of inheritance per gene, keeping "multiple reported modes" when applicable
omim_filtered <- omim_filtered %>%
  group_by(Gene) %>%
  mutate(
    final_mode_of_inheritance = case_when(
      "multiple reported modes of inheritance" %in% mode_of_inheritance ~ "multiple reported modes of inheritance",
      "XLR" %in% mode_of_inheritance ~ "XLR",
      "XLD" %in% mode_of_inheritance ~ "XLD",
      TRUE ~ "X-linked"
    )
  ) %>%
  ungroup() %>%
  select(Chromosome, Gene, final_mode_of_inheritance) %>%
  distinct() %>%
  filter(Gene %in% x_linked_sfari_genes)

# Add back in genes that were not in OMIM 
omim_filtered <- bind_rows(omim_filtered, data.frame(
  Chromosome = "chrX", 
  Gene = unmatched_genes,
  final_mode_of_inheritance = "NotInOMIM")
)
setequal(omim_filtered$Gene, x_linked_sfari_genes)

#############################
#CREATE PIE CHART
#############################
# Convert the table into a data frame
omim_filtered$final_mode_of_inheritance <- sapply(omim_filtered$final_mode_of_inheritance, function(x) ifelse(x == "multiple reported modes of inheritance", "X-linked", x))
moi_df <- table(omim_filtered$final_mode_of_inheritance) %>% as.data.frame() 

# Set the factor levels to control the order of pie slices
colnames(moi_df)[1] <- "final_mode_of_inheritance"
moi_df$final_mode_of_inheritance <- factor(moi_df$final_mode_of_inheritance,
                                            levels = c("XLD", "XLR", "X-linked", "NotInOMIM"))

# Create the pie chart
pdf("asd_moi_piechart.pdf")
ggplot(moi_df, aes(x = "", y = Freq, fill = final_mode_of_inheritance)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove the background grid
  scale_fill_manual(values = c("XLR" = "#ffd670",
                               "X-linked" = "#E195AB", 
                               "XLD" = "#ff9770", 
                               "NotInOMIM" = "#70d6ff")) + 
  labs(title = "Proportion of Modes of Inheritance") +
  theme(legend.position = "bottom")
dev.off()

# Write out supplementary table 
write.table(omim_filtered %>% arrange(Gene), "asd_moi_supp_table.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)