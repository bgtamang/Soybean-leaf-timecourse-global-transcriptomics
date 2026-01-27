#Code to summarize Salmon output from quant.sf files and meta_info.json files

#run on worker node and in directory containing all Salmon output folders (one per sample):
# $ srun --pty bash
# $ module load R/4.4.0-IGB-gcc-8.2.0 
# $ Rscript pathToFile/Salmon_summarizeOutput.R

#output is SalmonSummarizedOutput.RData file with R objects tx.all (contains abundances (TPM), counts and estimated lengths)
#and meta_info containing Salmon's output summary info including salmon version and parameters, read numbers input and mapped.
#See https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats for information on Salmon's outputs

#If you also prefer the script to output tab-delimited text files as well, run the script like this:

# $ Rscript pathToFile/Salmon_summarizeOutput.R outputText


#R codes start here:

args <- commandArgs(TRUE)

library(tximport)
library(readr)
library(jsonlite)
library(ggplot2)
library(scales)  # For better axis formatting

options(stringsAsFactors=F)


# Get the full paths to the quant.sf files in each subdirectory

files.quant <- dir(pattern = "quant.sf$", recursive = T, full.names = T)


# create wrapper around read_tsv function to automatically read in  the Name as character
# and everything as double format (needed if only integer values in first 1000 rows)

my_read_tsv <- function(x, ...) {
    read_tsv(x, col_types = cols(Name = "c", .default = "d"), ...)
}


tx.all <- tximport(files.quant, type = "salmon", txOut = TRUE, dropInfReps = TRUE, 
                   importer = my_read_tsv)
               
#add in cleaned sample names

if (length(grep("_L00", files.quant)) == length(files.quant)) {
    temp <- strsplit(files.quant, "_")
    quant.good <- sapply(temp, function(x) {
        Lpos <- grep("L00", x)
        paste(x[1:(Lpos-2)], collapse = "_")
    })
    } else {
    quant.good <- gsub("/quant.sf", "", files.quant)
    }

quant.good <- make.names(gsub("./", "", quant.good))

colnames(tx.all$abundance) <- quant.good
colnames(tx.all$counts) <- quant.good
colnames(tx.all$length) <- quant.good
#also read in and summarize the info in the meta_info.json: Get the full paths to the files in each subdirectory

files.info <- dir(pattern = "meta_info.json", recursive = T, full.names = T)


# read in the first file

temp <- fromJSON(files.info[1], flatten = T)

meta_info <- as.data.frame(temp[sapply(temp, length) == 1])

#Now do for the rest of the samples
#skip this step if there's only one sample
if(length(files.info)>1){
for (i in 2:length(files.info)) {
    meta_info <- rbind(meta_info, as.data.frame(fromJSON(files.info[i], flatten = T)[sapply(temp, length) == 1]))
}}

#add in cleaned sample names

if (length(grep("_L00", files.info)) == length(files.info)) {
    temp <- strsplit(files.info, "_")
    files.good <- sapply(temp, function(x) {
        Lpos <- grep("L00", x)
        paste(x[1:(Lpos-2)], collapse = "_")
    })
    } else {
    files.good <- gsub("/aux_info/meta_info.json", "", files.info)
    }

files.good <- make.names(gsub("./", "", files.good))

rownames(meta_info) <- files.good


#write out tx_all and meta_info objects

save(meta_info, tx.all, file = "SalmonSummarizedOutput.RData")

#also write out as tab-delimited text files if specified:

if ("outputText" %in% args) {
  write.table(data.frame(sample = rownames(meta_info), meta_info), file = "Salmon_MetaInfo.txt", 
              row.names = FALSE, sep = "\t")
  write.table(data.frame(transcript_id = rownames(tx.all$abundance), tx.all$abundance),
              file = "Salmon_Tx_TPMabundances.txt",row.names = FALSE, sep = "\t") 
  write.table(data.frame(transcript_id = rownames(tx.all$counts), tx.all$counts),
              file = "Salmon_Tx_NumReads.txt", row.names = FALSE, sep = "\t")
  write.table(data.frame(transcript_id = rownames(tx.all$length), tx.all$length),
              file = "Salmon_Tx_EffectiveLengths.txt", row.names = FALSE, sep = "\t")
}


#make publication-quality pct mapped plot

meta_info$Sample <- factor(rownames(meta_info), levels = unique(rownames(meta_info)))

# Calculate mean and median for reference lines
mean_pct <- mean(meta_info$percent_mapped)
median_pct <- median(meta_info$percent_mapped)

# Create color gradient based on mapping percentage
meta_info$color_group <- cut(meta_info$percent_mapped, 
                             breaks = c(0, 50, 70, 85, 100),
                             labels = c("Poor (<50%)", "Fair (50-70%)", "Good (70-85%)", "Excellent (>85%)"),
                             include.lowest = TRUE)

# Enhanced publication-quality plot
p <- ggplot(data = meta_info, aes(x = Sample, y = percent_mapped, fill = color_group)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 0.3) +
  
  # Add reference lines
  geom_hline(yintercept = mean_pct, linetype = "dashed", color = "#2C3E50", size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = median_pct, linetype = "dotted", color = "#34495E", size = 0.8, alpha = 0.7) +
  
  # Color scale
  scale_fill_manual(values = c("Poor (<50%)" = "#E74C3C",
                               "Fair (50-70%)" = "#F39C12", 
                               "Good (70-85%)" = "#3498DB",
                               "Excellent (>85%)" = "#27AE60"),
                    name = "Mapping Quality") +
  
  # Y-axis formatting
  scale_y_continuous(breaks = seq(0, 100, by = 10), 
                     limits = c(0, 105),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  
  # Labels and title
  labs(title = "Transcriptome Mapping Efficiency by Salmon",
       subtitle = paste0("Mean: ", round(mean_pct, 1), "% | Median: ", round(median_pct, 1), "%"),
       x = "Sample",
       y = "Reads Mapped (%)") +
  
  # Theme customization
  theme_classic(base_size = 14) +
  theme(
    # Title formatting
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    
    # Axis formatting
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    # Legend formatting
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    
    # Panel and grid
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Margins
    plot.margin = margin(20, 20, 20, 20)
  )

# Add sample size annotation
p <- p + annotate("text", 
                  x = length(unique(meta_info$Sample)) * 0.95, 
                  y = 102, 
                  label = paste("n =", nrow(meta_info)), 
                  size = 4, 
                  hjust = 1,
                  fontface = "italic")

# Save as high-resolution PNG (better for publications than JPEG)
ggsave("SalmonPctMapped.png", 
       plot = p,
       width = 12, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")

# Also save as PDF for vector graphics (best for publications)
ggsave("SalmonPctMapped.pdf", 
       plot = p,
       width = 12, 
       height = 8, 
       units = "in",
       device = cairo_pdf)

# Optional: Save as TIFF if required by journal
ggsave("SalmonPctMapped.tiff", 
       plot = p,
       width = 12, 
       height = 8, 
       units = "in", 
       dpi = 300,
       compression = "lzw",
       bg = "white")

# Print summary statistics to console
cat("\n=== Mapping Statistics ===\n")
cat("Total samples:", nrow(meta_info), "\n")
cat("Mean mapping %:", round(mean_pct, 2), "\n")
cat("Median mapping %:", round(median_pct, 2), "\n")
cat("Min mapping %:", round(min(meta_info$percent_mapped), 2), "\n")
cat("Max mapping %:", round(max(meta_info$percent_mapped), 2), "\n")
cat("SD:", round(sd(meta_info$percent_mapped), 2), "\n")
cat("\nSamples by mapping quality:\n")
print(table(meta_info$color_group))
cat("\n")
