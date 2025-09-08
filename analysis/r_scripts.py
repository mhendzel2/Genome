"""
R script templates for genomics analysis
"""

# Basic statistics R script template
BASIC_STATS_SCRIPT = """
library(GenomicRanges)
library(dplyr)

# Function to compute basic statistics for genomics data
compute_basic_stats <- function(data_file) {
    # Read the data
    if (grepl("\\.bed$|\\.bedgraph$", data_file, ignore.case = TRUE)) {
        data <- read.table(data_file, sep = "\\t", header = FALSE, comment.char = "#")
        colnames(data)[1:3] <- c("chr", "start", "end")
    } else {
        data <- read.table(data_file, sep = "\\t", header = TRUE, comment.char = "#")
    }
    
    stats <- list()
    stats$total_regions <- nrow(data)
    
    # If genomic coordinates are present
    if ("start" %in% colnames(data) && "end" %in% colnames(data)) {
        lengths <- data$end - data$start
        stats$mean_length <- mean(lengths, na.rm = TRUE)
        stats$median_length <- median(lengths, na.rm = TRUE)
        stats$max_length <- max(lengths, na.rm = TRUE)
        stats$min_length <- min(lengths, na.rm = TRUE)
    }
    
    # If score column is present
    score_cols <- grep("score|signal|value", colnames(data), ignore.case = TRUE)
    if (length(score_cols) > 0) {
        scores <- data[, score_cols[1]]
        scores <- as.numeric(scores)
        scores <- scores[!is.na(scores)]
        
        if (length(scores) > 0) {
            stats$mean_score <- mean(scores)
            stats$median_score <- median(scores)
            stats$max_score <- max(scores)
            stats$min_score <- min(scores)
            stats$sd_score <- sd(scores)
        }
    }
    
    return(stats)
}

# Main execution
results <- NULL
# Default: try processing a single input placeholder named "input"
try({
    results <- compute_basic_stats("{file:input}")
}, silent = TRUE)

# If default processing didn't populate results, fall back to injected processing code
if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results as JSON
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# Peak calling R script template
PEAK_CALLING_SCRIPT = """
library(GenomicRanges)
library(ChIPseeker)
library(dplyr)

# Function to call peaks from ChIP-seq data
call_peaks <- function(data_file, p_threshold = 0.05, fold_change = 2.0) {
    # Read ChIP-seq data
    data <- read.table(data_file, sep = "\\t", header = FALSE, comment.char = "#")
    
    # Assume BED format: chr, start, end, [name], [score]
    if (ncol(data) >= 3) {
        colnames(data)[1:3] <- c("chr", "start", "end")
        if (ncol(data) >= 5) {
            colnames(data)[5] <- "score"
        }
    }
    
    peaks <- list()
    
    if ("score" %in% colnames(data)) {
        # Filter peaks based on score threshold
        significant_peaks <- data[data$score > fold_change, ]
        
        if (nrow(significant_peaks) > 0) {
            for (i in 1:nrow(significant_peaks)) {
                peak <- list(
                    chr = significant_peaks$chr[i],
                    start = significant_peaks$start[i],
                    end = significant_peaks$end[i],
                    score = significant_peaks$score[i],
                    length = significant_peaks$end[i] - significant_peaks$start[i]
                )
                peaks[[i]] <- peak
            }
        }
    }
    
    return(peaks)
}

# Main execution
results <- NULL
# Default: try processing a single input placeholder named "input"
try({
    results <- call_peaks("{file:input}")
}, silent = TRUE)

if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# Differential expression R script template
DIFF_EXPRESSION_SCRIPT = """
library(DESeq2)
library(dplyr)
library(GenomicRanges)

# Function for differential expression analysis
perform_diff_expression <- function(count_matrix, sample_info, padj_threshold = 0.05, lfc_threshold = 1.0) {
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_info,
        design = ~ condition
    )
    
    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds, alpha = padj_threshold)
    
    # Filter significant genes
    sig_genes <- res[which(res$padj < padj_threshold & abs(res$log2FoldChange) > lfc_threshold), ]
    
    # Count up/down regulated genes
    upregulated <- sum(sig_genes$log2FoldChange > 0, na.rm = TRUE)
    downregulated <- sum(sig_genes$log2FoldChange < 0, na.rm = TRUE)
    
    results <- list(
        total_genes = nrow(res),
        significant_genes = nrow(sig_genes),
        upregulated = upregulated,
        downregulated = downregulated,
        padj_threshold = padj_threshold,
        lfc_threshold = lfc_threshold
    )
    
    return(results)
}

# Main execution
results <- NULL
# Default: expect two placeholders: counts and samples
try({
    count_matrix <- as.matrix(read.table("{file:counts}", sep = "\t", header = TRUE, row.names = 1))
    sample_info <- read.table("{file:samples}", sep = "\t", header = TRUE, row.names = 1)
    results <- perform_diff_expression(count_matrix, sample_info)
}, silent = TRUE)

if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# Tissue comparison R script template
TISSUE_COMPARISON_SCRIPT = """
library(GenomicRanges)
library(VennDiagram)
library(dplyr)

# Function to compare genomic features between tissues
compare_tissues <- function(tissue1_file, tissue2_file) {
    # Read tissue data
    tissue1_data <- read.table(tissue1_file, sep = "\\t", header = FALSE, comment.char = "#")
    tissue2_data <- read.table(tissue2_file, sep = "\\t", header = FALSE, comment.char = "#")
    
    # Create genomic ranges
    if (ncol(tissue1_data) >= 3) {
        colnames(tissue1_data)[1:3] <- c("chr", "start", "end")
        tissue1_gr <- GRanges(
            seqnames = tissue1_data$chr,
            ranges = IRanges(start = tissue1_data$start, end = tissue1_data$end)
        )
    }
    
    if (ncol(tissue2_data) >= 3) {
        colnames(tissue2_data)[1:3] <- c("chr", "start", "end")
        tissue2_gr <- GRanges(
            seqnames = tissue2_data$chr,
            ranges = IRanges(start = tissue2_data$start, end = tissue2_data$end)
        )
    }
    
    # Find overlaps
    overlaps <- findOverlaps(tissue1_gr, tissue2_gr)
    
    # Calculate statistics
    common_features <- length(unique(queryHits(overlaps)))
    unique_tissue1 <- length(tissue1_gr) - common_features
    unique_tissue2 <- length(tissue2_gr) - length(unique(subjectHits(overlaps)))
    
    results <- list(
        common_features = common_features,
        unique_tissue1 = unique_tissue1,
        unique_tissue2 = unique_tissue2,
        total_tissue1 = length(tissue1_gr),
        total_tissue2 = length(tissue2_gr),
        overlap_percentage = (common_features / length(tissue1_gr)) * 100
    )
    
    return(results)
}

# Main execution
results <- NULL
# Default: expect two placeholders named tissue1 and tissue2
try({
    results <- compare_tissues("{file:tissue1}", "{file:tissue2}")
}, silent = TRUE)

if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# HiC analysis R script template
HIC_ANALYSIS_SCRIPT = """
library(GenomicRanges)
library(HiCcompare)
library(dplyr)

# Function to analyze Hi-C interaction data
analyze_hic_interactions <- function(hic_file, resolution = 10000) {
    # Read Hi-C data (assuming BEDPE format)
    hic_data <- read.table(hic_file, sep = "\\t", header = FALSE, comment.char = "#")
    
    # Basic validation
    if (ncol(hic_data) < 6) {
        stop("Hi-C data should have at least 6 columns (BEDPE format)")
    }
    
    colnames(hic_data)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    
    # Calculate interaction distances
    intra_chr <- hic_data$chr1 == hic_data$chr2
    distances <- ifelse(intra_chr, 
                       abs((hic_data$start1 + hic_data$end1) / 2 - (hic_data$start2 + hic_data$end2) / 2),
                       NA)
    
    # Basic statistics
    results <- list(
        total_interactions = nrow(hic_data),
        intra_chromosomal = sum(intra_chr),
        inter_chromosomal = sum(!intra_chr),
        mean_distance = mean(distances, na.rm = TRUE),
        median_distance = median(distances, na.rm = TRUE),
        max_distance = max(distances, na.rm = TRUE),
        min_distance = min(distances, na.rm = TRUE)
    )
    
    # Identify long-range interactions (>1Mb)
    long_range <- distances > 1000000
    results$long_range_interactions <- sum(long_range, na.rm = TRUE)
    
    return(results)
}

# Main execution
results <- NULL
# Default: expect a placeholder named hic
try({
    results <- analyze_hic_interactions("{file:hic}")
}, silent = TRUE)

if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# Quality control R script template
QC_ANALYSIS_SCRIPT = """
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# Function to perform quality control on genomics data
perform_qc <- function(data_file, data_type = "bed") {
    # Read data
    data <- read.table(data_file, sep = "\\t", header = FALSE, comment.char = "#")
    
    qc_results <- list()
    qc_results$file_name <- basename(data_file)
    qc_results$total_lines <- nrow(data)
    
    # Basic format validation
    if (data_type %in% c("bed", "bedgraph")) {
        if (ncol(data) >= 3) {
            colnames(data)[1:3] <- c("chr", "start", "end")
            
            # Check coordinate validity
            valid_coords <- data$start < data$end & data$start >= 0
            qc_results$valid_coordinates <- sum(valid_coords)
            qc_results$invalid_coordinates <- sum(!valid_coords)
            qc_results$coordinate_error_rate <- mean(!valid_coords)
            
            # Check chromosome naming
            chr_pattern <- grepl("^(chr)?[0-9XYM]+$", data$chr, ignore.case = TRUE)
            qc_results$valid_chromosomes <- sum(chr_pattern)
            qc_results$invalid_chromosomes <- sum(!chr_pattern)
            
            # Feature length distribution
            lengths <- data$end - data$start
            qc_results$mean_feature_length <- mean(lengths, na.rm = TRUE)
            qc_results$median_feature_length <- median(lengths, na.rm = TRUE)
            qc_results$min_feature_length <- min(lengths, na.rm = TRUE)
            qc_results$max_feature_length <- max(lengths, na.rm = TRUE)
        }
        
        # Score column QC if present
        if (ncol(data) >= 4) {
            scores <- as.numeric(data[, 4])
            scores <- scores[!is.na(scores)]
            
            if (length(scores) > 0) {
                qc_results$score_mean <- mean(scores)
                qc_results$score_median <- median(scores)
                qc_results$score_range <- max(scores) - min(scores)
                qc_results$score_na_count <- sum(is.na(data[, 4]))
            }
        }
    }
    
    # Overall quality assessment
    if (qc_results$coordinate_error_rate < 0.01) {
        qc_results$overall_quality <- "PASS"
    } else if (qc_results$coordinate_error_rate < 0.05) {
        qc_results$overall_quality <- "WARNING"
    } else {
        qc_results$overall_quality <- "FAIL"
    }
    
    return(qc_results)
}

# Main execution
results <- NULL
# Default: try processing a single input placeholder named "input"
try({
    results <- perform_qc("{file:input}")
}, silent = TRUE)

if (is.null(results) || length(results) == 0) {
{file_processing_code}
}

# Output results
cat(jsonlite::toJSON(results, auto_unbox = TRUE))
"""

# Script templates dictionary
R_SCRIPT_TEMPLATES = {
    'basic_statistics': BASIC_STATS_SCRIPT,
    'peak_calling': PEAK_CALLING_SCRIPT,
    'differential_expression': DIFF_EXPRESSION_SCRIPT,
    'tissue_comparison': TISSUE_COMPARISON_SCRIPT,
    'hic_analysis': HIC_ANALYSIS_SCRIPT,
    'quality_control': QC_ANALYSIS_SCRIPT
}

def get_r_script(analysis_type: str, **kwargs) -> str:
    """Get R script template for specific analysis type"""
    if analysis_type not in R_SCRIPT_TEMPLATES:
        raise ValueError(f"Unknown analysis type: {analysis_type}")
    
    template = R_SCRIPT_TEMPLATES[analysis_type]
    
    # Replace placeholders with actual processing code
    file_processing_code = kwargs.get('file_processing_code', '# File processing code here')
    
    return template.format(file_processing_code=file_processing_code)

def generate_file_processing_code(files_dict: dict) -> str:
    """Generate R code for processing multiple files"""
    code_lines = []
    code_lines.append("results <- list()")
    
    for i, (filename, file_info) in enumerate(files_dict.items()):
        code_lines.append(f"# Processing file: {filename}")
        code_lines.append(f"file_{i} <- '{filename}'")
        code_lines.append(f"if (file.exists(file_{i})) {{")
        code_lines.append(f"    results[['{filename}']] <- compute_basic_stats(file_{i})")
        code_lines.append("} else {")
        code_lines.append(f"    results[['{filename}']] <- list(error = 'File not found')")
        code_lines.append("}")
        code_lines.append("")
    
    return "\n".join(code_lines)
