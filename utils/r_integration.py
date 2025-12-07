import os
import subprocess
import tempfile
import json
from typing import Dict, Any, List, Optional
import pandas as pd
import logging

class RIntegration:
    """Handles R integration for genomics analysis"""
    
    def __init__(self):
        self.r_available = self._check_r_availability()
        self.required_packages = [
            'GenomicRanges',
            'DESeq2', 
            'ChIPseeker',
            'HiCcompare',
            'BiocManager',
            'dplyr',
            'ggplot2'
        ]
        
        if self.r_available:
            self._install_required_packages()
    
    def _check_r_availability(self) -> bool:
        """Check if R is available on the system"""
        try:
            result = subprocess.run(['R', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def _install_required_packages(self):
        """Install required R packages"""
        if not self.r_available:
            return
        
        install_script = '''
        # Install BiocManager if not available
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos="https://cran.r-project.org")
        }
        
        # Required packages
        required_packages <- c("GenomicRanges", "DESeq2", "ChIPseeker", "dplyr", "ggplot2")
        
        # Install missing packages
        for (pkg in required_packages) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                if (pkg %in% c("GenomicRanges", "DESeq2", "ChIPseeker")) {
                    BiocManager::install(pkg, ask = FALSE, update = FALSE)
                } else {
                    install.packages(pkg, repos="https://cran.r-project.org")
                }
            }
        }
        
        cat("Package installation completed\\n")
        '''
        
        try:
            self._execute_r_script(install_script)
        except Exception as e:
            print(f"Warning: Could not install R packages: {e}")
    
    def _execute_r_script(self, script: str, data_files: Optional[Dict[str, str]] = None,
                         args_files: Optional[List[str]] = None, extra_args: Optional[List[str]] = None) -> str:
        """
        Execute R script safely with input validation.

        Args:
            script: The R script content.
            data_files: Dictionary mapping filenames to content (string).
            args_files: List of keys in data_files whose temporary paths should be passed as arguments.
            extra_args: List of additional string arguments to pass to the script.
        """
        if not self.r_available:
            raise RuntimeError("R is not available on this system")
        
        if not self._validate_r_script(script):
            raise ValueError("R script contains potentially unsafe content")
        
        logger = logging.getLogger(__name__)

        with tempfile.TemporaryDirectory() as temp_dir:
            # Write data files to temporary paths
            file_map = {}
            if data_files:
                for filename, content in data_files.items():
                    safe_filename = self._sanitize_filename(filename)
                    file_path = os.path.join(temp_dir, safe_filename)
                    with open(file_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                    file_map[filename] = file_path

            # Allow scripts to reference input files via a simple placeholder syntax {file:<name>}
            # (Deprecated but kept for backward compatibility if needed, though we prefer args now)
            for original_name, path in file_map.items():
                placeholder = f"{{file:{original_name}}}"
                script = script.replace(placeholder, path.replace('\\', '/'))

            # Write R script to temporary file with safe filename
            script_path = os.path.join(temp_dir, 'analysis.R')
            with open(script_path, 'w', encoding='utf-8') as f:
                f.write(script)

            # Construct command line arguments
            cmd_args = []
            if args_files:
                for key in args_files:
                    if key in file_map:
                        cmd_args.append(file_map[key])
                    else:
                        logger.warning(f"File key '{key}' not found in data_files map")

            if extra_args:
                cmd_args.extend([str(arg) for arg in extra_args])

            # Execute R script with restricted environment
            # Uses 'Rscript' to easily pass arguments
            try:
                cmd = ['Rscript', '--vanilla', script_path] + cmd_args

                result = subprocess.run(
                    cmd,
                    cwd=temp_dir,
                    capture_output=True,
                    text=True,
                    timeout=300,  # 5 minute timeout
                    env={'PATH': os.environ.get('PATH', ''), 'HOME': temp_dir}  # Restricted environment
                )

                if result.returncode != 0:
                    logger.error('R script stderr: %s', result.stderr)
                    raise RuntimeError(f"R script execution failed: {result.stderr}")

                return result.stdout

            except subprocess.TimeoutExpired:
                logger.exception('R script execution timed out')
                raise RuntimeError("R script execution timed out")
    
    def _validate_r_script(self, script: str) -> bool:
        """Validate R script for security issues"""
        dangerous_patterns = [
            'system(', 'shell(', 'Sys.setenv(', 'setwd(',
            'file.remove(', 'unlink(', 'file.create(',
            'source(', 'eval(', 'parse('
        ]
        
        script_lower = script.lower()
        for pattern in dangerous_patterns:
            if pattern.lower() in script_lower:
                return False
        
        return True
    
    def _sanitize_filename(self, filename: str) -> str:
        """Sanitize filename to prevent path traversal attacks"""
        import re
        safe_name = re.sub(r'[^\w\-_\.]', '_', filename)
        safe_name = re.sub(r'^[\.\-]+', '', safe_name)
        return safe_name[:100] if safe_name else 'data_file.txt'
    
    def basic_statistics(self, data_files: Dict[str, Any]) -> Dict[str, Any]:
        """Compute basic statistics using R"""
        if not self.r_available:
            return self._python_basic_statistics(data_files)
        
        results = {}
        
        for filename, file_info in data_files.items():
            try:
                # Get content
                content = self._get_file_content(file_info)

                r_script = '''
                library(jsonlite)
                args <- commandArgs(trailingOnly = TRUE)
                input_file <- args[1]

                # Read data from provided file path
                data <- read.table(file = input_file, sep = "\t", header = FALSE, comment.char = "#")

                stats <- list()
                stats$total_regions <- nrow(data)

                if (ncol(data) >= 4) {
                    scores <- as.numeric(data[,4])
                    scores <- scores[!is.na(scores)]
                    if (length(scores) > 0) {
                        stats$mean_score <- mean(scores)
                        stats$median_score <- median(scores)
                        stats$max_score <- max(scores)
                        stats$min_score <- min(scores)
                        stats$sd_score <- sd(scores)
                    }
                }

                if (ncol(data) >= 3) {
                    starts <- as.numeric(data[,2])
                    ends <- as.numeric(data[,3])
                    lengths <- ends - starts
                    lengths <- lengths[!is.na(lengths) & lengths > 0]
                    if (length(lengths) > 0) {
                        stats$mean_length <- mean(lengths)
                        stats$median_length <- median(lengths)
                        stats$max_length <- max(lengths)
                        stats$min_length <- min(lengths)
                    }
                }

                cat(toJSON(stats, auto_unbox = TRUE))
                '''

                output = self._execute_r_script(
                    r_script,
                    data_files={filename: content},
                    args_files=[filename]
                )

                try:
                    stats = json.loads(output.strip())
                    results[filename] = stats
                except Exception:
                    results[filename] = {'total_regions': 'N/A', 'error': 'Could not parse R output'}
                
            except Exception as e:
                results[filename] = {'error': f'Analysis failed: {str(e)}'}
        
        return results
    
    def peak_calling(self, data_files: Dict[str, Any], **params) -> List[Dict[str, Any]]:
        """Perform peak calling analysis"""
        if not self.r_available:
            return self._python_peak_calling(data_files, **params)

        logger = logging.getLogger(__name__)
        results = []
        threshold = float(params.get('fold_change', 2.0))

        for filename, file_info in data_files.items():
            if 'chip' in filename.lower() or file_info.get('type') == 'ChIP-seq':
                try:
                    content = self._get_file_content(file_info)

                    r_script = '''
                    library(jsonlite)
                    args <- commandArgs(trailingOnly = TRUE)
                    input_file <- args[1]
                    threshold <- as.numeric(args[2])

                    dat <- read.table(file = input_file, sep = "\t", header = FALSE, comment.char = "#")
                    peaks <- list()
                    if (ncol(dat) >= 4) {
                        # Limit to first 100 for example
                        for (i in seq_len(min(nrow(dat), 100))) {
                            row <- dat[i, ]
                            score <- as.numeric(row[[4]])
                            if (!is.na(score) && score > threshold) {
                                peaks[[length(peaks) + 1]] <- list(
                                    chr = as.character(row[[1]]),
                                    start = as.integer(row[[2]]),
                                    end = as.integer(row[[3]]),
                                    score = score
                                )
                            }
                        }
                    }
                    cat(toJSON(peaks, auto_unbox = TRUE))
                    '''

                    output = self._execute_r_script(
                        r_script,
                        data_files={filename: content},
                        args_files=[filename],
                        extra_args=[str(threshold)]
                    )

                    parsed = json.loads(output.strip()) if output and output.strip() else []
                    for p in parsed:
                        p['file'] = filename
                        if 'length' not in p and 'start' in p and 'end' in p:
                            p['length'] = p.get('end', 0) - p.get('start', 0)
                        results.append(p)
                except Exception as e:
                    logger.exception('Error calling peaks for %s: %s', filename, e)

        return results

    def differential_expression(self, data_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """Perform differential expression analysis"""
        if not self.r_available:
            return self._python_differential_expression(data_files, **params)

        logger = logging.getLogger(__name__)
        stats = {'upregulated': 0, 'downregulated': 0, 'total_genes': 0, 'significant_genes': 0}

        for filename, file_info in data_files.items():
            if 'expression' in filename.lower() or file_info.get('type') == 'Gene Expression':
                try:
                    content = self._get_file_content(file_info)

                    r_script = '''
                    library(jsonlite)
                    args <- commandArgs(trailingOnly = TRUE)
                    input_file <- args[1]

                    dat <- read.table(file = input_file, sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
                    res <- list()
                    res$total_genes <- nrow(dat)
                    if ('log2FoldChange' %in% colnames(dat)) {
                        up <- sum(dat$log2FoldChange > 0, na.rm = TRUE)
                        down <- sum(dat$log2FoldChange < 0, na.rm = TRUE)
                        res$upregulated <- up
                        res$downregulated <- down
                        res$significant_genes <- up + down
                    } else {
                        res$upregulated <- 0
                        res$downregulated <- 0
                        res$significant_genes <- 0
                    }
                    cat(toJSON(res, auto_unbox = TRUE))
                    '''

                    output = self._execute_r_script(
                        r_script,
                        data_files={filename: content},
                        args_files=[filename]
                    )

                    parsed = json.loads(output.strip())
                    stats.update(parsed)
                except Exception as e:
                    logger.exception('DE failed for %s: %s', filename, e)

        return stats
    
    def _python_differential_expression(self, data_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """
        Python implementation of differential expression analysis.
        
        Performs t-tests between control and treatment groups with FDR correction.
        Used as fallback when R is not available.
        
        Args:
            data_files: Dictionary of filename -> file_info containing expression data
            **params: Additional parameters:
                - fold_change_cutoff: Minimum fold change threshold (default 2.0)
                - p_value_cutoff: P-value significance threshold (default 0.05)
                - control_pattern: Pattern to identify control samples
                - treatment_pattern: Pattern to identify treatment samples
        
        Returns:
            Dictionary with differential expression results
        """
        import numpy as np
        from scipy import stats
        from statsmodels.stats.multitest import fdrcorrection
        
        logger = logging.getLogger(__name__)
        
        fold_change_cutoff = params.get('fold_change_cutoff', 2.0)
        p_value_cutoff = params.get('p_value_cutoff', 0.05)
        control_pattern = params.get('control_pattern', 'control')
        treatment_pattern = params.get('treatment_pattern', 'treatment')
        
        result = {
            'total_genes': 0,
            'tested_genes': 0,
            'upregulated': 0,
            'downregulated': 0,
            'significant_genes': 0,
            'method': 'python_ttest',
            'genes': [],
            'pvalue': [],
            'padj': [],
            'log2FoldChange': []
        }
        
        # Categorize files into control and treatment
        control_files = []
        treatment_files = []
        expression_data = {}
        
        for filename, file_info in data_files.items():
            if file_info.get('type') != 'Gene Expression':
                continue
            
            try:
                content = self._get_file_content(file_info)
                
                # Parse expression data
                lines = [l.strip() for l in content.split('\n') if l.strip() and not l.startswith('#')]
                
                if not lines:
                    continue
                
                # Detect header
                first_line = lines[0].split('\t')
                has_header = not first_line[1].replace('.', '').replace('-', '').isdigit() if len(first_line) > 1 else False
                
                # Parse data
                start_idx = 1 if has_header else 0
                genes = {}
                
                for line in lines[start_idx:]:
                    fields = line.split('\t')
                    if len(fields) >= 2:
                        gene = fields[0]
                        try:
                            value = float(fields[1])
                            genes[gene] = value
                        except ValueError:
                            continue
                
                # Determine if control or treatment based on filename
                filename_lower = filename.lower()
                if control_pattern.lower() in filename_lower:
                    control_files.append(filename)
                    expression_data[f'control_{len(control_files)}'] = genes
                elif treatment_pattern.lower() in filename_lower:
                    treatment_files.append(filename)
                    expression_data[f'treatment_{len(treatment_files)}'] = genes
                else:
                    # Default to treatment if can't determine
                    treatment_files.append(filename)
                    expression_data[f'treatment_{len(treatment_files)}'] = genes
                    
            except Exception as e:
                logger.warning(f"Error processing {filename}: {e}")
        
        # If only one file, check if it contains pre-computed DE results
        if len(data_files) == 1 and (not control_files or not treatment_files):
            for filename, file_info in data_files.items():
                try:
                    content = self._get_file_content(file_info)
                    lines = [l.strip() for l in content.split('\n') if l.strip() and not l.startswith('#')]
                    if lines:
                        header = lines[0].lower().split('\t')
                        
                        # Check for pre-computed DE results
                        if 'log2foldchange' in header or 'logfc' in header:
                            return self._parse_precomputed_de(lines, params)
                except Exception:
                    pass
        
        # Perform DE analysis if we have both groups
        if not control_files or not treatment_files:
            result['error'] = 'Need both control and treatment files for DE analysis'
            logger.warning("Insufficient files for DE: controls=%d, treatments=%d", 
                          len(control_files), len(treatment_files))
            return result
        
        # Find common genes
        all_genes = set()
        for genes in expression_data.values():
            all_genes.update(genes.keys())
        
        result['total_genes'] = len(all_genes)
        
        # Build expression matrix
        control_keys = [k for k in expression_data.keys() if k.startswith('control_')]
        treatment_keys = [k for k in expression_data.keys() if k.startswith('treatment_')]
        
        de_results = []
        
        for gene in all_genes:
            control_values = [expression_data[k].get(gene) for k in control_keys 
                            if expression_data[k].get(gene) is not None]
            treatment_values = [expression_data[k].get(gene) for k in treatment_keys 
                              if expression_data[k].get(gene) is not None]
            
            # Need at least 2 values in each group for t-test
            if len(control_values) < 2 or len(treatment_values) < 2:
                continue
            
            control_mean = np.mean(control_values)
            treatment_mean = np.mean(treatment_values)
            
            # Avoid division by zero and log of zero
            if control_mean <= 0:
                control_mean = 1e-10
            if treatment_mean <= 0:
                treatment_mean = 1e-10
            
            # Calculate log2 fold change
            log2fc = np.log2(treatment_mean / control_mean)
            
            # Perform t-test
            try:
                t_stat, p_value = stats.ttest_ind(control_values, treatment_values)
                if np.isnan(p_value):
                    continue
            except Exception:
                continue
            
            de_results.append({
                'gene': gene,
                'control_mean': control_mean,
                'treatment_mean': treatment_mean,
                'log2FoldChange': log2fc,
                'pvalue': p_value
            })
        
        result['tested_genes'] = len(de_results)
        
        if not de_results:
            return result
        
        # Apply FDR correction
        p_values = [r['pvalue'] for r in de_results]
        _, padj = fdrcorrection(p_values)
        
        for i, r in enumerate(de_results):
            r['padj'] = padj[i]
        
        # Count significant genes
        for r in de_results:
            if r['padj'] < p_value_cutoff and abs(r['log2FoldChange']) >= np.log2(fold_change_cutoff):
                if r['log2FoldChange'] > 0:
                    result['upregulated'] += 1
                else:
                    result['downregulated'] += 1
        
        result['significant_genes'] = result['upregulated'] + result['downregulated']
        
        # Sort by adjusted p-value and store top results
        de_results.sort(key=lambda x: x['padj'])
        
        result['genes'] = [r['gene'] for r in de_results[:1000]]
        result['pvalue'] = [r['pvalue'] for r in de_results[:1000]]
        result['padj'] = [r['padj'] for r in de_results[:1000]]
        result['log2FoldChange'] = [r['log2FoldChange'] for r in de_results[:1000]]
        
        logger.info(f"Python DE analysis: {result['tested_genes']} genes tested, "
                   f"{result['significant_genes']} significant (up={result['upregulated']}, down={result['downregulated']})")
        
        return result
    
    def _parse_precomputed_de(self, lines: list, params: dict) -> Dict[str, Any]:
        """Parse pre-computed differential expression results."""
        import numpy as np
        
        fold_change_cutoff = params.get('fold_change_cutoff', 2.0)
        p_value_cutoff = params.get('p_value_cutoff', 0.05)
        
        result = {
            'total_genes': 0,
            'upregulated': 0,
            'downregulated': 0,
            'significant_genes': 0,
            'method': 'precomputed',
            'genes': [],
            'pvalue': [],
            'padj': [],
            'log2FoldChange': []
        }
        
        header = lines[0].lower().split('\t')
        
        # Find relevant columns
        gene_col = None
        fc_col = None
        pval_col = None
        padj_col = None
        
        for i, h in enumerate(header):
            if h in ['gene', 'gene_id', 'symbol', 'name', 'id']:
                gene_col = i
            elif h in ['log2foldchange', 'logfc', 'log2fc', 'fc']:
                fc_col = i
            elif h in ['pvalue', 'pval', 'p_value', 'p.value']:
                pval_col = i
            elif h in ['padj', 'fdr', 'q_value', 'qvalue', 'adj.p.value', 'adjusted_pvalue']:
                padj_col = i
        
        if gene_col is None:
            gene_col = 0  # Assume first column is gene
        
        for line in lines[1:]:
            fields = line.split('\t')
            if len(fields) <= max(filter(None, [gene_col, fc_col, pval_col, padj_col]), default=0):
                continue
            
            gene = fields[gene_col] if gene_col is not None else ''
            
            try:
                log2fc = float(fields[fc_col]) if fc_col is not None else 0
                pval = float(fields[pval_col]) if pval_col is not None else 1
                padj = float(fields[padj_col]) if padj_col is not None else pval
            except (ValueError, IndexError):
                continue
            
            result['total_genes'] += 1
            result['genes'].append(gene)
            result['log2FoldChange'].append(log2fc)
            result['pvalue'].append(pval)
            result['padj'].append(padj)
            
            # Count significant
            if padj < p_value_cutoff and abs(log2fc) >= np.log2(fold_change_cutoff):
                if log2fc > 0:
                    result['upregulated'] += 1
                else:
                    result['downregulated'] += 1
        
        result['significant_genes'] = result['upregulated'] + result['downregulated']
        
        return result

    def tissue_comparison(self, tissues: List[str], comparison_type: str, data_files: Dict[str, Any]) -> Dict[str, Any]:
        """Perform tissue comparison analysis"""
        logger = logging.getLogger(__name__)
        results = {'tissues_compared': tissues, 'comparison_type': comparison_type, 'common_features': 0, 'unique_tissue1': 0, 'unique_tissue2': 0, 'total_features': 0}

        try:
            payload = {}
            for filename, file_info in data_files.items():
                try:
                    content = self._get_file_content(file_info)
                    payload[filename] = content
                except Exception as e:
                    logger.exception('Could not read %s: %s', filename, e)

            if len(payload) >= 2:
                names = list(payload.keys())
                f1 = names[0]
                f2 = names[1]

                r_script = '''
                    library(jsonlite)
                    args <- commandArgs(trailingOnly = TRUE)
                    file1 <- args[1]
                    file2 <- args[2]

                    a <- read.table(file = file1, sep = "\t", header = FALSE, comment.char = "#")
                    b <- read.table(file = file2, sep = "\t", header = FALSE, comment.char = "#")
                    fa <- apply(a[,1:3], 1, function(x) paste0(x[1], ':', x[2], '-', x[3]))
                    fb <- apply(b[,1:3], 1, function(x) paste0(x[1], ':', x[2], '-', x[3]))
                    common <- length(intersect(fa, fb))
                    res <- list(common_features = common, total_a = length(fa), total_b = length(fb))
                    cat(toJSON(res, auto_unbox = TRUE))
                '''

                output = self._execute_r_script(
                    r_script,
                    data_files=payload,
                    args_files=[f1, f2]
                )

                parsed = json.loads(output.strip())
                results.update({'common_features': parsed.get('common_features', 0), 'total_features': parsed.get('total_a', 0) + parsed.get('total_b', 0)})
        except Exception as e:
            logger.exception('Tissue comparison failed: %s', e)

        return results

    def _get_file_content(self, file_info):
        """Helper to extract content string from file info"""
        try:
            file_info['file'].seek(0)
        except Exception:
            pass
        raw = file_info['file'].read()
        if isinstance(raw, bytes):
            return raw.decode('utf-8')
        return raw

    def _python_basic_statistics(self, data_files: Dict[str, Any]) -> Dict[str, Any]:
        """Fallback Python-based basic statistics"""
        results = {}
        
        for filename, file_info in data_files.items():
            try:
                # Read file content
                file_info['file'].seek(0)
                content = file_info['file'].read().decode('utf-8')
                lines = [line for line in content.split('\n') 
                        if line.strip() and not line.startswith('#')]
                
                stats = {}
                stats['total_regions'] = len(lines)
                
                if lines:
                    # Parse first line to determine format
                    first_line = lines[0].split('\t')
                    
                    # Extract scores if available (4th column)
                    scores = []
                    lengths = []
                    
                    for line in lines:
                        fields = line.split('\t')
                        if len(fields) >= 4:
                            try:
                                score = float(fields[3])
                                scores.append(score)
                            except ValueError:
                                pass
                        
                        if len(fields) >= 3:
                            try:
                                start = int(fields[1])
                                end = int(fields[2])
                                if end > start:
                                    lengths.append(end - start)
                            except ValueError:
                                pass
                    
                    # Score statistics
                    if scores:
                        stats['mean_score'] = sum(scores) / len(scores)
                        stats['median_score'] = sorted(scores)[len(scores)//2]
                        stats['max_score'] = max(scores)
                        stats['min_score'] = min(scores)
                    
                    # Length statistics
                    if lengths:
                        stats['mean_length'] = sum(lengths) / len(lengths)
                        stats['median_length'] = sorted(lengths)[len(lengths)//2]
                        stats['max_length'] = max(lengths)
                        stats['min_length'] = min(lengths)
                
                results[filename] = stats
                
            except Exception as e:
                results[filename] = {'error': f'Analysis failed: {str(e)}'}
        
        return results

    def _python_peak_calling(self, data_files: Dict[str, Any], **params) -> List[Dict[str, Any]]:
        """Python fallback for peak calling"""
        peaks = []
        threshold = params.get('fold_change', 2.0)
        
        for filename, file_info in data_files.items():
            try:
                file_info['file'].seek(0)
                content = file_info['file'].read().decode('utf-8')
                lines = [line for line in content.split('\n') 
                        if line.strip() and not line.startswith('#')]
                
                for line in lines:
                    fields = line.split('\t')
                    if len(fields) >= 4:
                        try:
                            chr_name = fields[0]
                            start = int(fields[1])
                            end = int(fields[2])
                            score = float(fields[3])
                            
                            if score > threshold:
                                peaks.append({
                                    'chr': chr_name,
                                    'start': start,
                                    'end': end,
                                    'score': score,
                                    'length': end - start,
                                    'file': filename
                                })
                        except (ValueError, IndexError):
                            continue
            
            except Exception as e:
                print(f"Error in peak calling for {filename}: {e}")
        
        return peaks
