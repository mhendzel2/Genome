import os
import subprocess
import tempfile
import json
from typing import Dict, Any, List, Optional
import pandas as pd

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

    def _execute_r_script(self, script: str, data_files: Optional[Dict[str, str]] = None) -> str:
        """Execute R script safely with input validation"""
        if not self.r_available:
            raise RuntimeError("R is not available on this system")

        if not self._validate_r_script(script):
            raise ValueError("R script contains potentially unsafe content")

        import logging
        logger = logging.getLogger(__name__)

        with tempfile.TemporaryDirectory() as temp_dir:
            # Write data files to temporary paths and prepare replacements
            file_map = {}
            if data_files:
                for filename, content in data_files.items():
                    safe_filename = self._sanitize_filename(filename)
                    file_path = os.path.join(temp_dir, safe_filename)
                    with open(file_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                    file_map[filename] = file_path

            # Allow scripts to reference input files via a simple placeholder syntax {file:<name>}
            # Replace placeholders in the script with the temp file paths where possible.
            for original_name, path in file_map.items():
                placeholder = f"{{file:{original_name}}}"
                script = script.replace(placeholder, path.replace('\\', '/'))

            # Write R script to temporary file with safe filename
            script_path = os.path.join(temp_dir, 'analysis.R')
            with open(script_path, 'w', encoding='utf-8') as f:
                f.write(script)

            # Execute R script with restricted environment
            try:
                result = subprocess.run(
                    ['R', '--vanilla', '--slave', '-f', script_path],
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
            # Fallback to Python-based statistics
            return self._python_basic_statistics(data_files)

        results = {}

        for filename, file_info in data_files.items():
            try:
                # Prepare file content mapping; support file-like objects
                if hasattr(file_info['file'], 'getvalue'):
                    file_content = file_info['file'].getvalue()
                else:
                    try:
                        file_info['file'].seek(0)
                    except Exception:
                        pass
                    raw = file_info['file'].read()
                    if isinstance(raw, bytes):
                        file_content = raw.decode('utf-8')
                    else:
                        file_content = raw

                data_files_payload = {filename: file_content}

                # Create R script that reads from a file path placeholder
                r_script = '''
                library(jsonlite)
                # Read data from provided file path (placeholder will be replaced)
                data <- read.table(file = "{file:%s}", sep = "\t", header = FALSE, comment.char = "#")

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
                ''' % (filename)

                # Execute R script, letting _execute_r_script write the file content to temp
                output = self._execute_r_script(r_script, data_files=data_files_payload)

                # Parse JSON output
                try:
                    import json as _json
                    stats = _json.loads(output.strip())
                    results[filename] = stats
                except Exception:
                    # Fallback to basic parsing
                    results[filename] = {'total_regions': 'N/A', 'error': 'Could not parse R output'}

            except Exception as e:
                results[filename] = {'error': f'Analysis failed: {str(e)}'}

        return results

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

    def peak_calling(self, data_files: Dict[str, Any], **params) -> List[Dict[str, Any]]:
        """Perform peak calling analysis"""
        if not self.r_available:
            return self._python_peak_calling(data_files, **params)

        import logging
        logger = logging.getLogger(__name__)

        results = []
        threshold = float(params.get('fold_change', 2.0))

        for filename, file_info in data_files.items():
            if 'chip' in filename.lower() or file_info.get('type') == 'ChIP-seq':
                try:
                    try:
                        file_info['file'].seek(0)
                    except Exception:
                        pass
                    raw = file_info['file'].read()
                    content = raw.decode('utf-8') if isinstance(raw, bytes) else raw

                    payload = {filename: content}

                    r_script = '''
                    library(jsonlite)
                    dat <- read.table(file = "{file:%s}", sep = "\t", header = FALSE, comment.char = "#")
                    peaks <- list()
                    if (ncol(dat) >= 4) {
                        for (i in seq_len(min(nrow(dat), 100))) {
                            row <- dat[i, ]
                            score <- as.numeric(row[[4]])
                            if (!is.na(score) && score > %f) {
                                peaks[[length(peaks) + 1]] <- list(chr = as.character(row[[1]]), start = as.integer(row[[2]]), end = as.integer(row[[3]]), score = score)
                            }
                        }
                    }
                    cat(toJSON(peaks, auto_unbox = TRUE))
                    ''' % (filename, threshold)

                    out = self._execute_r_script(r_script, data_files=payload)
                    import json as _json
                    parsed = _json.loads(out.strip()) if out and out.strip() else []
                    for p in parsed:
                        p['file'] = filename
                        if 'length' not in p and 'start' in p and 'end' in p:
                            p['length'] = p.get('end', 0) - p.get('start', 0)
                        results.append(p)
                except Exception as e:
                    logger.exception('Error calling peaks for %s: %s', filename, e)

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

    def differential_expression(self, data_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """Perform differential expression analysis"""
        import logging
        logger = logging.getLogger(__name__)

        if not self.r_available:
            return self._python_basic_statistics(data_files)

        """Perform differential expression analysis and return the full results table."""
        import logging
        logger = logging.getLogger(__name__)

        # This implementation assumes we process one file for simplicity.
        # A more robust implementation would handle control vs. treatment groups.
        first_file_key = next(iter(data_files), None)
        if not first_file_key:
            return pd.DataFrame()

        file_info = data_files[first_file_key]

        try:
            file_info['file'].seek(0)
            content = file_info['file'].read().decode('utf-8')
            payload = {first_file_key: content}

            # This R script now reads a table and outputs it as JSON
            r_script = '''
            library(jsonlite)
            # Assuming the input file is a standard DE result table
            dat <- read.table(file = "{file:%s}", sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)

            # Check for required columns for a volcano plot
            if (!("log2FoldChange" %in% colnames(dat)) || !("pvalue" %in% colnames(dat))) {
                stop("Input file must contain 'log2FoldChange' and 'pvalue' columns.")
            }

            # Output the entire dataframe as JSON
            cat(toJSON(dat))
            ''' % (first_file_key)

            json_output = self._execute_r_script(r_script, data_files=payload)

            # The output is a JSON string, which will be parsed into a DataFrame
            # in the calling function (in GenomicsAnalyzer)
            return json_output

        except Exception as e:
            logger.exception('Differential expression analysis failed for %s: %s', first_file_key, e)
            return pd.DataFrame()

    def tissue_comparison(self, tissues: List[str], comparison_type: str, data_files: Dict[str, Any]) -> Dict[str, Any]:
        """Perform tissue comparison analysis"""
        import logging
        logger = logging.getLogger(__name__)

        results = {'tissues_compared': tissues, 'comparison_type': comparison_type, 'common_features': 0, 'unique_tissue1': 0, 'unique_tissue2': 0, 'total_features': 0}

        try:
            payload = {}
            for filename, file_info in data_files.items():
                try:
                    try:
                        file_info['file'].seek(0)
                    except Exception:
                        pass
                    raw = file_info['file'].read()
                    content = raw.decode('utf-8') if isinstance(raw, bytes) else raw
                    payload[filename] = content
                except Exception as e:
                    logger.exception('Could not read %s: %s', filename, e)

            if len(payload) >= 2:
                names = list(payload.keys())
                f1 = names[0]
                f2 = names[1]
                r_script = '''
                    library(jsonlite)
                    a <- read.table(file = "{file:%s}", sep = "\t", header = FALSE, comment.char = "#")
                    b <- read.table(file = "{file:%s}", sep = "\t", header = FALSE, comment.char = "#")
                    fa <- apply(a[,1:3], 1, function(x) paste0(x[1], ':', x[2], '-', x[3]))
                    fb <- apply(b[,1:3], 1, function(x) paste0(x[1], ':', x[2], '-', x[3]))
                    common <- length(intersect(fa, fb))
                    res <- list(common_features = common, total_a = length(fa), total_b = length(fb))
                    cat(toJSON(res, auto_unbox = TRUE))
                ''' % (f1, f2)

                out = self._execute_r_script(r_script, data_files=payload)
                import json as _json
                parsed = _json.loads(out.strip())
                results.update({'common_features': parsed.get('common_features', 0), 'total_features': parsed.get('total_a', 0) + parsed.get('total_b', 0)})
        except Exception as e:
            logger.exception('Tissue comparison failed: %s', e)

        return results
