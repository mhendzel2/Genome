"""
DNA methylation analysis module for bisulfite sequencing data.
Supports methylation quantification, differential methylation, and region analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from scipy import stats
import logging


logger = logging.getLogger(__name__)


class MethylationAnalyzer:
    """
    Analyzer for DNA methylation data from bisulfite sequencing.
    """
    
    def __init__(self):
        self.supported_formats = ['bedGraph', 'BiSeq', 'Bismark', 'MethylKit']
    
    def load_bismark_coverage(self, file_path: str) -> pd.DataFrame:
        """
        Load Bismark coverage file.
        
        Format: chr, start, end, methylation%, count_methylated, count_unmethylated
        
        Args:
            file_path: Path to Bismark coverage file
        
        Returns:
            DataFrame with methylation data
        """
        try:
            df = pd.read_csv(file_path, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'methylation_pct', 
                                 'count_meth', 'count_unmeth'])
            
            # Calculate total coverage
            df['coverage'] = df['count_meth'] + df['count_unmeth']
            
            # Convert methylation percentage to fraction
            df['methylation_frac'] = df['methylation_pct'] / 100.0
            
            logger.info(f"Loaded Bismark coverage: {len(df)} CpG sites")
            return df
            
        except Exception as e:
            logger.exception("Error loading Bismark coverage: %s", e)
            raise ValueError(f"Failed to load Bismark data: {e}")
    
    def load_bedgraph_methylation(self, file_path: str) -> pd.DataFrame:
        """
        Load methylation data from bedGraph format.
        
        Args:
            file_path: Path to bedGraph file
        
        Returns:
            DataFrame with methylation data
        """
        try:
            df = pd.read_csv(file_path, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'methylation_value'])
            
            logger.info(f"Loaded bedGraph: {len(df)} methylation sites")
            return df
            
        except Exception as e:
            logger.exception("Error loading bedGraph: %s", e)
            raise ValueError(f"Failed to load bedGraph data: {e}")
    
    def filter_by_coverage(self, df: pd.DataFrame, 
                          min_coverage: int = 10,
                          max_coverage: Optional[int] = None) -> pd.DataFrame:
        """
        Filter methylation sites by coverage depth.
        
        Args:
            df: DataFrame with methylation data
            min_coverage: Minimum coverage threshold
            max_coverage: Maximum coverage threshold (to filter PCR duplicates)
        
        Returns:
            Filtered DataFrame
        """
        if 'coverage' not in df.columns:
            logger.warning("Coverage column not found, skipping filtering")
            return df
        
        df_filtered = df[df['coverage'] >= min_coverage].copy()
        
        if max_coverage is not None:
            df_filtered = df_filtered[df_filtered['coverage'] <= max_coverage].copy()
        
        logger.info(f"Filtered by coverage: {len(df_filtered)}/{len(df)} sites retained")
        return df_filtered
    
    def calculate_methylation_statistics(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Calculate global methylation statistics.
        
        Args:
            df: DataFrame with methylation data
        
        Returns:
            Dictionary with methylation statistics
        """
        meth_col = 'methylation_frac' if 'methylation_frac' in df.columns else 'methylation_value'
        
        if meth_col not in df.columns:
            return {'error': 'Methylation values not found'}
        
        methylation = df[meth_col].dropna()
        
        stats_dict = {
            'n_sites': len(df),
            'mean_methylation': float(methylation.mean()),
            'median_methylation': float(methylation.median()),
            'std_methylation': float(methylation.std()),
            'min_methylation': float(methylation.min()),
            'max_methylation': float(methylation.max()),
            'hypermethylated_sites': int((methylation > 0.75).sum()),  # > 75% methylated
            'hypomethylated_sites': int((methylation < 0.25).sum()),   # < 25% methylated
            'intermediate_sites': int(((methylation >= 0.25) & (methylation <= 0.75)).sum())
        }
        
        # Distribution by chromosome
        if 'chr' in df.columns:
            chr_stats = df.groupby('chr')[meth_col].agg(['mean', 'count'])
            stats_dict['methylation_by_chr'] = chr_stats.to_dict()
        
        return stats_dict
    
    def differential_methylation(self, control_files: List[pd.DataFrame],
                                treatment_files: List[pd.DataFrame],
                                min_diff: float = 0.25,
                                alpha: float = 0.05) -> pd.DataFrame:
        """
        Identify differentially methylated sites between two groups.
        
        Args:
            control_files: List of DataFrames for control samples
            treatment_files: List of DataFrames for treatment samples
            min_diff: Minimum methylation difference threshold
            alpha: Significance level
        
        Returns:
            DataFrame with differential methylation results
        """
        # Merge all samples on genomic coordinates
        # This is simplified - full implementation would handle replicates properly
        
        results = []
        
        # For each site, compare control vs treatment
        # This is a placeholder for the actual implementation
        logger.info("Performing differential methylation analysis...")
        
        # Would need to:
        # 1. Merge all dataframes on chr/start/end
        # 2. Calculate mean methylation for each group
        # 3. Perform statistical test (Fisher's exact test or chi-square)
        # 4. Apply multiple testing correction
        
        return pd.DataFrame(results)
    
    def identify_dmr(self, df: pd.DataFrame,
                    window_size: int = 1000,
                    min_sites: int = 3,
                    methylation_diff_threshold: float = 0.25) -> pd.DataFrame:
        """
        Identify Differentially Methylated Regions (DMRs).
        
        Args:
            df: DataFrame with differential methylation results
            window_size: Size of genomic window for DMR detection
            min_sites: Minimum number of CpG sites in a DMR
            methylation_diff_threshold: Minimum methylation difference
        
        Returns:
            DataFrame with DMR coordinates and statistics
        """
        if 'methylation_diff' not in df.columns:
            logger.warning("methylation_diff column not found")
            return pd.DataFrame()
        
        # Filter significant sites
        sig_sites = df[np.abs(df['methylation_diff']) >= methylation_diff_threshold].copy()
        
        # Sort by chromosome and position
        sig_sites = sig_sites.sort_values(['chr', 'start'])
        
        dmrs = []
        
        # Group by chromosome
        for chr_name, chr_df in sig_sites.groupby('chr'):
            # Find clusters of nearby sites
            chr_df = chr_df.reset_index(drop=True)
            
            current_dmr = []
            for idx, row in chr_df.iterrows():
                if not current_dmr:
                    current_dmr.append(row)
                else:
                    # Check if within window
                    if row['start'] - current_dmr[-1]['start'] <= window_size:
                        current_dmr.append(row)
                    else:
                        # Save current DMR if it has enough sites
                        if len(current_dmr) >= min_sites:
                            dmr_df = pd.DataFrame(current_dmr)
                            dmrs.append({
                                'chr': chr_name,
                                'start': int(dmr_df['start'].min()),
                                'end': int(dmr_df['end'].max()),
                                'n_sites': len(current_dmr),
                                'mean_diff': float(dmr_df['methylation_diff'].mean()),
                                'median_diff': float(dmr_df['methylation_diff'].median())
                            })
                        
                        # Start new DMR
                        current_dmr = [row]
            
            # Don't forget the last DMR
            if len(current_dmr) >= min_sites:
                dmr_df = pd.DataFrame(current_dmr)
                dmrs.append({
                    'chr': chr_name,
                    'start': int(dmr_df['start'].min()),
                    'end': int(dmr_df['end'].max()),
                    'n_sites': len(current_dmr),
                    'mean_diff': float(dmr_df['methylation_diff'].mean()),
                    'median_diff': float(dmr_df['methylation_diff'].median())
                })
        
        dmr_df = pd.DataFrame(dmrs)
        logger.info(f"Identified {len(dmr_df)} DMRs")
        return dmr_df
    
    def annotate_cpg_islands(self, df: pd.DataFrame,
                            cpg_islands_file: Optional[str] = None) -> pd.DataFrame:
        """
        Annotate methylation sites with CpG island information.
        
        Args:
            df: DataFrame with methylation data
            cpg_islands_file: Path to CpG islands BED file (optional)
        
        Returns:
            DataFrame with CpG island annotations
        """
        if cpg_islands_file is None:
            logger.warning("CpG islands file not provided, skipping annotation")
            df['cpg_island'] = 'unknown'
            return df
        
        try:
            # Load CpG islands
            cpg_islands = pd.read_csv(cpg_islands_file, sep='\t', header=None,
                                     names=['chr', 'start', 'end'])
            
            # Annotate each site
            df['cpg_island'] = 'intergenic'
            
            for idx, island in cpg_islands.iterrows():
                mask = (df['chr'] == island['chr']) & \
                       (df['start'] >= island['start']) & \
                       (df['end'] <= island['end'])
                df.loc[mask, 'cpg_island'] = 'island'
                
                # Shores (2kb flanking)
                shore_mask = (df['chr'] == island['chr']) & \
                            (((df['start'] >= island['start'] - 2000) & (df['start'] < island['start'])) |
                             ((df['end'] > island['end']) & (df['end'] <= island['end'] + 2000)))
                df.loc[shore_mask, 'cpg_island'] = 'shore'
            
            logger.info(f"Annotated {(df['cpg_island'] != 'intergenic').sum()} sites in CpG islands/shores")
            return df
            
        except Exception as e:
            logger.exception("Error annotating CpG islands: %s", e)
            df['cpg_island'] = 'error'
            return df
    
    def calculate_methylation_entropy(self, df: pd.DataFrame,
                                     window_size: int = 100) -> pd.DataFrame:
        """
        Calculate methylation entropy (epigenetic disorder) for genomic regions.
        
        Args:
            df: DataFrame with methylation data
            window_size: Window size for entropy calculation
        
        Returns:
            DataFrame with entropy values
        """
        meth_col = 'methylation_frac' if 'methylation_frac' in df.columns else 'methylation_value'
        
        if meth_col not in df.columns:
            logger.warning("Methylation values not found")
            return df
        
        df = df.sort_values(['chr', 'start']).reset_index(drop=True)
        
        entropy_values = []
        
        for chr_name, chr_df in df.groupby('chr'):
            chr_df = chr_df.reset_index(drop=True)
            
            for idx in range(len(chr_df)):
                # Get sites within window
                window_start = chr_df.loc[idx, 'start']
                window_end = window_start + window_size
                
                window_sites = chr_df[(chr_df['start'] >= window_start) & 
                                     (chr_df['start'] < window_end)][meth_col]
                
                if len(window_sites) > 0:
                    # Calculate entropy
                    probs = window_sites.value_counts(normalize=True)
                    entropy = -np.sum(probs * np.log2(probs + 1e-10))
                    entropy_values.append(entropy)
                else:
                    entropy_values.append(np.nan)
        
        df['methylation_entropy'] = entropy_values
        
        return df
    
    def compare_methylation_patterns(self, sample1: pd.DataFrame,
                                    sample2: pd.DataFrame,
                                    correlation_method: str = 'pearson') -> Dict[str, Any]:
        """
        Compare methylation patterns between two samples.
        
        Args:
            sample1: First sample methylation data
            sample2: Second sample methylation data
            correlation_method: Correlation method ('pearson' or 'spearman')
        
        Returns:
            Dictionary with comparison statistics
        """
        # Merge samples on genomic coordinates
        merged = pd.merge(sample1, sample2, on=['chr', 'start', 'end'],
                         suffixes=('_1', '_2'), how='inner')
        
        meth_col1 = 'methylation_frac_1' if 'methylation_frac_1' in merged.columns else 'methylation_value_1'
        meth_col2 = 'methylation_frac_2' if 'methylation_frac_2' in merged.columns else 'methylation_value_2'
        
        # Calculate correlation
        if correlation_method == 'pearson':
            corr, pval = stats.pearsonr(merged[meth_col1], merged[meth_col2])
        else:
            corr, pval = stats.spearmanr(merged[meth_col1], merged[meth_col2])
        
        # Calculate absolute differences
        diffs = np.abs(merged[meth_col1] - merged[meth_col2])
        
        return {
            'n_common_sites': len(merged),
            'correlation': float(corr),
            'correlation_pvalue': float(pval),
            'mean_absolute_difference': float(diffs.mean()),
            'median_absolute_difference': float(diffs.median()),
            'sites_with_large_difference': int((diffs > 0.25).sum())  # >25% difference
        }
