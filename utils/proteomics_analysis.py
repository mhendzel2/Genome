"""
Proteomics data analysis module for mass spectrometry data.
Supports protein identification, quantification, and differential expression.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from scipy import stats
import logging


logger = logging.getLogger(__name__)


class ProteomicsAnalyzer:
    """
    Analyzer for proteomics data from mass spectrometry experiments.
    """
    
    def __init__(self):
        self.supported_formats = ['mzTab', 'MaxQuant', 'Proteome Discoverer', 'Generic']
    
    def load_maxquant_data(self, file_path: str, 
                          intensity_columns: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Load MaxQuant proteinGroups.txt file.
        
        Args:
            file_path: Path to proteinGroups.txt
            intensity_columns: List of intensity column names to extract
        
        Returns:
            DataFrame with protein data
        """
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            # Filter out contaminants and reverse hits
            if 'Potential contaminant' in df.columns:
                df = df[df['Potential contaminant'] != '+']
            if 'Reverse' in df.columns:
                df = df[df['Reverse'] != '+']
            if 'Only identified by site' in df.columns:
                df = df[df['Only identified by site'] != '+']
            
            # Extract intensity columns if specified
            if intensity_columns:
                id_cols = ['Protein IDs', 'Gene names', 'Protein names']
                available_id_cols = [col for col in id_cols if col in df.columns]
                available_intensity_cols = [col for col in intensity_columns if col in df.columns]
                
                df = df[available_id_cols + available_intensity_cols]
            
            logger.info(f"Loaded MaxQuant data: {df.shape[0]} proteins, {df.shape[1]} columns")
            return df
            
        except Exception as e:
            logger.exception("Error loading MaxQuant data: %s", e)
            raise ValueError(f"Failed to load MaxQuant data: {e}")
    
    def normalize_intensities(self, df: pd.DataFrame, 
                             intensity_columns: List[str],
                             method: str = 'median') -> pd.DataFrame:
        """
        Normalize protein intensities.
        
        Args:
            df: DataFrame with protein intensities
            intensity_columns: List of column names containing intensities
            method: Normalization method ('median', 'quantile', 'log2')
        
        Returns:
            DataFrame with normalized intensities
        """
        df_norm = df.copy()
        
        for col in intensity_columns:
            if col not in df.columns:
                continue
            
            intensities = df[col].replace(0, np.nan)
            
            if method == 'median':
                # Median normalization
                median_val = intensities.median()
                global_median = df[intensity_columns].replace(0, np.nan).median().median()
                df_norm[col] = intensities * (global_median / median_val)
            
            elif method == 'quantile':
                # Quantile normalization (simplified)
                sorted_vals = np.sort(intensities.dropna())
                ranks = intensities.rank(method='average')
                df_norm[col] = ranks.map(lambda r: sorted_vals[int(r)-1] if not np.isnan(r) else np.nan)
            
            elif method == 'log2':
                # Log2 transformation
                df_norm[col] = np.log2(intensities)
            
            else:
                raise ValueError(f"Unknown normalization method: {method}")
        
        logger.info(f"Normalized {len(intensity_columns)} intensity columns using {method}")
        return df_norm
    
    def differential_expression_proteomics(self, df: pd.DataFrame,
                                          control_columns: List[str],
                                          treatment_columns: List[str],
                                          protein_id_col: str = 'Protein IDs',
                                          log_transform: bool = True,
                                          alpha: float = 0.05) -> pd.DataFrame:
        """
        Perform differential expression analysis on proteomics data.
        
        Args:
            df: DataFrame with protein intensities
            control_columns: Column names for control samples
            treatment_columns: Column names for treatment samples
            protein_id_col: Column name for protein identifiers
            log_transform: Whether to log2-transform intensities
            alpha: Significance threshold
        
        Returns:
            DataFrame with differential expression results
        """
        results = []
        
        for idx, row in df.iterrows():
            protein_id = row[protein_id_col]
            
            control_vals = row[control_columns].replace(0, np.nan).dropna().values
            treatment_vals = row[treatment_columns].replace(0, np.nan).dropna().values
            
            # Need at least 2 values in each group
            if len(control_vals) < 2 or len(treatment_vals) < 2:
                continue
            
            # Log transform if requested
            if log_transform:
                control_vals = np.log2(control_vals)
                treatment_vals = np.log2(treatment_vals)
            
            # Calculate fold change
            mean_control = np.mean(control_vals)
            mean_treatment = np.mean(treatment_vals)
            
            if log_transform:
                log2_fc = mean_treatment - mean_control
                fc = 2 ** log2_fc
            else:
                fc = mean_treatment / mean_control if mean_control != 0 else np.inf
                log2_fc = np.log2(fc)
            
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(treatment_vals, control_vals)
            
            results.append({
                'Protein ID': protein_id,
                'Gene': row.get('Gene names', ''),
                'Mean Control': mean_control,
                'Mean Treatment': mean_treatment,
                'Fold Change': fc,
                'Log2 Fold Change': log2_fc,
                't-statistic': t_stat,
                'p-value': p_value,
                'n_control': len(control_vals),
                'n_treatment': len(treatment_vals)
            })
        
        results_df = pd.DataFrame(results)
        
        # FDR correction
        from statsmodels.stats.multitest import fdrcorrection
        _, results_df['q-value'] = fdrcorrection(results_df['p-value'])
        
        # Add significance flag
        results_df['Significant'] = (results_df['q-value'] < alpha) & (np.abs(results_df['Log2 Fold Change']) > 1)
        
        # Sort by p-value
        results_df = results_df.sort_values('p-value')
        
        logger.info(f"DE analysis: {results_df['Significant'].sum()} significant proteins out of {len(results_df)}")
        return results_df
    
    def protein_abundance_analysis(self, df: pd.DataFrame,
                                   intensity_columns: List[str]) -> Dict[str, Any]:
        """
        Analyze protein abundance distribution.
        
        Args:
            df: DataFrame with protein intensities
            intensity_columns: Column names containing intensities
        
        Returns:
            Dictionary with abundance statistics
        """
        # Extract intensity data
        intensity_data = df[intensity_columns].replace(0, np.nan)
        
        # Calculate statistics
        stats_dict = {
            'n_proteins': len(df),
            'n_samples': len(intensity_columns),
            'mean_abundance': intensity_data.mean().mean(),
            'median_abundance': intensity_data.median().median(),
            'dynamic_range': np.log10(intensity_data.max().max() / intensity_data.min().min()),
            'missing_value_percentage': (intensity_data.isna().sum().sum() / intensity_data.size) * 100,
            'proteins_detected_per_sample': intensity_data.notna().sum().tolist(),
            'sample_correlations': intensity_data.corr().values.tolist()
        }
        
        return stats_dict
    
    def impute_missing_values(self, df: pd.DataFrame,
                             intensity_columns: List[str],
                             method: str = 'min') -> pd.DataFrame:
        """
        Impute missing values in proteomics data.
        
        Args:
            df: DataFrame with protein intensities
            intensity_columns: Column names containing intensities
            method: Imputation method ('min', 'mean', 'median', 'knn')
        
        Returns:
            DataFrame with imputed values
        """
        df_imputed = df.copy()
        
        for col in intensity_columns:
            intensities = df[col].replace(0, np.nan)
            
            if method == 'min':
                # Replace with minimum value (typical for left-censored data)
                min_val = intensities.min()
                df_imputed[col] = intensities.fillna(min_val * 0.1)  # 10% of minimum
            
            elif method == 'mean':
                df_imputed[col] = intensities.fillna(intensities.mean())
            
            elif method == 'median':
                df_imputed[col] = intensities.fillna(intensities.median())
            
            elif method == 'knn':
                # Simplified KNN imputation
                from sklearn.impute import KNNImputer
                imputer = KNNImputer(n_neighbors=5)
                df_imputed[intensity_columns] = imputer.fit_transform(df[intensity_columns].replace(0, np.nan))
                break  # KNN imputes all columns at once
            
            else:
                raise ValueError(f"Unknown imputation method: {method}")
        
        logger.info(f"Imputed missing values using {method} method")
        return df_imputed
    
    def identify_modified_peptides(self, df: pd.DataFrame,
                                   modification_col: str = 'Modifications') -> Dict[str, Any]:
        """
        Analyze post-translational modifications (PTMs).
        
        Args:
            df: DataFrame with peptide data
            modification_col: Column containing modification information
        
        Returns:
            Dictionary with PTM statistics
        """
        if modification_col not in df.columns:
            return {'error': 'Modification column not found'}
        
        # Count modification types
        all_mods = []
        for mods in df[modification_col].dropna():
            if isinstance(mods, str):
                all_mods.extend(mods.split(';'))
        
        mod_counts = pd.Series(all_mods).value_counts()
        
        return {
            'total_modified_proteins': df[modification_col].notna().sum(),
            'modification_types': mod_counts.to_dict(),
            'most_common_modification': mod_counts.index[0] if len(mod_counts) > 0 else None,
            'modification_frequency': len(all_mods) / len(df)
        }
