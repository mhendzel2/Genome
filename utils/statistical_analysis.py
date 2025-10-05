"""
Advanced statistical analysis methods for genomics data.
Provides robust statistical tests, multiple testing corrections, and survival analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from statsmodels.stats.multitest import multipletests, fdrcorrection
import warnings


class StatisticalAnalyzer:
    """Comprehensive statistical analysis for genomics data"""
    
    def __init__(self):
        self.scaler = StandardScaler()
    
    def t_test_analysis(self, group1: np.ndarray, group2: np.ndarray, 
                       paired: bool = False, alternative: str = 'two-sided') -> Dict[str, Any]:
        """
        Perform t-test between two groups.
        
        Args:
            group1: First group of samples
            group2: Second group of samples
            paired: Whether to perform paired t-test
            alternative: 'two-sided', 'less', or 'greater'
        
        Returns:
            Dictionary with t-statistic, p-value, and effect size
        """
        if paired:
            if len(group1) != len(group2):
                raise ValueError("Groups must have equal length for paired t-test")
            statistic, pvalue = stats.ttest_rel(group1, group2, alternative=alternative)
        else:
            statistic, pvalue = stats.ttest_ind(group1, group2, alternative=alternative)
        
        # Calculate Cohen's d (effect size)
        pooled_std = np.sqrt((np.var(group1) + np.var(group2)) / 2)
        cohens_d = (np.mean(group1) - np.mean(group2)) / pooled_std if pooled_std > 0 else 0
        
        return {
            't_statistic': float(statistic),
            'p_value': float(pvalue),
            'cohens_d': float(cohens_d),
            'mean_group1': float(np.mean(group1)),
            'mean_group2': float(np.mean(group2)),
            'std_group1': float(np.std(group1)),
            'std_group2': float(np.std(group2))
        }
    
    def anova_analysis(self, *groups: np.ndarray) -> Dict[str, Any]:
        """
        Perform one-way ANOVA across multiple groups.
        
        Args:
            *groups: Variable number of sample groups
        
        Returns:
            Dictionary with F-statistic and p-value
        """
        if len(groups) < 2:
            raise ValueError("At least 2 groups required for ANOVA")
        
        f_statistic, pvalue = stats.f_oneway(*groups)
        
        # Calculate effect size (eta-squared)
        grand_mean = np.mean(np.concatenate(groups))
        ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
        ss_total = sum(np.sum((g - grand_mean)**2) for g in groups)
        eta_squared = ss_between / ss_total if ss_total > 0 else 0
        
        return {
            'f_statistic': float(f_statistic),
            'p_value': float(pvalue),
            'eta_squared': float(eta_squared),
            'num_groups': len(groups),
            'group_means': [float(np.mean(g)) for g in groups]
        }
    
    def multiple_testing_correction(self, pvalues: List[float], 
                                   method: str = 'fdr_bh', 
                                   alpha: float = 0.05) -> Dict[str, Any]:
        """
        Apply multiple testing correction to p-values.
        
        Args:
            pvalues: List of p-values
            method: Correction method ('bonferroni', 'fdr_bh', 'fdr_by', 'holm', etc.)
            alpha: Significance level
        
        Returns:
            Dictionary with corrected p-values and rejected hypotheses
        """
        pvalues_array = np.array(pvalues)
        
        # Remove NaN values and track their positions
        valid_mask = ~np.isnan(pvalues_array)
        valid_pvalues = pvalues_array[valid_mask]
        
        if len(valid_pvalues) == 0:
            return {
                'method': method,
                'corrected_pvalues': pvalues_array.tolist(),
                'rejected': [False] * len(pvalues),
                'num_significant': 0
            }
        
        # Apply correction
        rejected, corrected_pvals, _, _ = multipletests(valid_pvalues, alpha=alpha, method=method)
        
        # Reconstruct full arrays with NaN in original positions
        full_corrected = np.full_like(pvalues_array, np.nan)
        full_rejected = np.full(len(pvalues_array), False)
        full_corrected[valid_mask] = corrected_pvals
        full_rejected[valid_mask] = rejected
        
        return {
            'method': method,
            'alpha': alpha,
            'corrected_pvalues': full_corrected.tolist(),
            'rejected': full_rejected.tolist(),
            'num_significant': int(np.sum(rejected)),
            'num_tests': len(valid_pvalues),
            'proportion_significant': float(np.sum(rejected) / len(valid_pvalues))
        }
    
    def fold_change_analysis(self, group1: np.ndarray, group2: np.ndarray, 
                           log_transform: bool = True) -> Dict[str, Any]:
        """
        Calculate fold change between two groups.
        
        Args:
            group1: First group (e.g., control)
            group2: Second group (e.g., treatment)
            log_transform: Whether to return log2 fold change
        
        Returns:
            Dictionary with fold change metrics
        """
        mean1 = np.mean(group1)
        mean2 = np.mean(group2)
        
        # Avoid division by zero
        if mean1 == 0:
            mean1 = np.finfo(float).eps
        
        fc = mean2 / mean1
        log2fc = np.log2(fc) if log_transform else fc
        
        return {
            'fold_change': float(fc),
            'log2_fold_change': float(log2fc),
            'mean_group1': float(mean1),
            'mean_group2': float(mean2),
            'direction': 'up' if fc > 1 else 'down'
        }
    
    def correlation_analysis(self, x: np.ndarray, y: np.ndarray, 
                           method: str = 'pearson') -> Dict[str, Any]:
        """
        Calculate correlation between two variables.
        
        Args:
            x: First variable
            y: Second variable
            method: 'pearson', 'spearman', or 'kendall'
        
        Returns:
            Dictionary with correlation coefficient and p-value
        """
        if method == 'pearson':
            corr, pvalue = stats.pearsonr(x, y)
        elif method == 'spearman':
            corr, pvalue = stats.spearmanr(x, y)
        elif method == 'kendall':
            corr, pvalue = stats.kendalltau(x, y)
        else:
            raise ValueError(f"Unknown correlation method: {method}")
        
        return {
            'method': method,
            'correlation': float(corr),
            'p_value': float(pvalue),
            'n_samples': len(x)
        }
    
    def hierarchical_clustering(self, data: np.ndarray, 
                              method: str = 'ward', 
                              metric: str = 'euclidean',
                              n_clusters: Optional[int] = None) -> Dict[str, Any]:
        """
        Perform hierarchical clustering on data.
        
        Args:
            data: Data matrix (samples x features)
            method: Linkage method ('ward', 'complete', 'average', etc.)
            metric: Distance metric
            n_clusters: Number of clusters to extract (optional)
        
        Returns:
            Dictionary with linkage matrix and cluster assignments
        """
        # Standardize data
        data_scaled = self.scaler.fit_transform(data)
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(data_scaled, method=method, metric=metric)
        
        result = {
            'method': method,
            'metric': metric,
            'linkage_matrix': linkage_matrix.tolist()
        }
        
        # Extract clusters if requested
        if n_clusters is not None:
            clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            result['cluster_labels'] = clusters.tolist()
            result['n_clusters'] = n_clusters
        
        return result
    
    def kmeans_clustering(self, data: np.ndarray, 
                         n_clusters: int = 3,
                         n_init: int = 10,
                         random_state: int = 42) -> Dict[str, Any]:
        """
        Perform K-means clustering.
        
        Args:
            data: Data matrix (samples x features)
            n_clusters: Number of clusters
            n_init: Number of initializations
            random_state: Random seed
        
        Returns:
            Dictionary with cluster assignments and metrics
        """
        # Standardize data
        data_scaled = self.scaler.fit_transform(data)
        
        # Perform K-means
        kmeans = KMeans(n_clusters=n_clusters, n_init=n_init, random_state=random_state)
        cluster_labels = kmeans.fit_predict(data_scaled)
        
        return {
            'cluster_labels': cluster_labels.tolist(),
            'n_clusters': n_clusters,
            'cluster_centers': kmeans.cluster_centers_.tolist(),
            'inertia': float(kmeans.inertia_),
            'n_iterations': int(kmeans.n_iter_)
        }
    
    def pca_analysis(self, data: np.ndarray, 
                    n_components: Optional[int] = None) -> Dict[str, Any]:
        """
        Perform Principal Component Analysis.
        
        Args:
            data: Data matrix (samples x features)
            n_components: Number of components (default: min(n_samples, n_features))
        
        Returns:
            Dictionary with PCA results
        """
        # Standardize data
        data_scaled = self.scaler.fit_transform(data)
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        transformed = pca.fit_transform(data_scaled)
        
        return {
            'transformed_data': transformed.tolist(),
            'explained_variance': pca.explained_variance_.tolist(),
            'explained_variance_ratio': pca.explained_variance_ratio_.tolist(),
            'cumulative_variance_ratio': np.cumsum(pca.explained_variance_ratio_).tolist(),
            'n_components': pca.n_components_,
            'components': pca.components_.tolist()
        }
    
    def survival_analysis_logrank(self, times1: np.ndarray, events1: np.ndarray,
                                 times2: np.ndarray, events2: np.ndarray) -> Dict[str, Any]:
        """
        Perform log-rank test for survival analysis.
        
        Args:
            times1: Survival times for group 1
            events1: Event indicators for group 1 (1=event, 0=censored)
            times2: Survival times for group 2
            events2: Event indicators for group 2
        
        Returns:
            Dictionary with test statistic and p-value
        """
        # Combine data
        all_times = np.concatenate([times1, times2])
        all_events = np.concatenate([events1, events2])
        all_groups = np.concatenate([np.zeros(len(times1)), np.ones(len(times2))])
        
        # Get unique event times
        unique_times = np.unique(all_times[all_events == 1])
        
        # Calculate log-rank statistic
        observed = np.zeros(2)
        expected = np.zeros(2)
        
        for t in unique_times:
            # At risk in each group
            at_risk = [np.sum((all_times >= t) & (all_groups == g)) for g in [0, 1]]
            
            # Events in each group
            events = [np.sum((all_times == t) & (all_events == 1) & (all_groups == g)) for g in [0, 1]]
            
            total_at_risk = sum(at_risk)
            total_events = sum(events)
            
            if total_at_risk > 0:
                for g in [0, 1]:
                    observed[g] += events[g]
                    expected[g] += (at_risk[g] / total_at_risk) * total_events
        
        # Calculate chi-square statistic
        chi_square = np.sum((observed - expected)**2 / (expected + 1e-10))
        pvalue = 1 - stats.chi2.cdf(chi_square, df=1)
        
        return {
            'chi_square': float(chi_square),
            'p_value': float(pvalue),
            'observed_events': observed.tolist(),
            'expected_events': expected.tolist()
        }
    
    def mann_whitney_u_test(self, group1: np.ndarray, group2: np.ndarray,
                           alternative: str = 'two-sided') -> Dict[str, Any]:
        """
        Perform Mann-Whitney U test (non-parametric alternative to t-test).
        
        Args:
            group1: First group of samples
            group2: Second group of samples
            alternative: 'two-sided', 'less', or 'greater'
        
        Returns:
            Dictionary with U statistic and p-value
        """
        statistic, pvalue = stats.mannwhitneyu(group1, group2, alternative=alternative)
        
        return {
            'u_statistic': float(statistic),
            'p_value': float(pvalue),
            'median_group1': float(np.median(group1)),
            'median_group2': float(np.median(group2)),
            'n_group1': len(group1),
            'n_group2': len(group2)
        }
    
    def kruskal_wallis_test(self, *groups: np.ndarray) -> Dict[str, Any]:
        """
        Perform Kruskal-Wallis H test (non-parametric alternative to ANOVA).
        
        Args:
            *groups: Variable number of sample groups
        
        Returns:
            Dictionary with H statistic and p-value
        """
        if len(groups) < 2:
            raise ValueError("At least 2 groups required for Kruskal-Wallis test")
        
        statistic, pvalue = stats.kruskal(*groups)
        
        return {
            'h_statistic': float(statistic),
            'p_value': float(pvalue),
            'num_groups': len(groups),
            'group_medians': [float(np.median(g)) for g in groups]
        }
    
    def wilcoxon_signed_rank_test(self, group1: np.ndarray, group2: np.ndarray,
                                  alternative: str = 'two-sided') -> Dict[str, Any]:
        """
        Perform Wilcoxon signed-rank test (paired non-parametric test).
        
        Args:
            group1: First group of paired samples
            group2: Second group of paired samples
            alternative: 'two-sided', 'less', or 'greater'
        
        Returns:
            Dictionary with test statistic and p-value
        """
        if len(group1) != len(group2):
            raise ValueError("Groups must have equal length for paired test")
        
        statistic, pvalue = stats.wilcoxon(group1, group2, alternative=alternative)
        
        return {
            'statistic': float(statistic),
            'p_value': float(pvalue),
            'median_difference': float(np.median(group1 - group2)),
            'n_pairs': len(group1)
        }
    
    def chi_square_test(self, contingency_table: np.ndarray) -> Dict[str, Any]:
        """
        Perform chi-square test of independence.
        
        Args:
            contingency_table: 2D contingency table
        
        Returns:
            Dictionary with test results
        """
        chi2, pvalue, dof, expected = stats.chi2_contingency(contingency_table)
        
        return {
            'chi_square': float(chi2),
            'p_value': float(pvalue),
            'degrees_of_freedom': int(dof),
            'expected_frequencies': expected.tolist()
        }
    
    def normality_test(self, data: np.ndarray, method: str = 'shapiro') -> Dict[str, Any]:
        """
        Test for normality of data distribution.
        
        Args:
            data: Data array
            method: 'shapiro' or 'kstest'
        
        Returns:
            Dictionary with test results
        """
        if method == 'shapiro':
            statistic, pvalue = stats.shapiro(data)
            test_name = 'Shapiro-Wilk'
        elif method == 'kstest':
            statistic, pvalue = stats.kstest(data, 'norm')
            test_name = 'Kolmogorov-Smirnov'
        else:
            raise ValueError(f"Unknown normality test method: {method}")
        
        return {
            'test': test_name,
            'statistic': float(statistic),
            'p_value': float(pvalue),
            'is_normal': pvalue > 0.05,
            'interpretation': 'Data appears normally distributed' if pvalue > 0.05 else 'Data does not appear normally distributed'
        }
    
    def variance_test(self, group1: np.ndarray, group2: np.ndarray) -> Dict[str, Any]:
        """
        Test for equality of variances (Levene's test).
        
        Args:
            group1: First group of samples
            group2: Second group of samples
        
        Returns:
            Dictionary with test results
        """
        statistic, pvalue = stats.levene(group1, group2)
        
        return {
            'test': "Levene's test",
            'statistic': float(statistic),
            'p_value': float(pvalue),
            'equal_variances': pvalue > 0.05,
            'variance_group1': float(np.var(group1)),
            'variance_group2': float(np.var(group2))
        }
