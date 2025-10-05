"""
Unit tests for statistical analysis module
"""

import pytest
import numpy as np
import pandas as pd
from utils.statistical_analysis import StatisticalAnalyzer


@pytest.fixture
def analyzer():
    return StatisticalAnalyzer()


@pytest.fixture
def sample_data():
    np.random.seed(42)
    return {
        'group1': np.random.normal(10, 2, 30),
        'group2': np.random.normal(12, 2, 30),
        'group3': np.random.normal(11, 2, 30),
        'matrix': np.random.randn(50, 10)
    }


def test_t_test_analysis(analyzer, sample_data):
    """Test t-test analysis"""
    result = analyzer.t_test_analysis(sample_data['group1'], sample_data['group2'])
    
    assert 't_statistic' in result
    assert 'p_value' in result
    assert 'cohens_d' in result
    assert isinstance(result['p_value'], float)
    assert 0 <= result['p_value'] <= 1


def test_t_test_paired(analyzer, sample_data):
    """Test paired t-test"""
    result = analyzer.t_test_analysis(
        sample_data['group1'], 
        sample_data['group2'], 
        paired=True
    )
    
    assert 't_statistic' in result
    assert 'p_value' in result
    assert result['mean_group1'] != result['mean_group2']


def test_anova_analysis(analyzer, sample_data):
    """Test ANOVA"""
    result = analyzer.anova_analysis(
        sample_data['group1'],
        sample_data['group2'],
        sample_data['group3']
    )
    
    assert 'f_statistic' in result
    assert 'p_value' in result
    assert 'eta_squared' in result
    assert result['num_groups'] == 3
    assert len(result['group_means']) == 3


def test_multiple_testing_correction(analyzer):
    """Test multiple testing correction"""
    pvalues = [0.001, 0.01, 0.05, 0.1, 0.5]
    
    result = analyzer.multiple_testing_correction(pvalues, method='fdr_bh', alpha=0.05)
    
    assert 'corrected_pvalues' in result
    assert 'rejected' in result
    assert len(result['corrected_pvalues']) == len(pvalues)
    assert len(result['rejected']) == len(pvalues)
    assert isinstance(result['num_significant'], int)


def test_multiple_testing_with_nan(analyzer):
    """Test multiple testing correction with NaN values"""
    pvalues = [0.001, np.nan, 0.05, 0.1, np.nan]
    
    result = analyzer.multiple_testing_correction(pvalues, method='bonferroni')
    
    assert len(result['corrected_pvalues']) == len(pvalues)
    assert result['num_tests'] == 3  # Only valid values


def test_fold_change_analysis(analyzer, sample_data):
    """Test fold change calculation"""
    result = analyzer.fold_change_analysis(
        sample_data['group1'],
        sample_data['group2'],
        log_transform=True
    )
    
    assert 'fold_change' in result
    assert 'log2_fold_change' in result
    assert 'direction' in result
    assert result['direction'] in ['up', 'down']


def test_correlation_analysis(analyzer, sample_data):
    """Test correlation analysis"""
    x = sample_data['group1']
    y = sample_data['group2']
    
    # Pearson
    result = analyzer.correlation_analysis(x, y, method='pearson')
    assert 'correlation' in result
    assert 'p_value' in result
    assert -1 <= result['correlation'] <= 1
    
    # Spearman
    result = analyzer.correlation_analysis(x, y, method='spearman')
    assert result['method'] == 'spearman'
    
    # Kendall
    result = analyzer.correlation_analysis(x, y, method='kendall')
    assert result['method'] == 'kendall'


def test_hierarchical_clustering(analyzer, sample_data):
    """Test hierarchical clustering"""
    result = analyzer.hierarchical_clustering(
        sample_data['matrix'],
        n_clusters=3
    )
    
    assert 'linkage_matrix' in result
    assert 'cluster_labels' in result
    assert 'n_clusters' in result
    assert len(result['cluster_labels']) == sample_data['matrix'].shape[0]


def test_kmeans_clustering(analyzer, sample_data):
    """Test K-means clustering"""
    result = analyzer.kmeans_clustering(
        sample_data['matrix'],
        n_clusters=3
    )
    
    assert 'cluster_labels' in result
    assert 'cluster_centers' in result
    assert 'inertia' in result
    assert len(result['cluster_labels']) == sample_data['matrix'].shape[0]
    assert result['n_clusters'] == 3


def test_pca_analysis(analyzer, sample_data):
    """Test PCA analysis"""
    result = analyzer.pca_analysis(sample_data['matrix'], n_components=3)
    
    assert 'transformed_data' in result
    assert 'explained_variance_ratio' in result
    assert 'cumulative_variance_ratio' in result
    assert result['n_components'] == 3
    assert len(result['transformed_data']) == sample_data['matrix'].shape[0]


def test_mann_whitney_u_test(analyzer, sample_data):
    """Test Mann-Whitney U test"""
    result = analyzer.mann_whitney_u_test(
        sample_data['group1'],
        sample_data['group2']
    )
    
    assert 'u_statistic' in result
    assert 'p_value' in result
    assert 'median_group1' in result
    assert 'median_group2' in result


def test_kruskal_wallis_test(analyzer, sample_data):
    """Test Kruskal-Wallis H test"""
    result = analyzer.kruskal_wallis_test(
        sample_data['group1'],
        sample_data['group2'],
        sample_data['group3']
    )
    
    assert 'h_statistic' in result
    assert 'p_value' in result
    assert result['num_groups'] == 3


def test_wilcoxon_signed_rank_test(analyzer, sample_data):
    """Test Wilcoxon signed-rank test"""
    result = analyzer.wilcoxon_signed_rank_test(
        sample_data['group1'],
        sample_data['group2']
    )
    
    assert 'statistic' in result
    assert 'p_value' in result
    assert 'median_difference' in result


def test_chi_square_test(analyzer):
    """Test chi-square test"""
    contingency_table = np.array([[10, 20], [15, 25]])
    
    result = analyzer.chi_square_test(contingency_table)
    
    assert 'chi_square' in result
    assert 'p_value' in result
    assert 'degrees_of_freedom' in result
    assert 'expected_frequencies' in result


def test_normality_test(analyzer, sample_data):
    """Test normality tests"""
    # Shapiro-Wilk
    result = analyzer.normality_test(sample_data['group1'], method='shapiro')
    assert result['test'] == 'Shapiro-Wilk'
    assert 'is_normal' in result
    
    # Kolmogorov-Smirnov
    result = analyzer.normality_test(sample_data['group1'], method='kstest')
    assert result['test'] == 'Kolmogorov-Smirnov'


def test_variance_test(analyzer, sample_data):
    """Test variance equality test"""
    result = analyzer.variance_test(
        sample_data['group1'],
        sample_data['group2']
    )
    
    assert 'statistic' in result
    assert 'p_value' in result
    assert 'equal_variances' in result
    assert 'variance_group1' in result
    assert 'variance_group2' in result


def test_survival_analysis_logrank(analyzer):
    """Test log-rank test for survival analysis"""
    np.random.seed(42)
    times1 = np.random.exponential(100, 50)
    events1 = np.random.binomial(1, 0.7, 50)
    times2 = np.random.exponential(120, 50)
    events2 = np.random.binomial(1, 0.6, 50)
    
    result = analyzer.survival_analysis_logrank(times1, events1, times2, events2)
    
    assert 'chi_square' in result
    assert 'p_value' in result
    assert 'observed_events' in result
    assert 'expected_events' in result


def test_invalid_inputs(analyzer):
    """Test error handling with invalid inputs"""
    with pytest.raises(ValueError):
        analyzer.t_test_analysis(
            np.array([1, 2]),
            np.array([1, 2, 3]),
            paired=True
        )
    
    with pytest.raises(ValueError):
        analyzer.anova_analysis(np.array([1, 2, 3]))  # Only one group
    
    with pytest.raises(ValueError):
        analyzer.correlation_analysis(
            np.array([1, 2, 3]),
            np.array([1, 2, 3]),
            method='invalid'
        )
