"""
Unit tests for GEO client
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from utils.geo_client import GEOClient


@pytest.fixture
def geo_client():
    return GEOClient()


@pytest.fixture
def mock_search_response():
    return {
        'esearchresult': {
            'idlist': ['200012345', '200012346']
        }
    }


@pytest.fixture
def mock_summary_response():
    return {
        'result': {
            '200012345': {
                'accession': 'GSE12345',
                'title': 'Test Gene Expression Study',
                'summary': 'A test study of gene expression',
                'taxon': 'Homo sapiens',
                'gpl': 'GPL570',
                'n_samples': 10,
                'pdat': '2020-01-01',
                'entrytype': 'GSE'
            },
            '200012346': {
                'accession': 'GSE12346',
                'title': 'Another Test Study',
                'summary': 'Another test study',
                'taxon': 'Mus musculus',
                'gpl': 'GPL1261',
                'n_samples': 5,
                'pdat': '2020-02-01',
                'entrytype': 'GSE'
            }
        }
    }


def test_geo_client_initialization(geo_client):
    """Test GEO client initialization"""
    assert geo_client.base_url == "https://www.ncbi.nlm.nih.gov/geo"
    assert geo_client.eutils_base == "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@patch('utils.geo_client.requests.get')
def test_search_datasets(mock_get, geo_client, mock_search_response, mock_summary_response):
    """Test searching for GEO datasets"""
    # Mock the search and summary requests
    mock_search = Mock()
    mock_search.json.return_value = mock_search_response
    mock_search.raise_for_status = Mock()
    
    mock_summary = Mock()
    mock_summary.json.return_value = mock_summary_response
    mock_summary.raise_for_status = Mock()
    
    mock_get.side_effect = [mock_search, mock_summary]
    
    results = geo_client.search_datasets('cancer', dataset_type='gse', max_results=20)
    
    assert len(results) == 2
    assert results[0]['accession'] == 'GSE12345'
    assert results[0]['title'] == 'Test Gene Expression Study'
    assert results[0]['database'] == 'GEO'
    assert results[1]['organism'] == 'Mus musculus'


@patch('utils.geo_client.requests.get')
def test_search_datasets_no_results(mock_get, geo_client):
    """Test searching with no results"""
    mock_response = Mock()
    mock_response.json.return_value = {'esearchresult': {'idlist': []}}
    mock_response.raise_for_status = Mock()
    mock_get.return_value = mock_response
    
    results = geo_client.search_datasets('nonexistent')
    
    assert len(results) == 0


@patch('utils.geo_client.requests.get')
def test_search_datasets_connection_error(mock_get, geo_client):
    """Test handling of connection errors"""
    mock_get.side_effect = Exception("Connection error")
    
    results = geo_client.search_datasets('test')
    
    assert len(results) == 0


@patch('utils.geo_client.requests.get')
def test_download_series_matrix(mock_get, geo_client, tmp_path):
    """Test downloading series matrix file"""
    mock_response = Mock()
    mock_response.iter_content = Mock(return_value=[b'test content'])
    mock_response.raise_for_status = Mock()
    mock_get.return_value = mock_response
    
    download_path = geo_client.download_series_matrix('GSE12345', download_dir=str(tmp_path))
    
    assert download_path.endswith('GSE12345_series_matrix.txt.gz')
    assert 'GSE12345_series_matrix.txt.gz' in download_path


def test_parse_series_matrix(geo_client, tmp_path):
    """Test parsing series matrix file"""
    # Create a mock series matrix file
    matrix_content = """!Series_title\t"Test Series"
!Series_summary\t"Test summary"
!dataset_table_begin
ID_REF\tGSM1\tGSM2
GENE1\t10.5\t12.3
GENE2\t8.2\t9.1
!dataset_table_end
"""
    
    matrix_file = tmp_path / "test_series_matrix.txt"
    matrix_file.write_text(matrix_content)
    
    result = geo_client.parse_series_matrix(str(matrix_file))
    
    assert 'metadata' in result
    assert 'expression_data' in result
    assert result['metadata']['title'] == 'Test Series'
    assert result['expression_data']['n_rows'] == 2
    assert result['expression_data']['n_columns'] == 3


@patch('utils.geo_client.requests.get')
def test_get_dataset_details(mock_get, geo_client):
    """Test getting dataset details"""
    mock_response = Mock()
    mock_response.content = b'<xml><Item Name="title">Test Title</Item></xml>'
    mock_response.raise_for_status = Mock()
    mock_get.return_value = mock_response
    
    details = geo_client.get_dataset_details('GSE12345')
    
    assert details is not None
    assert 'accession' in details


def test_geo_base_url_format(geo_client):
    """Test that base URLs are properly formatted"""
    assert not geo_client.base_url.endswith('/')
    assert geo_client.eutils_base.startswith('https://')
    assert geo_client.ftp_base.startswith('https://')
