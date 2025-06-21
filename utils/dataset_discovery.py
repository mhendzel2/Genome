import random
from typing import Dict, List, Any, Optional

class DatasetDiscovery:
    """Mock dataset discovery service for genomics databases"""
    
    def __init__(self):
        # Mock database of public datasets
        self.mock_datasets = self._generate_mock_datasets()
    
    def _generate_mock_datasets(self) -> List[Dict[str, Any]]:
        """Generate mock genomics datasets"""
        datasets = []
        
        # Define realistic dataset templates
        data_types = ["Histone Marks", "Gene Expression", "ChIP-seq", "HiC"]
        tissues = ["Liver", "Heart", "Brain", "Lung", "Kidney", "Muscle", "Skin", "Blood", "Spleen", "Pancreas", "Intestine", "Stomach", "Thymus", "Bone Marrow"]
        cell_types = ["Any", "Primary cells", "Cell lines", "Stem cells", "Immune cells"]
        organisms = ["Human", "Mouse"]
        databases = ["ENCODE", "GEO", "ArrayExpress", "TCGA"]
        
        # Histone mark datasets
        histone_marks = ["H3K4me1", "H3K4me3", "H3K27me3", "H3K36me3", "H3K9me3", "H3K27ac"]
        
        # Generate datasets
        for i in range(100):  # Generate 100 mock datasets
            organism = random.choice(organisms)
            tissue = random.choice(tissues)
            data_type = random.choice(data_types)
            database = random.choice(databases)
            cell_type = random.choice(cell_types)
            
            # Generate titles and descriptions based on data type
            if data_type == "Histone Marks":
                mark = random.choice(histone_marks)
                title = f"{mark} ChIP-seq in {organism} {tissue}"
                description = f"Chromatin immunoprecipitation sequencing of {mark} histone modification in {tissue.lower()} tissue from {organism.lower()}."
            elif data_type == "Gene Expression":
                title = f"RNA-seq in {organism} {tissue}"
                description = f"RNA sequencing analysis of gene expression in {tissue.lower()} tissue from {organism.lower()}."
            elif data_type == "ChIP-seq":
                proteins = ["CTCF", "p53", "MYC", "STAT1", "NFkB", "RNA Pol II"]
                protein = random.choice(proteins)
                title = f"{protein} ChIP-seq in {organism} {tissue}"
                description = f"Chromatin immunoprecipitation sequencing of {protein} binding sites in {tissue.lower()} tissue from {organism.lower()}."
            else:  # HiC
                title = f"Hi-C chromatin interaction in {organism} {tissue}"
                description = f"Chromosome conformation capture (Hi-C) analysis of chromatin interactions in {tissue.lower()} tissue from {organism.lower()}."
            
            # Generate realistic accession numbers
            if database == "ENCODE":
                accession = f"ENCSR{random.randint(100000, 999999)}"
            elif database == "GEO":
                accession = f"GSE{random.randint(10000, 99999)}"
            elif database == "ArrayExpress":
                accession = f"E-MTAB-{random.randint(1000, 9999)}"
            else:  # TCGA
                accession = f"TCGA-{random.choice(['BRCA', 'LUAD', 'COAD', 'STAD'])}-{random.randint(1000, 9999)}"
            
            # Select specific cell line if cell type is "Cell lines"
            if cell_type == "Cell lines":
                specific_cell_line = random.choice(["HeLa", "K562", "HepG2", "MCF-7", "A549", "293T", "Jurkat"])
            else:
                specific_cell_line = None

            dataset = {
                'id': f"dataset_{i}",
                'title': title,
                'description': description,
                'data_type': data_type,
                'tissue': tissue,
                'organism': organism,
                'database': database,
                'accession': accession,
                'cell_type': cell_type,
                'cell_line': specific_cell_line,
                'file_count': random.randint(1, 10),
                'file_size_gb': round(random.uniform(0.1, 50.0), 2),
                'publication_date': f"2023-{random.randint(1, 12):02d}-{random.randint(1, 28):02d}",
                'treatment': random.choice([None, "Control", "LPS", "TNF-α", "Dexamethasone"]) if random.random() > 0.7 else None
            }
            
            datasets.append(dataset)
        
        return datasets
    
    def search_datasets(self, 
                       search_term: str = "", 
                       data_types: Optional[List[str]] = None, 
                       tissues: Optional[List[str]] = None,
                       cell_types: Optional[List[str]] = None,
                       organism: str = "Human") -> List[Dict[str, Any]]:
        """Search for datasets based on criteria"""
        
        results = []
        
        for dataset in self.mock_datasets:
            # Filter by organism
            if dataset['organism'] != organism:
                continue
            
            # Filter by data types
            if data_types and dataset['data_type'] not in data_types:
                continue
            
            # Filter by tissues
            if tissues and dataset['tissue'] not in tissues:
                continue
            
            # Filter by cell types
            if cell_types and cell_types != ["Any"] and dataset.get('cell_type') not in cell_types:
                continue
            
            # Filter by search term - improved search logic
            if search_term and search_term.strip():
                search_lower = search_term.lower().strip()
                searchable_text = (
                    dataset['title'].lower() + " " +
                    dataset['description'].lower() + " " +
                    dataset['tissue'].lower() + " " +
                    dataset['data_type'].lower() + " " +
                    dataset.get('cell_type', '').lower() + " " +
                    (dataset.get('cell_line', '') or '').lower() + " " +
                    dataset['accession'].lower()
                )
                
                # Check if any search term words are found
                search_words = search_lower.split()
                found_words = sum(1 for word in search_words if word in searchable_text)
                
                # Require at least half of the search words to match
                if len(search_words) > 0 and found_words / len(search_words) < 0.5:
                    continue
            
            results.append(dataset.copy())
        
        # Sort by relevance (improved scoring)
        if search_term:
            for result in results:
                score = 0
                search_lower = search_term.lower()
                search_words = search_lower.split()
                
                # Title match gets highest score
                title_matches = sum(1 for word in search_words if word in result['title'].lower())
                score += title_matches * 10
                
                # Description match
                desc_matches = sum(1 for word in search_words if word in result['description'].lower())
                score += desc_matches * 5
                
                # Tissue/data type match
                if any(word in result['tissue'].lower() for word in search_words):
                    score += 8
                if any(word in result['data_type'].lower() for word in search_words):
                    score += 8
                
                # Cell type match
                if any(word in result.get('cell_type', '').lower() for word in search_words):
                    score += 6
                
                # Accession match
                if any(word in result['accession'].lower() for word in search_words):
                    score += 15  # Exact accession matches are very relevant
                
                result['relevance_score'] = score
            
            results.sort(key=lambda x: x.get('relevance_score', 0), reverse=True)
        
        # Limit results
        return results[:50]  # Return top 50 results
    
    def get_dataset_details(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """Get detailed information about a specific dataset"""
        for dataset in self.mock_datasets:
            if dataset['id'] == dataset_id:
                # Add additional details
                detailed_dataset = dataset.copy()
                detailed_dataset.update({
                    'files': self._generate_mock_files(dataset),
                    'metadata': self._generate_mock_metadata(dataset),
                    'download_urls': self._generate_mock_urls(dataset)
                })
                return detailed_dataset
        
        return None
    
    def _generate_mock_files(self, dataset: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate mock file list for a dataset"""
        files = []
        file_count = dataset.get('file_count', 3)
        
        data_type = dataset['data_type']
        
        for i in range(file_count):
            if data_type == "Histone Marks" or data_type == "ChIP-seq":
                file_types = ['bed', 'bedgraph', 'bigwig', 'narrowpeak']
                file_type = random.choice(file_types)
                filename = f"{dataset['accession']}_peaks_{i+1}.{file_type}"
                description = f"Peak calls and signal tracks"
            
            elif data_type == "Gene Expression":
                file_types = ['tsv', 'csv', 'fpkm']
                file_type = random.choice(file_types)
                filename = f"{dataset['accession']}_expression_{i+1}.{file_type}"
                description = f"Gene expression quantification"
            
            elif data_type == "HiC":
                file_types = ['bedpe', 'txt', 'hic']
                file_type = random.choice(file_types)
                filename = f"{dataset['accession']}_interactions_{i+1}.{file_type}"
                description = f"Chromatin interaction data"
            
            else:
                file_type = 'bed'
                filename = f"{dataset['accession']}_data_{i+1}.bed"
                description = "Genomics data file"
            
            files.append({
                'filename': filename,
                'file_type': file_type,
                'size_mb': round(random.uniform(10, 1000), 2),
                'description': description,
                'md5_checksum': f"{''.join(random.choices('0123456789abcdef', k=32))}"
            })
        
        return files
    
    def _generate_mock_metadata(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        """Generate mock metadata for a dataset"""
        metadata = {
            'experiment_type': dataset['data_type'],
            'assembly': random.choice(['GRCh38', 'GRCh37', 'mm10', 'mm9']),
            'read_length': random.choice([50, 75, 100, 150]),
            'sequencing_platform': random.choice(['Illumina HiSeq', 'Illumina NovaSeq', 'Illumina MiSeq']),
            'library_strategy': 'ChIP-Seq' if 'ChIP' in dataset['data_type'] else 'RNA-Seq',
            'library_source': 'GENOMIC' if 'ChIP' in dataset['data_type'] else 'TRANSCRIPTOMIC',
            'total_reads': random.randint(10000000, 100000000),
            'mapped_reads': random.randint(8000000, 90000000),
            'quality_score': round(random.uniform(85, 99), 1)
        }
        
        # Add specific metadata based on data type
        if dataset['data_type'] in ["Histone Marks", "ChIP-seq"]:
            metadata.update({
                'antibody': random.choice(['Anti-H3K4me3', 'Anti-H3K27me3', 'Anti-CTCF', 'Anti-p53']),
                'peak_count': random.randint(1000, 50000),
                'fdr_threshold': 0.05,
                'fold_enrichment': round(random.uniform(2.0, 10.0), 2)
            })
        
        elif dataset['data_type'] == "Gene Expression":
            metadata.update({
                'normalization_method': random.choice(['TPM', 'FPKM', 'DESeq2']),
                'gene_count': random.randint(15000, 25000),
                'expressed_genes': random.randint(10000, 20000),
                'batch_correction': random.choice([True, False])
            })
        
        return metadata
    
    def _generate_mock_urls(self, dataset: Dict[str, Any]) -> Dict[str, str]:
        """Generate mock download URLs"""
        base_urls = {
            'ENCODE': 'https://www.encodeproject.org/files/',
            'GEO': 'https://ftp.ncbi.nlm.nih.gov/geo/series/',
            'ArrayExpress': 'https://www.ebi.ac.uk/arrayexpress/files/',
            'TCGA': 'https://portal.gdc.cancer.gov/files/'
        }
        
        database = dataset['database']
        base_url = base_urls.get(database, 'https://example.com/files/')
        
        return {
            'download_page': f"{base_url}{dataset['accession']}/",
            'ftp_download': f"ftp://ftp.example.com/{dataset['accession']}/",
            'api_endpoint': f"{base_url}api/{dataset['accession']}/metadata"
        }
    
    def get_popular_datasets(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get popular/featured datasets"""
        # Return a curated list of "popular" datasets
        popular = random.sample(self.mock_datasets, min(limit, len(self.mock_datasets)))
        
        # Add popularity metrics
        for dataset in popular:
            dataset['download_count'] = random.randint(100, 10000)
            dataset['citation_count'] = random.randint(10, 500)
            dataset['featured'] = random.random() > 0.7
        
        # Sort by popularity
        popular.sort(key=lambda x: x['download_count'], reverse=True)
        
        return popular
    
    def get_datasets_by_tissue(self, tissue: str) -> List[Dict[str, Any]]:
        """Get all datasets for a specific tissue"""
        return [d for d in self.mock_datasets if d['tissue'].lower() == tissue.lower()]
    
    def get_datasets_by_data_type(self, data_type: str) -> List[Dict[str, Any]]:
        """Get all datasets for a specific data type"""
        return [d for d in self.mock_datasets if d['data_type'] == data_type]
