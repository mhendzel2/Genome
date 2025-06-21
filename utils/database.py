import os
import pandas as pd
import json
from datetime import datetime
from typing import Dict, List, Any, Optional
from sqlalchemy import create_engine, Column, Integer, String, Text, DateTime, Float, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.postgresql import JSONB

Base = declarative_base()

class Dataset(Base):
    __tablename__ = 'datasets'
    
    id = Column(Integer, primary_key=True)
    title = Column(String(500), nullable=False)
    description = Column(Text)
    data_type = Column(String(100), nullable=False)
    tissue = Column(String(100))
    organism = Column(String(100))
    database = Column(String(100))
    accession = Column(String(100), unique=True)
    cell_type = Column(String(100))
    cell_line = Column(String(100))
    file_count = Column(Integer)
    file_size_gb = Column(Float)
    publication_date = Column(String(20))
    treatment = Column(String(100))
    dataset_metadata = Column(JSONB)
    download_urls = Column(JSONB)
    created_at = Column(DateTime, default=datetime.utcnow)
    
class AnalysisResult(Base):
    __tablename__ = 'analysis_results'
    
    id = Column(Integer, primary_key=True)
    analysis_type = Column(String(100), nullable=False)
    dataset_accessions = Column(JSONB)  # List of dataset accessions used
    parameters = Column(JSONB)
    results = Column(JSONB)
    created_at = Column(DateTime, default=datetime.utcnow)
    user_session = Column(String(100))

class UploadedFile(Base):
    __tablename__ = 'uploaded_files'
    
    id = Column(Integer, primary_key=True)
    filename = Column(String(255), nullable=False)
    data_type = Column(String(100), nullable=False)
    file_size_bytes = Column(Integer)
    validation_status = Column(String(50))
    validation_info = Column(Text)
    file_content_hash = Column(String(64))
    uploaded_at = Column(DateTime, default=datetime.utcnow)
    user_session = Column(String(100))

class DatabaseManager:
    """Manages PostgreSQL database connections and operations"""
    
    def __init__(self):
        self.database_url = os.getenv('DATABASE_URL')
        if not self.database_url:
            raise ValueError("DATABASE_URL environment variable not found")
        
        self.engine = create_engine(self.database_url)
        self.SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=self.engine)
        self.create_tables()
    
    def create_tables(self):
        """Create all tables if they don't exist"""
        Base.metadata.create_all(bind=self.engine)
    
    def get_session(self):
        """Get a database session"""
        return self.SessionLocal()
    
    def add_dataset(self, dataset_data: Dict[str, Any]) -> int:
        """Add a new dataset to the database"""
        session = self.get_session()
        try:
            # Check if dataset already exists
            existing = session.query(Dataset).filter_by(accession=dataset_data['accession']).first()
            if existing:
                return existing.id
            
            dataset = Dataset(
                title=dataset_data['title'],
                description=dataset_data['description'],
                data_type=dataset_data['data_type'],
                tissue=dataset_data['tissue'],
                organism=dataset_data['organism'],
                database=dataset_data['database'],
                accession=dataset_data['accession'],
                cell_type=dataset_data.get('cell_type'),
                cell_line=dataset_data.get('cell_line'),
                file_count=dataset_data.get('file_count'),
                file_size_gb=dataset_data.get('file_size_gb'),
                publication_date=dataset_data.get('publication_date'),
                treatment=dataset_data.get('treatment'),
                dataset_metadata=dataset_data.get('metadata', {}),
                download_urls=dataset_data.get('download_urls', {})
            )
            
            session.add(dataset)
            session.commit()
            return dataset.id
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.close()
    
    def search_datasets(self, 
                       search_term: str = "", 
                       data_types: Optional[List[str]] = None,
                       tissues: Optional[List[str]] = None,
                       cell_types: Optional[List[str]] = None,
                       organism: str = "Human",
                       limit: int = 50) -> List[Dict[str, Any]]:
        """Search datasets in the database"""
        session = self.get_session()
        try:
            query = session.query(Dataset)
            
            # Filter by organism
            if organism:
                query = query.filter(Dataset.organism == organism)
            
            # Filter by data types
            if data_types:
                query = query.filter(Dataset.data_type.in_(data_types))
            
            # Filter by tissues
            if tissues:
                query = query.filter(Dataset.tissue.in_(tissues))
            
            # Filter by cell types
            if cell_types and cell_types != ["Any"]:
                query = query.filter(Dataset.cell_type.in_(cell_types))
            
            # Search term filtering
            if search_term and search_term.strip():
                search_term = search_term.strip()
                query = query.filter(
                    Dataset.title.ilike(f'%{search_term}%') |
                    Dataset.description.ilike(f'%{search_term}%') |
                    Dataset.accession.ilike(f'%{search_term}%') |
                    Dataset.tissue.ilike(f'%{search_term}%') |
                    Dataset.data_type.ilike(f'%{search_term}%')
                )
            
            # Order by creation date (newest first)
            query = query.order_by(Dataset.created_at.desc())
            
            # Limit results
            datasets = query.limit(limit).all()
            
            # Convert to dictionaries
            results = []
            for dataset in datasets:
                result = {
                    'id': f"db_{dataset.id}",
                    'title': dataset.title,
                    'description': dataset.description,
                    'data_type': dataset.data_type,
                    'tissue': dataset.tissue,
                    'organism': dataset.organism,
                    'database': dataset.database,
                    'accession': dataset.accession,
                    'cell_type': dataset.cell_type,
                    'cell_line': dataset.cell_line,
                    'file_count': dataset.file_count,
                    'file_size_gb': dataset.file_size_gb,
                    'publication_date': dataset.publication_date,
                    'treatment': dataset.treatment,
                    'metadata': dataset.dataset_metadata or {},
                    'download_urls': dataset.download_urls or {}
                }
                results.append(result)
            
            return results
        
        finally:
            session.close()
    
    def save_analysis_result(self, analysis_type: str, dataset_accessions: List[str], 
                           parameters: Dict[str, Any], results: Dict[str, Any], 
                           user_session: str) -> int:
        """Save analysis results to the database"""
        session = self.get_session()
        try:
            analysis = AnalysisResult(
                analysis_type=analysis_type,
                dataset_accessions=dataset_accessions,
                parameters=parameters,
                results=results,
                user_session=user_session
            )
            
            session.add(analysis)
            session.commit()
            return analysis.id
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.close()
    
    def get_analysis_history(self, user_session: str) -> List[Dict[str, Any]]:
        """Get analysis history for a user session"""
        session = self.get_session()
        try:
            analyses = session.query(AnalysisResult).filter_by(
                user_session=user_session
            ).order_by(AnalysisResult.created_at.desc()).all()
            
            results = []
            for analysis in analyses:
                result = {
                    'id': analysis.id,
                    'analysis_type': analysis.analysis_type,
                    'dataset_accessions': analysis.dataset_accessions,
                    'parameters': analysis.parameters,
                    'results': analysis.results,
                    'created_at': analysis.created_at
                }
                results.append(result)
            
            return results
        
        finally:
            session.close()
    
    def save_uploaded_file(self, filename: str, data_type: str, file_size: int,
                          validation_status: str, validation_info: str,
                          file_hash: str, user_session: str) -> int:
        """Save uploaded file metadata to the database"""
        session = self.get_session()
        try:
            uploaded_file = UploadedFile(
                filename=filename,
                data_type=data_type,
                file_size_bytes=file_size,
                validation_status=validation_status,
                validation_info=validation_info,
                file_content_hash=file_hash,
                user_session=user_session
            )
            
            session.add(uploaded_file)
            session.commit()
            return uploaded_file.id
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.close()
    
    def populate_sample_datasets(self):
        """Populate the database with sample genomics datasets"""
        from utils.dataset_discovery import DatasetDiscovery
        
        discovery = DatasetDiscovery()
        sample_datasets = discovery._generate_mock_datasets()
        
        session = self.get_session()
        try:
            # Check if datasets already exist
            existing_count = session.query(Dataset).count()
            if existing_count > 0:
                return existing_count
            
            # Add sample datasets
            added_count = 0
            for dataset_data in sample_datasets:
                try:
                    dataset = Dataset(
                        title=dataset_data['title'],
                        description=dataset_data['description'],
                        data_type=dataset_data['data_type'],
                        tissue=dataset_data['tissue'],
                        organism=dataset_data['organism'],
                        database=dataset_data['database'],
                        accession=dataset_data['accession'],
                        cell_type=dataset_data.get('cell_type'),
                        cell_line=dataset_data.get('cell_line'),
                        file_count=dataset_data.get('file_count'),
                        file_size_gb=dataset_data.get('file_size_gb'),
                        publication_date=dataset_data.get('publication_date'),
                        treatment=dataset_data.get('treatment'),
                        dataset_metadata={},
                        download_urls={}
                    )
                    
                    session.add(dataset)
                    added_count += 1
                except Exception as e:
                    print(f"Error adding dataset {dataset_data.get('accession', 'unknown')}: {e}")
                    continue
            
            session.commit()
            return added_count
        
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.close()
    
    def get_database_stats(self) -> Dict[str, Any]:
        """Get database statistics"""
        session = self.get_session()
        try:
            stats = {
                'total_datasets': session.query(Dataset).count(),
                'total_analyses': session.query(AnalysisResult).count(),
                'total_uploaded_files': session.query(UploadedFile).count(),
                'datasets_by_organism': {},
                'datasets_by_data_type': {},
                'datasets_by_database': {}
            }
            
            # Group by organism
            organisms = session.query(Dataset.organism, Dataset.id).all()
            for organism, _ in organisms:
                stats['datasets_by_organism'][organism] = stats['datasets_by_organism'].get(organism, 0) + 1
            
            # Group by data type
            data_types = session.query(Dataset.data_type, Dataset.id).all()
            for data_type, _ in data_types:
                stats['datasets_by_data_type'][data_type] = stats['datasets_by_data_type'].get(data_type, 0) + 1
            
            # Group by database
            databases = session.query(Dataset.database, Dataset.id).all()
            for database, _ in databases:
                stats['datasets_by_database'][database] = stats['datasets_by_database'].get(database, 0) + 1
            
            return stats
        
        finally:
            session.close()