from setuptools import setup, find_packages

setup(
    name="genomics-analysis-platform",
    version="1.0.0",
    description="A comprehensive Streamlit-based web application for genomics data analysis",
    author="Genomics Analysis Platform",
    packages=find_packages(),
    install_requires=[
        "streamlit>=1.28.0",
        "pandas>=2.0.0",
        "numpy>=1.24.0",
        "plotly>=5.15.0",
        "psycopg2-binary>=2.9.0",
        "sqlalchemy>=2.0.0",
    ],
    python_requires=">=3.11",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.11",
    ],
    entry_points={
        "console_scripts": [
            "genomics-platform=app:main",
        ],
    },
)