# setup.py
from setuptools import setup, find_packages

setup(
    name="mmseqspy",
    version="0.1.0",
    description="Python tools for protein sequence embeddings",
    author="Singh Lab",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.5.0",
        "scikit-learn>=1.0.0",
        "h5py>=3.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "pre-commit>=2.20.0",
            "ruff>=0.10.0",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
