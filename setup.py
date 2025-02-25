# setup.py
from setuptools import setup, find_packages

setup(
    name="mmseqspy",
    version="0.1.7",
    description="A utility library for clustering sequences and splitting data with MMseqs2 and pandas",
    packages=find_packages(),
    install_requires=[
        "pandas",
    ],
    python_requires=">=3.6",
)