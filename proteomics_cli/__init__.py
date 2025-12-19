"""
Proteomics CLI - Unified PRM Data Processing

A command-line tool for automated proteomics PRM (Parallel Reaction Monitoring) 
data analysis with automatic format detection.
"""

__version__ = "3.0.0"
__author__ = "Neev"

from .process_prm_data import process_prm_data

__all__ = ['process_prm_data']
