"""
Proteomics CLI - Unified PRM Data Processing

A command-line tool for automated proteomics PRM (Parallel Reaction Monitoring) 
data analysis with automatic format detection.
"""

__version__ = "2.0.0"
__author__ = "Neev"

from .processor import process_prm_unified

__all__ = ['process_prm_unified']
