#!/usr/bin/env python3
"""
Configuration Management for PRM Analysis

Handles loading and validation of configuration files for replicate grouping
and parsing logic.
"""

from typing import Dict, List, Tuple, Any
from pathlib import Path
import yaml


# Default configuration (matches original hard-coded behavior)
DEFAULT_CONFIG = {
    'replicate_parsing': {
        'delimiter': '_',
        'number_position': -1,
        'dilution_position': 0
    },
    'replicate_groups': [
        [7, 14],
        [17, 24],
        [27, 34],
        [37, 44]
    ],
    'special_cases': {
        'baseline_prefixes': ['col', 'control']
    }
}


class PRMConfig:
    """
    Configuration manager for PRM analysis.
    
    Handles loading from YAML files and provides convenient access to config values.
    """
    
    def __init__(self, config_file: str = None):
        """
        Initialize configuration.
        
        Args:
            config_file: Path to YAML config file. If None, uses default config.
        """
        if config_file and Path(config_file).exists():
            with open(config_file, 'r') as f:
                self.config = yaml.safe_load(f)
            print(f"Loaded configuration from: {config_file}")
        else:
            self.config = DEFAULT_CONFIG.copy()
            if config_file:
                print(f"Warning: Config file not found: {config_file}")
                print("Using default configuration")
    
    @property
    def delimiter(self) -> str:
        """Get replicate name delimiter."""
        return self.config['replicate_parsing']['delimiter']
    
    @property
    def number_position(self) -> int:
        """Get position of replicate number in name."""
        return self.config['replicate_parsing']['number_position']
    
    @property
    def dilution_position(self) -> int:
        """Get position of dilution identifier in name."""
        return self.config['replicate_parsing']['dilution_position']
    
    @property
    def group_ranges(self) -> List[Tuple[int, int]]:
        """Get list of (start, end) tuples for each replicate group."""
        return [tuple(r) for r in self.config['replicate_groups']]
    
    @property
    def baseline_prefixes(self) -> List[str]:
        """Get list of prefixes that indicate baseline/control samples."""
        return self.config['special_cases']['baseline_prefixes']
    
    def get_replicate_number(self, replicate_name: str) -> int:
        """
        Extract replicate number from replicate name using configured position.
        
        Args:
            replicate_name: Full replicate name (e.g., "D1_2ul_24")
            
        Returns:
            Replicate number as integer
        """
        try:
            parts = replicate_name.split(self.delimiter)
            return int(parts[self.number_position])
        except (ValueError, IndexError):
            return 0
    
    def get_dilution(self, replicate_name: str) -> str:
        """
        Extract dilution identifier from replicate name.
        
        Args:
            replicate_name: Full replicate name (e.g., "D1_2ul_24")
            
        Returns:
            Dilution identifier (e.g., "D1")
        """
        parts = replicate_name.split(self.delimiter)
        
        # Check for special baseline prefixes
        first_part = parts[0]
        for prefix in self.baseline_prefixes:
            if first_part.startswith(prefix):
                return 'D0'
        
        # Extract dilution from configured position
        try:
            return parts[self.dilution_position]
        except IndexError:
            return 'D0'
    
    def get_replicate_group(self, replicate_number: int) -> int:
        """
        Assign replicate to a group based on configured ranges.
        
        Args:
            replicate_number: The replicate number
            
        Returns:
            Group number (1-indexed), or 0 if no group matches
        """
        for idx, (start, end) in enumerate(self.group_ranges, 1):
            if start <= replicate_number <= end:
                return idx
        
        # Fallback: try to infer pattern
        if len(self.group_ranges) > 0:
            first_start = self.group_ranges[0][0]
            group_size = self.group_ranges[0][1] - self.group_ranges[0][0] + 1
            gap = self.group_ranges[1][0] - self.group_ranges[0][1] if len(self.group_ranges) > 1 else 10
            
            # Estimate which group this might belong to
            offset = replicate_number - first_start
            if offset >= 0:
                estimated_group = (offset // (group_size + gap)) + 1
                if 1 <= estimated_group <= len(self.group_ranges) + 2:
                    return estimated_group
        
        return 0
    
    def print_summary(self):
        """Print configuration summary."""
        print("\n" + "="*70)
        print("CONFIGURATION SUMMARY")
        print("="*70)
        print(f"Replicate delimiter: '{self.delimiter}'")
        print(f"Number position: {self.number_position}")
        print(f"Dilution position: {self.dilution_position}")
        print(f"\nReplicate groups ({len(self.group_ranges)}):")
        for idx, (start, end) in enumerate(self.group_ranges, 1):
            print(f"  Group {idx}: {start}-{end}")
        print(f"\nBaseline prefixes: {', '.join(self.baseline_prefixes)}")
        print("="*70)


def load_config(config_file: str = None) -> PRMConfig:
    """
    Load configuration from file or use defaults.
    
    Args:
        config_file: Path to YAML config file
        
    Returns:
        PRMConfig instance
    """
    return PRMConfig(config_file)
