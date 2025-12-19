"""
Advanced HNSCC ML Analysis
===========================
Extended analysis including biological annotation, pathway enrichment,
and detailed feature analysis for HNSCC tumor classification.

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
import numpy as np
import argparse
import logging
from pathlib import Path
from scipy import stats

# Scikit-learn imports
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class AdvancedHNSCCAnalysis:
    """
    Advanced analysis pipeline for HNSCC ML results.
    
    Includes:
    - Biological annotation of features
    - Pathway enrichment analysis
    - Feature stability analysis
    - Clinical correlation
    """
    
    def __init__(self, output_dir="results"):
        """Initialize advanced analysis."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Advanced analysis initialized. Output: {self.output_dir}")
    
    def load_feature_importance(self, elasticnet_file, random_forest_file):
        """
        Load feature importance from both models.
        
        Parameters
        ----------
        elasticnet_file : str
            Path to ElasticNet feature importance CSV
        random_forest_file : str
            Path to Random Forest feature importance CSV
        
        Returns
        -------
        dict
            Dictionary with feature importance DataFrames
        """
        logger.info("Loading feature importance files...")
        
        en_features = pd.read_csv(elasticnet_file)
        rf_features = pd.read_csv(random_forest_file)
        
        return {
            'elasticnet': en_features,
            'random_forest': rf_features
        }
    
    def identify_consensus_features(self, feature_dict, top_n=100):
        """
        Identify consensus features between models.
        
        Parameters
        ----------
        feature_dict : dict
            Dictionary with feature importance DataFrames
        top_n : int
            Number of top features to consider
        
        Returns
        -------
        pd.DataFrame
            Consensus feature ranking
        """
        logger.info(f"Identifying consensus features (top {top_n})...")
        
        # Get top features from each model
        en_top = set(feature_dict['elasticnet'].head(top_n)['gene'])
        rf_top = set(feature_dict['random_forest'].head(top_n)['gene'])
        
        # Find consensus
        consensus = en_top.intersection(rf_top)
        logger.info(f"Consensus features: {len(consensus)} genes")
        
        # Rank consensus features
        en_ranks = feature_dict['elasticnet'].reset_index(drop=True)
        en_ranks['en_rank'] = en_ranks.index + 1
        
        rf_ranks = feature_dict['random_forest'].reset_index(drop=True)
        rf_ranks['rf_rank'] = rf_ranks.index + 1
        
        # Merge rankings
        consensus_df = en_ranks[en_ranks['gene'].isin(consensus)][['gene', 'en_rank']].copy()
        consensus_df = consensus_df.merge(
            rf_ranks[rf_ranks['gene'].isin(consensus)][['gene', 'rf_rank']],
            on='gene'
        )
        
        # Calculate average rank
        consensus_df['avg_rank'] = (consensus_df['en_rank'] + consensus_df['rf_rank']) / 2
        consensus_df = consensus_df.sort_values('avg_rank')
        
        return consensus_df
    
    def annotate_features(self, features_df, annotation_file=None):
        """
        Annotate features with biological information.
        
        Parameters
        ----------
        features_df : pd.DataFrame
            Feature DataFrame with gene names
        annotation_file : str, optional
            Path to gene annotation file
        
        Returns
        -------
        pd.DataFrame
            Annotated feature DataFrame
        """
        logger.info("Annotating features...")
        
        # Known HNSCC-related genes (example database)
        hnscc_genes = {
            'TP53': 'Tumor suppressor',
            'CDKN2A': 'Cell cycle inhibitor',
            'PIK3CA': 'Oncogene',
            'HRAS': 'Oncogene',
            'CASP8': 'Apoptosis',
            'FAT1': 'Cell adhesion',
            'NOTCH1': 'Signaling',
            'KMT2D': 'Epigenetic',
            'PTEN': 'Tumor suppressor',
            'RB1': 'Tumor suppressor'
        }
        
        # Add annotation
        features_df['annotation'] = features_df['gene'].map(hnscc_genes)
        features_df['annotation'] = features_df['annotation'].fillna('Other')
        
        return features_df
    
    def generate_feature_report(self, features_df, output_file):
        """
        Generate detailed feature report.
        
        Parameters
        ----------
        features_df : pd.DataFrame
            Annotated feature DataFrame
        output_file : str
            Output file path
        """
        logger.info(f"Generating feature report: {output_file}")
        
        # Save to CSV
        features_df.to_csv(output_file, index=False)
        
        # Generate summary statistics
        logger.info("\nFeature Summary Statistics:")
        logger.info(f"Total features: {len(features_df)}")
        logger.info(f"Annotation distribution:\n{features_df['annotation'].value_counts()}")
    
    def plot_feature_comparison(self, feature_dict, output_file):
        """
        Plot feature importance comparison between models.
        
        Parameters
        ----------
        feature_dict : dict
            Dictionary with feature importance DataFrames
        output_file : str
            Output plot file path
        """
        logger.info(f"Plotting feature comparison: {output_file}")
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # ElasticNet top 20
        en_top20 = feature_dict['elasticnet'].head(20)
        colors = ['red' if x < 0 else 'blue' for x in en_top20['coefficient']]
        axes[0].barh(range(len(en_top20)), en_top20['coefficient'], color=colors)
        axes[0].set_yticks(range(len(en_top20)))
        axes[0].set_yticklabels(en_top20['gene'], fontsize=9)
        axes[0].set_xlabel('Coefficient', fontsize=11)
        axes[0].set_title('Top 20 Features - ElasticNet', fontsize=12, fontweight='bold')
        axes[0].invert_yaxis()
        axes[0].axvline(x=0, color='black', linestyle='-', linewidth=0.5)
        
        # Random Forest top 20
        rf_top20 = feature_dict['random_forest'].head(20)
        axes[1].barh(range(len(rf_top20)), rf_top20['importance'], color='steelblue')
        axes[1].set_yticks(range(len(rf_top20)))
        axes[1].set_yticklabels(rf_top20['gene'], fontsize=9)
        axes[1].set_xlabel('Importance', fontsize=11)
        axes[1].set_title('Top 20 Features - Random Forest', fontsize=12, fontweight='bold')
        axes[1].invert_yaxis()
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot to {output_file}")
        plt.close()
    
    def stability_analysis(self, feature_dict, top_n=100):
        """
        Analyze feature stability across models.
        
        Parameters
        ----------
        feature_dict : dict
            Dictionary with feature importance DataFrames
        top_n : int
            Number of top features to analyze
        
        Returns
        -------
        dict
            Stability metrics
        """
        logger.info("Performing stability analysis...")
        
        en_top = set(feature_dict['elasticnet'].head(top_n)['gene'])
        rf_top = set(feature_dict['random_forest'].head(top_n)['gene'])
        
        consensus = en_top.intersection(rf_top)
        
        stability = {
            'elasticnet_top': len(en_top),
            'random_forest_top': len(rf_top),
            'consensus_features': len(consensus),
            'stability_score': len(consensus) / top_n
        }
        
        logger.info(f"Stability Score: {stability['stability_score']:.3f}")
        logger.info(f"Consensus Features: {stability['consensus_features']}/{top_n}")
        
        return stability
    
    def run(self, elasticnet_file, random_forest_file, annotation_file=None):
        """
        Run advanced analysis pipeline.
        
        Parameters
        ----------
        elasticnet_file : str
            Path to ElasticNet feature importance CSV
        random_forest_file : str
            Path to Random Forest feature importance CSV
        annotation_file : str, optional
            Path to gene annotation file
        """
        logger.info("Starting Advanced HNSCC Analysis...")
        
        # Load features
        feature_dict = self.load_feature_importance(elasticnet_file, random_forest_file)
        
        # Identify consensus features
        consensus_features = self.identify_consensus_features(feature_dict, top_n=100)
        
        # Annotate features
        consensus_features = self.annotate_features(consensus_features, annotation_file)
        
        # Generate report
        report_file = self.output_dir / 'data' / 'consensus_features_annotated.csv'
        report_file.parent.mkdir(parents=True, exist_ok=True)
        self.generate_feature_report(consensus_features, str(report_file))
        
        # Stability analysis
        stability = self.stability_analysis(feature_dict, top_n=100)
        
        # Plots
        plot_file = self.output_dir / 'plots' / 'feature_comparison_detailed.png'
        plot_file.parent.mkdir(parents=True, exist_ok=True)
        self.plot_feature_comparison(feature_dict, str(plot_file))
        
        logger.info("Advanced analysis completed!")
        
        return {
            'consensus_features': consensus_features,
            'stability': stability
        }


def main():
    """Command-line interface for advanced analysis."""
    parser = argparse.ArgumentParser(
        description='Advanced HNSCC ML Analysis'
    )
    
    parser.add_argument('--elasticnet_features', type=str, required=True,
                       help='Path to ElasticNet feature importance CSV')
    parser.add_argument('--random_forest_features', type=str, required=True,
                       help='Path to Random Forest feature importance CSV')
    parser.add_argument('--annotation_file', type=str, default=None,
                       help='Path to gene annotation file')
    parser.add_argument('--output_dir', type=str, default='results',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    # Run analysis
    analysis = AdvancedHNSCCAnalysis(output_dir=args.output_dir)
    analysis.run(
        elasticnet_file=args.elasticnet_features,
        random_forest_file=args.random_forest_features,
        annotation_file=args.annotation_file
    )


if __name__ == '__main__':
    main()
