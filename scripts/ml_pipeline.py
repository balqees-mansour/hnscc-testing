"""
HNSCC Machine Learning Pipeline
================================
Binary classification of HNSCC samples as Primary Tumor vs. Normal tissue
using gene expression data with ElasticNet and Random Forest models.

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
import numpy as np
import argparse
import pickle
import logging
from pathlib import Path

# Scikit-learn imports
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, 
    roc_auc_score, 
    confusion_matrix, 
    classification_report,
    roc_curve,
    auc
)

# Visualization imports
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class HNSCCMLPipeline:
    """
    Machine learning pipeline for HNSCC tumor classification.
    
    This class handles:
    - Data loading and preprocessing
    - Model training and evaluation
    - Feature importance extraction
    - Results visualization and saving
    """
    
    def __init__(self, output_dir="results", random_state=42):
        """
        Initialize the pipeline.
        
        Parameters
        ----------
        output_dir : str
            Directory to save results and models
        random_state : int
            Random seed for reproducibility
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.random_state = random_state
        self.models = {}
        self.results = {}
        self.feature_importance = {}
        
        logger.info(f"Pipeline initialized. Output directory: {self.output_dir}")
    
    def load_data(self, expression_file, metadata_file):
        """
        Load expression matrix and metadata.
        
        Parameters
        ----------
        expression_file : str
            Path to expression matrix (genes x samples)
        metadata_file : str
            Path to metadata CSV with 'sample_type' column
        
        Returns
        -------
        X : pd.DataFrame
            Expression matrix (samples x genes)
        y : pd.Series
            Binary labels (0=Normal, 1=Tumor)
        """
        logger.info("Loading data...")
        
        # Load expression matrix
        X = pd.read_csv(expression_file, index_col=0)
        
        # If genes are in rows, transpose
        if X.shape[0] > X.shape[1]:
            logger.info("Transposing expression matrix (genes in rows)")
            X = X.T
        
        # Load metadata
        metadata = pd.read_csv(metadata_file, index_col=0)
        
        # Align samples
        common_samples = X.index.intersection(metadata.index)
        X = X.loc[common_samples]
        metadata = metadata.loc[common_samples]
        
        logger.info(f"Loaded {X.shape[0]} samples and {X.shape[1]} genes")
        
        # Create binary labels
        # 1 = Primary Tumor, 0 = Normal
        y = (metadata['sample_type'] == 'Primary Tumor').astype(int)
        
        logger.info(f"Class distribution: {y.value_counts().to_dict()}")
        
        return X, y
    
    def generate_placeholder_data(self, n_samples=150, n_genes=500):
        """
        Generate placeholder data for testing without real data.
        
        Parameters
        ----------
        n_samples : int
            Number of samples
        n_genes : int
            Number of genes
        
        Returns
        -------
        X : pd.DataFrame
            Expression matrix
        y : pd.Series
            Binary labels
        """
        logger.info(f"Generating placeholder data: {n_samples} samples, {n_genes} genes")
        
        # Generate random expression data
        np.random.seed(self.random_state)
        X = pd.DataFrame(
            np.random.randn(n_samples, n_genes),
            columns=[f"ENSG{i:06d}" for i in range(n_genes)]
        )
        
        # Create imbalanced labels (more tumors than normal)
        n_normal = int(n_samples * 0.25)
        y = pd.Series([0] * n_normal + [1] * (n_samples - n_normal))
        y = y.sample(frac=1, random_state=self.random_state).reset_index(drop=True)
        
        logger.info(f"Placeholder data class distribution: {y.value_counts().to_dict()}")
        
        return X, y
    
    def preprocess(self, X, y, test_size=0.2):
        """
        Preprocess data: train-test split and feature scaling.
        
        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix
        y : pd.Series
            Binary labels
        test_size : float
            Proportion of data for testing
        
        Returns
        -------
        dict
            Dictionary with train/test data and scaler
        """
        logger.info("Preprocessing data...")
        
        # Stratified train-test split (maintains class distribution)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, 
            test_size=test_size, 
            random_state=self.random_state,
            stratify=y
        )
        
        logger.info(f"Train set: {X_train.shape[0]} samples")
        logger.info(f"Test set: {X_test.shape[0]} samples")
        
        # Feature scaling (important for ElasticNet)
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Convert back to DataFrame for easier handling
        X_train_scaled = pd.DataFrame(X_train_scaled, columns=X.columns, index=X_train.index)
        X_test_scaled = pd.DataFrame(X_test_scaled, columns=X.columns, index=X_test.index)
        
        return {
            'X_train': X_train_scaled,
            'X_test': X_test_scaled,
            'y_train': y_train,
            'y_test': y_test,
            'scaler': scaler,
            'feature_names': X.columns
        }
    
    def train_models(self, data_dict, cv_folds=5):
        """
        Train ElasticNet and Random Forest models.
        
        Parameters
        ----------
        data_dict : dict
            Dictionary with train/test data
        cv_folds : int
            Number of cross-validation folds
        """
        logger.info("Training models...")
        
        X_train = data_dict['X_train']
        X_test = data_dict['X_test']
        y_train = data_dict['y_train']
        y_test = data_dict['y_test']
        
        # Define models
        models_dict = {
            'ElasticNet': LogisticRegression(
                penalty='elasticnet',
                solver='saga',
                l1_ratio=0.5,
                max_iter=1000,
                random_state=self.random_state
            ),
            'RandomForest': RandomForestClassifier(
                n_estimators=100,
                max_depth=15,
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=self.random_state,
                n_jobs=-1
            )
        }
        
        # Train and evaluate each model
        for model_name, model in models_dict.items():
            logger.info(f"Training {model_name}...")
            
            # Cross-validation
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
            cv_scores = cross_val_score(model, X_train, y_train, cv=cv, scoring='roc_auc')
            
            # Train on full training set
            model.fit(X_train, y_train)
            
            # Predictions
            y_pred = model.predict(X_test)
            y_pred_proba = model.predict_proba(X_test)[:, 1]
            
            # Metrics
            accuracy = accuracy_score(y_test, y_pred)
            auc_score = roc_auc_score(y_test, y_pred_proba)
            cm = confusion_matrix(y_test, y_pred)
            
            # Store results
            self.models[model_name] = model
            self.results[model_name] = {
                'cv_auc': cv_scores,
                'cv_auc_mean': cv_scores.mean(),
                'cv_auc_std': cv_scores.std(),
                'test_accuracy': accuracy,
                'test_auc': auc_score,
                'confusion_matrix': cm,
                'y_pred': y_pred,
                'y_pred_proba': y_pred_proba,
                'y_test': y_test
            }
            
            logger.info(f"{model_name} - CV AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
            logger.info(f"{model_name} - Test AUC: {auc_score:.4f}, Accuracy: {accuracy:.4f}")
    
    def extract_feature_importance(self, data_dict):
        """
        Extract feature importance from trained models.
        
        Parameters
        ----------
        data_dict : dict
            Dictionary with feature names
        """
        logger.info("Extracting feature importance...")
        
        feature_names = data_dict['feature_names']
        
        # ElasticNet coefficients
        en_model = self.models['ElasticNet']
        en_importance = pd.DataFrame({
            'gene': feature_names,
            'coefficient': en_model.coef_[0],
            'abs_coefficient': np.abs(en_model.coef_[0])
        }).sort_values('abs_coefficient', ascending=False)
        
        self.feature_importance['ElasticNet'] = en_importance
        
        # Random Forest feature importances
        rf_model = self.models['RandomForest']
        rf_importance = pd.DataFrame({
            'gene': feature_names,
            'importance': rf_model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        self.feature_importance['RandomForest'] = rf_importance
        
        logger.info("Feature importance extracted")
    
    def print_results(self):
        """Print model performance results."""
        logger.info("\n" + "="*60)
        logger.info("MODEL PERFORMANCE SUMMARY")
        logger.info("="*60)
        
        for model_name, results in self.results.items():
            logger.info(f"\n{model_name}:")
            logger.info(f"  Cross-Validation AUC: {results['cv_auc_mean']:.4f} ± {results['cv_auc_std']:.4f}")
            logger.info(f"  Test Set AUC: {results['test_auc']:.4f}")
            logger.info(f"  Test Set Accuracy: {results['test_accuracy']:.4f}")
            logger.info(f"  Confusion Matrix:\n{results['confusion_matrix']}")
        
        logger.info("\n" + "="*60)
        logger.info("TOP 20 FEATURES")
        logger.info("="*60)
        
        for model_name, importance_df in self.feature_importance.items():
            logger.info(f"\n{model_name}:")
            logger.info(importance_df.head(20).to_string())
    
    def save_results(self):
        """Save models, results, and feature importance to disk."""
        logger.info("Saving results...")
        
        # Save models
        models_dir = self.output_dir / 'models'
        models_dir.mkdir(exist_ok=True)
        
        for model_name, model in self.models.items():
            model_path = models_dir / f"{model_name.lower()}_model.pkl"
            with open(model_path, 'wb') as f:
                pickle.dump(model, f)
            logger.info(f"Saved {model_name} model to {model_path}")
        
        # Save feature importance
        data_dir = self.output_dir / 'data'
        data_dir.mkdir(exist_ok=True)
        
        for model_name, importance_df in self.feature_importance.items():
            csv_path = data_dir / f"{model_name.lower()}_feature_importance.csv"
            importance_df.to_csv(csv_path, index=False)
            logger.info(f"Saved {model_name} feature importance to {csv_path}")
        
        # Save performance metrics
        metrics_list = []
        for model_name, results in self.results.items():
            metrics_list.append({
                'Model': model_name,
                'CV_AUC_Mean': results['cv_auc_mean'],
                'CV_AUC_Std': results['cv_auc_std'],
                'Test_AUC': results['test_auc'],
                'Test_Accuracy': results['test_accuracy']
            })
        
        metrics_df = pd.DataFrame(metrics_list)
        metrics_path = data_dir / 'model_performance_metrics.csv'
        metrics_df.to_csv(metrics_path, index=False)
        logger.info(f"Saved performance metrics to {metrics_path}")
    
    def plot_results(self):
        """Generate and save visualization plots."""
        logger.info("Generating plots...")
        
        plots_dir = self.output_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        # Plot 1: Confusion Matrices
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        for idx, (model_name, results) in enumerate(self.results.items()):
            cm = results['confusion_matrix']
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[idx],
                       xticklabels=['Normal', 'Tumor'],
                       yticklabels=['Normal', 'Tumor'])
            axes[idx].set_title(f'{model_name}\nConfusion Matrix')
            axes[idx].set_ylabel('True Label')
            axes[idx].set_xlabel('Predicted Label')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'confusion_matrices.png', dpi=300, bbox_inches='tight')
        logger.info("Saved confusion matrices plot")
        plt.close()
        
        # Plot 2: ROC Curves
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for model_name, results in self.results.items():
            y_test = results['y_test']
            y_pred_proba = results['y_pred_proba']
            
            fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
            roc_auc = auc(fpr, tpr)
            
            ax.plot(fpr, tpr, label=f'{model_name} (AUC = {roc_auc:.3f})', linewidth=2)
        
        ax.plot([0, 1], [0, 1], 'k--', label='Random Classifier', linewidth=1)
        ax.set_xlabel('False Positive Rate', fontsize=12)
        ax.set_ylabel('True Positive Rate', fontsize=12)
        ax.set_title('ROC Curves - HNSCC Tumor Classification', fontsize=14, fontweight='bold')
        ax.legend(loc='lower right', fontsize=11)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'roc_curves.png', dpi=300, bbox_inches='tight')
        logger.info("Saved ROC curves plot")
        plt.close()
        
        # Plot 3: Feature Importance (Top 20)
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        
        for idx, (model_name, importance_df) in enumerate(self.feature_importance.items()):
            top_features = importance_df.head(20)
            
            if model_name == 'ElasticNet':
                colors = ['red' if x < 0 else 'blue' for x in top_features['coefficient']]
                axes[idx].barh(range(len(top_features)), top_features['coefficient'], color=colors)
            else:
                axes[idx].barh(range(len(top_features)), top_features['importance'], color='steelblue')
            
            axes[idx].set_yticks(range(len(top_features)))
            axes[idx].set_yticklabels(top_features['gene'], fontsize=9)
            axes[idx].set_xlabel('Importance', fontsize=11)
            axes[idx].set_title(f'Top 20 Features - {model_name}', fontsize=12, fontweight='bold')
            axes[idx].invert_yaxis()
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'feature_importance.png', dpi=300, bbox_inches='tight')
        logger.info("Saved feature importance plot")
        plt.close()
    
    def run(self, expression_file=None, metadata_file=None, use_placeholder=True, 
            test_size=0.2, cv_folds=5):
        """
        Run the complete ML pipeline.
        
        Parameters
        ----------
        expression_file : str, optional
            Path to expression matrix
        metadata_file : str, optional
            Path to metadata CSV
        use_placeholder : bool
            Use placeholder data if True
        test_size : float
            Proportion for test set
        cv_folds : int
            Number of cross-validation folds
        """
        logger.info("Starting HNSCC ML Pipeline...")
        
        # Load data
        if use_placeholder or (expression_file is None or metadata_file is None):
            logger.info("Using placeholder data")
            X, y = self.generate_placeholder_data()
        else:
            X, y = self.load_data(expression_file, metadata_file)
        
        # Preprocess
        data_dict = self.preprocess(X, y, test_size=test_size)
        
        # Train models
        self.train_models(data_dict, cv_folds=cv_folds)
        
        # Extract features
        self.extract_feature_importance(data_dict)
        
        # Print results
        self.print_results()
        
        # Save results
        self.save_results()
        
        # Plot results
        self.plot_results()
        
        logger.info("Pipeline completed successfully!")


def main():
    """Command-line interface for the pipeline."""
    parser = argparse.ArgumentParser(
        description='HNSCC Machine Learning Pipeline for Tumor Classification'
    )
    
    parser.add_argument('--expression_matrix', type=str, default=None,
                       help='Path to expression matrix CSV (genes x samples)')
    parser.add_argument('--metadata', type=str, default=None,
                       help='Path to metadata CSV with sample_type column')
    parser.add_argument('--output_dir', type=str, default='results',
                       help='Output directory for results')
    parser.add_argument('--test_size', type=float, default=0.2,
                       help='Proportion of data for testing (0-1)')
    parser.add_argument('--cv_folds', type=int, default=5,
                       help='Number of cross-validation folds')
    parser.add_argument('--use_placeholder', action='store_true', default=True,
                       help='Use placeholder data (default: True)')
    parser.add_argument('--random_state', type=int, default=42,
                       help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Initialize and run pipeline
    pipeline = HNSCCMLPipeline(output_dir=args.output_dir, random_state=args.random_state)
    
    pipeline.run(
        expression_file=args.expression_matrix,
        metadata_file=args.metadata,
        use_placeholder=args.use_placeholder,
        test_size=args.test_size,
        cv_folds=args.cv_folds
    )


if __name__ == '__main__':
    main()
