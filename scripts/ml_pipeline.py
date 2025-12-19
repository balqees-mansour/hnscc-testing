import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix

# --- 1. Load and Prepare Data (Placeholder) ---
def load_placeholder_data():
    """
    Generates a placeholder dataset with random gene expression data.
    In a real scenario, this function would load the actual expression matrix and metadata.
    """
    n_samples = 150
    n_genes = 500
    X = pd.DataFrame(pd.np.random.rand(n_samples, n_genes), 
                     columns=[f"gene_{i}" for i in range(n_genes)])
    y = pd.Series(pd.np.random.randint(0, 2, n_samples), name="condition")
    return X, y

# --- 2. Preprocessing ---
def preprocess_data(X, y):
    """
    Preprocesses the data, including train-test split and feature scaling.
    """
    # Train-test split (80/20 split, stratified)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Feature scaling
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    return X_train_scaled, X_test_scaled, y_train, y_test, scaler

# --- 3. Model Training and Evaluation ---
def train_and_evaluate(X_train, y_train, X_test, y_test):
    """
    Trains and evaluates ElasticNet and Random Forest models.
    """
    models = {
        "ElasticNet": LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, random_state=42),
        "RandomForest": RandomForestClassifier(random_state=42)
    }

    results = {}

    for name, model in models.items():
        # 5-fold cross-validation
        cv_scores = cross_val_score(model, X_train, y_train, cv=5, scoring='roc_auc')

        # Train the model
        model.fit(X_train, y_train)

        # Make predictions
        y_pred = model.predict(X_test)

        # Evaluate performance
        accuracy = accuracy_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_pred)
        cm = confusion_matrix(y_test, y_pred)

        results[name] = {
            "model": model,
            "cv_auc_mean": cv_scores.mean(),
            "accuracy": accuracy,
            "auc": auc,
            "confusion_matrix": cm
        }

    return results

# --- 4. Feature Importance ---
def get_feature_importance(results):
    """
    Extracts and returns feature importance from the trained models.
    """
    feature_importance = {}

    # ElasticNet coefficients
    en_model = results["ElasticNet"]["model"]
    feature_importance["ElasticNet"] = pd.DataFrame(
        {'feature': X.columns, 'importance': en_model.coef_[0]}
    ).sort_values(by='importance', ascending=False)

    # Random Forest feature importances
    rf_model = results["RandomForest"]["model"]
    feature_importance["RandomForest"] = pd.DataFrame(
        {'feature': X.columns, 'importance': rf_model.feature_importances_}
    ).sort_values(by='importance', ascending=False)

    return feature_importance

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Load data
    X, y = load_placeholder_data()

    # 2. Preprocess data
    X_train_scaled, X_test_scaled, y_train, y_test, scaler = preprocess_data(X, y)

    # 3. Train and evaluate models
    results = train_and_evaluate(X_train_scaled, y_train, X_test_scaled, y_test)

    # 4. Get feature importance
    feature_importance = get_feature_importance(results)

    # --- Print Results ---
    print("--- Model Performance ---")
    for name, res in results.items():
        print(f"\n--- {name} ---")
        print(f"Cross-Validation AUC: {res['cv_auc_mean']:.4f}")
        print(f"Test Set Accuracy: {res['accuracy']:.4f}")
        print(f"Test Set AUC: {res['auc']:.4f}")
        print(f"Confusion Matrix:\n{res['confusion_matrix']}")

    print("\n--- Feature Importance ---")
    for name, importance_df in feature_importance.items():
        print(f"\n--- Top 10 Features for {name} ---")
        print(importance_df.head(10))
