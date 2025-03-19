import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
import joblib
import os

def random_forest_sp(file_name):

# def load_data(data_path):

    # data = pd.read_csv(os.path.join(data_path, file_name))
    data = file_name

    # Separate features and labels
    X = data.drop(columns=['tree', 'first_eigenvalue', 'sum_eigenvalues'])

    # Normalize input features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    
    # Load pretained model
    
    pretrained_rf_model = joblib.load("C:/pruebas/random_forest_pretrained_model.joblib")
    
    # Optionally, use the pretrained model for prediction on the test set
    pretrained_rf_predictions = pretrained_rf_model.predict(X)
    
    return pretrained_rf_predictions

