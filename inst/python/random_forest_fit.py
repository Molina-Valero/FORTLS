import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
import joblib
import os

def random_forest_fit(data_path, file_name):

# def load_data(data_path):

    data = pd.read_csv(os.path.join(data_path, file_name))
    
    # Separate features and labels
    X = data.drop(columns=['species', 'first_eigenvalue', 'sum_eigenvalues'])
    y = data['species']
    
    # Normalize input features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    
    # Convert labels to numerical format
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y)
    
#   return X, y
  
  

# def split_data(X, y):
    
    x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1, random_state=42)
    
#     return x_train, x_val, x_test, y_train, y_val, y_test
  
  

# def create_and_evaluate_model(x_train, x_val, x_test, y_train, y_val, y_test):
    
    model_rf = RandomForestClassifier(n_estimators=100, random_state=42)
    
    # Model training
    model_rf.fit(x_train, y_train)
    
    # Model evaluation
    train_accuracy = accuracy_score(y_train, model_rf.predict(x_train))
    val_accuracy = accuracy_score(y_val, model_rf.predict(x_val))
    test_accuracy = accuracy_score(y_test, model_rf.predict(x_test))
    
    print("Training Accuracy:", train_accuracy)
    print("Validation Accuracy:", val_accuracy)
    print("Test Accuracy:", test_accuracy)
    
#     return model_rf
  
  

# def save_model(model, filename):
    
    joblib.dump(model_rf, os.path.join(data_path, "random_forest_pretrained_model.joblib"))
#   print(f"Model saved to {filename}")

# Main script

# def random_forest_fit(data_path):

# data_path = "C:/pruebas/SpeciesClassification_2.csv"
# X, y = load_data(data_path)
# x_train, x_val, x_test, y_train, y_val, y_test = split_data(X, y)
# model_rf = create_and_evaluate_model(x_train, x_val, x_test, y_train, y_val, y_test)
# save_model(model_rf, "C:/pruebas/random_forest_pretrained_model.joblib")
