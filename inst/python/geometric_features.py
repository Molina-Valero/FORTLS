
import numpy as np
import pandas as pd
from jakteristics import compute_features, FEATURE_NAMES


def get_feature_names():
    """
    Get the names of all available geometric features.
    """
    return FEATURE_NAMES


def compute_geom_features(input_df: pd.DataFrame, 
                          dist: float,
                          features: list,
                          threads: int = -1):
    """
    compute geometric features for a given las file

    input_df : pandas DataFrame that includes the columns 'x', 'y', 'z'
    dist : distance for the computation of geometric features
    features : list of geometric features to be computed
    threads : number of threads to be used for the computation (default is -1, which means use all available cores)

    return : pandas DataFrame with the computed geometric features

    Example :

        res_out = compute_geom_features(input_df = df, dist = 0.1, features = ['verticality', 'linearity', 'planarity'])
    """

    available_features = get_feature_names()
    
    # exception if not all features are available
    for feature in features:
        if feature not in available_features:
            raise ValueError(f"Feature '{feature}' is not available. Available features are {available_features}")
        
    # convert the dataframe to a numpy array
    required_cols = ['x', 'y', 'z']

    # exception if not all required columns are present in the dataframe
    for col in required_cols:
        if col not in input_df.columns:
            raise ValueError(f"Column '{col}' is missing from the input dataframe!")

    points = input_df[required_cols].to_numpy()

    # Compute features (returns a numeric array)
    try:
        computed_features = compute_features(points, 
                                             search_radius=dist, 
                                             feature_names=features, 
                                             num_threads=threads)
    except Exception as e:
        print(f"Error in compute_features: {e}")
        print(f"Points shape: {points.shape}")
        print(f"Features requested: {features}")
        raise

    # Ensure both arrays have compatible dtypes
    points = points.astype(np.float64)
    computed_features = computed_features.astype(np.float64)
    
    # Horizontally stack points and features
    combined_data = np.hstack([points, computed_features])

    # Create DataFrame with proper column names
    df_res = pd.DataFrame(combined_data, columns=['x', 'y', 'z'] + features)

    # Include the 'point' column if it exists in the input_df
    if 'point' in input_df.columns:
        df_res['point'] = input_df['point']

    return df_res