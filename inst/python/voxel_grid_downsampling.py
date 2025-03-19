import numpy as np


def voxel_grid_downsampling(points, voxel_size):
    voxel_grid = {}
    for point_idx, point in enumerate(points):
        voxel_coord = tuple(np.floor(point / voxel_size).astype(int))
        if voxel_coord not in voxel_grid:
            voxel_grid[voxel_coord] = [point]
        else:
            voxel_grid[voxel_coord].append(point)
    
    voxel_centers = []
    for voxel_coord, voxel_points in voxel_grid.items():
        voxel_center = np.mean(voxel_points, axis=0)  # Calculate mean of coordinates
        voxel_centers.append(voxel_center)
    
    return np.array(voxel_centers)
