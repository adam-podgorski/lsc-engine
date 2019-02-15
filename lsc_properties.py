import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

def create_initial_vectors(lsc_definition_array):

    vectors = [
    lsc_definition_array[1] - lsc_definition_array[0],
    lsc_definition_array[2] - lsc_definition_array[0],
    lsc_definition_array[3] - lsc_definition_array[0]
    ]

    return vectors

def create_lsc_definition_array(lsc_definition):

    lsc_definition_array = [
        np.array(list(item))
        for item in lsc_definition
    ]

    return lsc_definition_array

def create_points(lsc_definition_array,vectors):

    points = []
    points += lsc_definition_array

    points += [lsc_definition_array[0] + vectors[0] + vectors[1]]
    points += [lsc_definition_array[0] + vectors[0] + vectors[2]]
    points += [lsc_definition_array[0] + vectors[1] + vectors[2]]
    points += [lsc_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    return points

def create_edges(points):
    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]
    return edges

def create_faces(edges):

    faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))

    return faces


def plot_lsc(lsc_definition,fig,ax):

    lsc_definition_array = create_lsc_definition_array(lsc_definition)
    vectors = create_initial_vectors(lsc_definition_array)
    points = create_points(lsc_definition_array,vectors)
    edges = create_edges(points)
    faces = create_faces(edges)
    
    ax.add_collection3d(faces)

    # Plot the points themselves to force the scaling of the axes
    ax.scatter(points[:,0], points[:,1], points[:,2], s=0)
    ax.set_aspect('equal')
