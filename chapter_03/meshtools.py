#%%
import numpy as np
from scipy.spatial import Delaunay

def rectangle_mesh(p, q, N, M):
    """Generates a Delaunay triangulation for the rectangle
    defined by corner points p and q and (N+1)*(M+1) mesh points
    defined by N and M uniform subintervals in x and y
    direction. Returns (p,t) consisting of the 
    points and connectivity matrix"""

    #TODO: Add kwargs to pass stuff like marker size etc on
    # to plot function
    # TODO: Remove automatic scaling and add fix coordinate system

    x_coord = np.linspace(p[0], q[0], N+1)
    y_coord = np.linspace(p[1], q[1], M+1)

    x, y = np.meshgrid(x_coord, y_coord)
    points = np.array([x.flatten(), y.flatten()]).T

    mesh = Delaunay(points)
    return  (points, mesh.simplices)
  
def unitsquare_mesh(N):
    return rectangle_mesh((0, 0), (1,1), N, N)  

if __name__ == "__main__":

    import plottools as pt
   
    # Create unitsquare mesh
    N = 2
    P, T = unitsquare_mesh(N)

    # Plot mesh
    pt.plot_mesh_2d(P, T)

    # Define rectangle mesh with corner points (0,0) and (1,2)
    p = (0,0)
    q = (1,2)
    N, M = 5, 10
    P, T = rectangle_mesh(p, q, N, M)
    
    print(P)
    print(T)

    # Plot mesh
    pt.plot_mesh_2d(P, T)
    
# %%
