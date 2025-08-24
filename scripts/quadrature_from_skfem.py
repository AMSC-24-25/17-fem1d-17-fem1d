# import skfem
import numpy as np

# Method 1: Use specific quadrature functions
from skfem.quadrature import get_quadrature_line, get_quadrature_tri, get_quadrature_tet

# Try different orders for line
for norder in [1, 2, 4, 6]:
    try:
        points, weights = get_quadrature_line(norder=norder)

        print(f"\nLine quadrature (order {norder}):")
        print(f"Number of points: {points.shape[1]}")
        
        print(f"Detailed points for order {norder}:")
        for i in range(points.shape[1]):
            x = points[:, i]
            w = weights[i]
            # Convert to barycentric coordinates
            bary = [1-x[0], x[0]]
            print(f"Bary=({bary[0]:.10f}, {bary[1]:.10f}), w={w:.10f}")

    except Exception as e:
        print(f"Order {norder}: Error - {e}")

quit(0)
# Try different orders for triangles
for norder in [1, 2, 4]:
    try:
        points, weights = get_quadrature_tri(norder=norder)

        print(f"\nTriangle quadrature (order {norder}):")
        print(f"Number of points: {points.shape[1]}")
        
        print(f"Detailed points for order {norder}:")
        for i in range(points.shape[1]):
            x, y = points[:, i]
            w = weights[i]
            # Convert to barycentric coordinates
            bary = [1-x-y, x, y]
            print(f"Bary=({bary[0]:.10f}, {bary[1]:.10f}, {bary[2]:.10f}), w={w:.10f}")

    except Exception as e:
        print(f"Order {norder}: Error - {e}")

# Try different orders for tetrahedra
for norder in [1, 2, 4]:
    try:
        points, weights = get_quadrature_tet(norder=norder)

        print(f"\nTetrahedron quadrature (order {norder}):")
        print(f"Number of points: {points.shape[1]}")

        print(f"Detailed points for order {norder}:")
        for i in range(points.shape[1]):
            x, y, z = points[:, i]
            w = weights[i]
            # Convert to barycentric coordinates
            bary = [1-x-y-z, x, y, z]
            print(f"Bary=({bary[0]:.10f}, {bary[1]:.10f}, {bary[2]:.10f}, {bary[3]:.10f}), w={w:.10f}")

    except Exception as e:
        print(f"Order {norder}: Error - {e}")