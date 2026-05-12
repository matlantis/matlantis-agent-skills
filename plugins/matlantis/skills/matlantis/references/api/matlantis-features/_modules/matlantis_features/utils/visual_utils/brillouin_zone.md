# Source code for matlantis_features.utils.visual_utils.brillouin_zone
    
    
    import numpy as np
    import plotly.graph_objects as go
    from ase.cell import Cell
    from ase.lattice import identify_lattice
    from scipy.spatial import Voronoi
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.utils.visual_utils.brillouin_zone.get_brillouin_zone_3d.html#matlantis_features.utils.visual_utils.brillouin_zone.get_brillouin_zone_3d)def get_brillouin_zone_3d(
        cell: Cell,
        show_reciprocal_basis: bool = False,
        show_special_k_points: bool = False,
        tolerance: float = 2e-4,
    ) -> go.Figure:
        reciprocal_cell = cell.reciprocal()[:]
    
        px, py, pz = np.tensordot(reciprocal_cell.T, np.mgrid[-1:2, -1:2, -1:2], axes=(0, 0))
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
    
        vor = Voronoi(points)
    
        edges = []
        for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
            if pid[0] == 13 or pid[1] == 13:
                edges.append(vor.vertices[np.r_[rid, [rid[0]]]])
    
        fig = go.Figure()
        for e in edges:
            fig.add_trace(
                go.Scatter3d(
                    x=e[:, 0],
                    y=e[:, 1],
                    z=e[:, 2],
                    line=dict(color="black", width=2),
                    mode="lines",
                    showlegend=False,
                )
            )
        if show_reciprocal_basis:
            for b in reciprocal_cell:
                fig.add_trace(
                    go.Scatter3d(
                        x=[0, b[0] * 0.7],
                        y=[0, b[1] * 0.7],
                        z=[0, b[2] * 0.7],
                        line=dict(color="red", width=2),
                        mode="lines",
                        showlegend=False,
                    )
                )
        if show_special_k_points:
            lattice = identify_lattice(cell, eps=tolerance)[0]
            coordinates = lattice.get_special_points_array().dot(cell.reciprocal()[:])
            labels = [la + " " + cord.__str__() for la, cord in lattice.get_special_points().items()]
            fig.add_trace(
                go.Scatter3d(
                    x=coordinates[:, 0],
                    y=coordinates[:, 1],
                    z=coordinates[:, 2],
                    mode="markers",
                    marker=dict(size=5, color="green"),
                    text=labels,
                    showlegend=False,
                )
            )
    
        fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)
        return fig
    
    
    
