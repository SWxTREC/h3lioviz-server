"""
Close the longitude seam of the spherical model grids.

The models cover every longitude, so the first and last planes of grid points
(0 and 360 degrees) sit on top of each other. They are separate points though,
so the cells on either side of the seam are not connected. Converting cell data
to point data then only averages the cells on one side of it, which shows up as
a discontinuity in anything cutting through the seam. The two planes also only
agree to within round-off, so a slice lying exactly along the seam can fall
through the gap between them and lose half of its surface.
"""

import numpy as np
import paraview.simple as pvs
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonCore import VTK_DOUBLE, VTK_FLOAT, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid

# How close (relative to the size of the grid) the first and last planes of
# points need to be for them to be treated as the same longitude
TOLERANCE = 1e-5


def close_seam(data, registrationName=None):
    """
    Close the periodic seam of a structured grid with point data.

    Parameters
    ----------
    data : Paraview source
        Structured grid with point data, e.g. the output of CellDatatoPointData
    registrationName : str
        Name to register the filter under

    Returns
    -------
    Paraview filter with the seam closed. Grids without a seam pass through
    unchanged.
    """
    seam_filter = pvs.ProgrammableFilter(registrationName=registrationName, Input=data)
    seam_filter.Script = (
        "import seam\n"
        "seam.close_grid_seam(self.GetInputDataObject(0, 0), "
        "self.GetOutputDataObject(0))"
    )
    return seam_filter


def close_grid_seam(grid, output):
    """
    Copy *grid* to *output*, joining the two sides of its seam.

    The last plane of seam points is moved onto the first, and the point data
    on both planes is set to the average of the two sides. That average is what
    cell data -> point data would have produced had the cells been connected,
    since each side of a seam point touches the same number of cells.

    Parameters
    ----------
    grid : vtkDataObject
        Grid to close the seam of
    output : vtkDataObject
        Grid to store the result in
    """
    output.ShallowCopy(grid)
    if not isinstance(grid, vtkStructuredGrid) or grid.GetNumberOfPoints() == 0:
        return

    # numpy ordering is (k, j, i) for VTK's (i, j, k) dimensions
    shape = grid.GetDimensions()[::-1]
    points = vtk_to_numpy(grid.GetPoints().GetData()).reshape(*shape, 3)
    axis = _find_seam_axis(points, TOLERANCE * grid.GetLength())
    if axis is None:
        return

    # Copy anything we modify so the upstream filter's output is left alone
    new_points = vtkPoints()
    new_points.DeepCopy(grid.GetPoints())
    planes = _seam_planes(new_points.GetData(), shape, axis)
    planes[-1] = planes[0]
    output.SetPoints(new_points)

    point_data = output.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(i)
        # Only average physical quantities, not masks or string arrays
        if array is None or array.GetDataType() not in (VTK_FLOAT, VTK_DOUBLE):
            continue
        new_array = array.NewInstance()
        new_array.DeepCopy(array)
        planes = _seam_planes(new_array, shape, axis)
        planes[0] = planes[-1] = 0.5 * (planes[0] + planes[-1])
        # Replaces the original array of the same name
        point_data.AddArray(new_array)


def _seam_planes(array, shape, axis):
    """Writable numpy view of a VTK array with the seam axis first"""
    return np.moveaxis(vtk_to_numpy(array).reshape(*shape, -1), axis, 0)


def _find_seam_axis(points, tolerance):
    """The axis whose first and last planes of *points* coincide, or None"""
    for axis in range(3):
        planes = np.moveaxis(points, axis, 0)
        if len(planes) > 2 and np.abs(planes[0] - planes[-1]).max() <= tolerance:
            return axis
    return None
