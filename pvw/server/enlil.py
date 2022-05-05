import datetime
import glob
import math
import os

import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc

from evolution import Evolution

# Global definitions of variables
# Range for each lookup table
LUT_RANGE = {'Vr': [300, 900],
             'Density': [0, 30],
             'Pressure': [1e5, 1e7],
             'T': [1e4, 1e6],
             'Br': [-10, 10],
             'Bx': [-10, 10],
             'By': [-10, 10],
             'Bz': [-10, 10]}

# Control points for the opacity mapping
# Can be either 2 or 3 values
# 2: Min/max opacity corresponding to the min/max data
# 3: Min, middle, max opacity corresponding to min/center/max data
OPACITY_VALUES = {'Vr': [0.2, 0.9],
                  'Density': [0.2, 0.9],
                  'Pressure': [0.2, 0.9],
                  'T': [0.2, 0.9],
                  'Br': [0.9, 0.2, 0.9],
                  'Bx': [0.9, 0.2, 0.9],
                  'By': [0.9, 0.2, 0.9],
                  'Bz': [0.9, 0.2, 0.9]}

# Default colormaps to use for the variables
DEFAULT_CMAP = {'Vr': "Plasma (matplotlib)",
                'Density': "Viridis (matplotlib)",
                'Pressure': "Viridis (matplotlib)",
                'T': "Inferno (matplotlib)",
                'Br': "Cool to Warm",
                'Bx': "Cool to Warm",
                'By': "Cool to Warm",
                'Bz': "Cool to Warm"}

# Name mapping from frontend to variables in dataset
VARIABLE_MAP = {'velocity': 'Vr',
                'density': 'Density',
                'pressure': 'Pressure',
                'temperature': 'T',
                'b': 'Br',
                'bx': 'Bx',
                'by': 'By',
                'bz': 'Bz'}

VARIABLE_LABEL = {'velocity': 'Velocity (km/s)',
                  'density': 'Density (r$^2$N/cm$^3$)',
                  'pressure': 'Ram pressure (r$^2$N/cm$^3$ * km$^2$/s$^2$)',
                  'temperature': 'Temperature (K)',
                  'b': 'Br (nT)',
                  'bx': 'Bx (nT)',
                  'by': 'By (nT)',
                  'bz': 'Bz (nT)'}

# List of satellite colors
SATELLITE_COLORS = {"earth": [0.0, 0.3333333333333333, 0.0],
                    "stereoa": [177/255, 138/255, 142/255],
                    "stereob": [94/255, 96/255, 185/255]}


class EnlilDataset(pv_protocols.ParaViewWebProtocol):
    def __init__(self, dirname):
        """
        Enlil 4D dataset representation in Paraview.

        This class will read in the given NetCDF file that contains
        the Enlil output. It is designed to enable an easy storage
        and access layer to the data.
        """
        # Initialize the PV web protocols
        super().__init__()
        # Save the data directory
        self._data_dir = dirname
        self.evolutions = {x.name: x for x in load_evolution_files(dirname)}
        # create a new 'NetCDF Reader' from the full data path
        fname = os.path.join(dirname, "pv-data-3d.nc")
        self.data = pvs.NetCDFReader(
            registrationName='enlil-data', FileName=[fname])
        self.data.Dimensions = '(longitude, latitude, radius)'

        self.time_string = pvs.Text(registrationName='Time')
        # Don't add in any text right now
        self.time_string.Text = ""
        # Need to keep track of whether the CME/Threshold should be visible
        # so that we can hide 0 value returns
        self._CME_VISIBLE = True
        self._THRESHOLD_VISIBLE = False

        # create a new 'Threshold' to represent the CME
        self.threshold_cme = pvs.Threshold(registrationName='CME',
                                           Input=self.data)
        # We really only want a minimum value, so just set the maximum high
        self.threshold_cme.ThresholdRange = [1e-5, 1e5]
        # DP is the variable name in Enlil
        self.threshold_cme.Scalars = ['CELLS', 'DP']
        # This resamples the CME to a uniform grid to make Volume rendering
        # work better and faster
        self.cme = pvs.ResampleToImage(registrationName='resampled_cme',
                                       Input=self.threshold_cme)

        # Move to point arrays for contouring
        self._cme_points = pvs.CellDatatoPointData(
            registrationName=f"CME-points",
            Input=self.threshold_cme)
        self._cme_points.ProcessAllArrays = 1

        self.cme_contours = pvs.Contour(
            registrationName='CME-contour',
            Input=self._cme_points)
        self.cme_contours.ContourBy = ['POINTS', 'Density']
        self.cme_contours.Isosurfaces = []
        self.cme_contours.PointMergeMethod = 'Uniform Binning'

        # Create a threshold that can be modified by the user
        self.threshold_data = pvs.Threshold(registrationName='Threshold',
                                            Input=self.data)
        # We really only want a minimum value, so just set the maximum high
        self.threshold_data.ThresholdRange = [10, 500]
        # The quantity of interest
        self.threshold_data.Scalars = ['CELLS', 'Density']
        self.threshold = pvs.ResampleToImage(
            registrationName='resampled_threshold',
            Input=self.threshold_data)

        # Create a Longitude slice
        self.lon_slice = pvs.Slice(
            registrationName='Longitude', Input=self.data)
        self.lon_slice.SliceType = 'Plane'
        self.lon_slice.HyperTreeGridSlicer = 'Plane'
        self.lon_slice.SliceOffsetValues = [0.0]
        self.lon_slice.SliceType.Origin = [0, 0, 0]
        self.lon_slice.SliceType.Normal = [0.0, 0.0, 1.0]
        self._add_streamlines("lon")

        # Create a Latitude slice
        self.lat_slice = pvs.Slice(
            registrationName='Latitude', Input=self.data)
        self.lat_slice.SliceType = 'Plane'
        self.lat_slice.HyperTreeGridSlicer = 'Plane'
        self.lat_slice.SliceOffsetValues = [0.0]
        self.lat_slice.SliceType.Origin = [0, 0, 0]
        self.lat_slice.SliceType.Normal = [0.0, 1.0, 0.0]
        self._add_streamlines("lat")

        # Dictionary mapping of string names to the object
        self.objs = {s: getattr(self, s) for s in (
            "lon_slice", "lat_slice", "cme", "data",
            "lon_arrows", "lon_streamlines", "lat_arrows", "lat_streamlines",
            "threshold", "cme_contours")}
        # Initialize an empty dictionary to store the displays of the objects
        self.displays = {}
        self._setup_views()
        self._setup_satellites()
        self.update(None, None)

        # Call the update function every time the TimeKeeper gets modified
        # The final 1.0 is optional, but sets it as high priority to first
        # do this before other rendering.
        pvs.GetAnimationScene().TimeKeeper.AddObserver(
            "PropertyModifiedEvent", self.update, 1.0)

        # Set the default to have the slice go through Earth
        self.snap_solar_plane("ecliptic")

    def _setup_views(self):
        """Setup the rendering view."""
        # disable automatic camera reset on 'Show'
        pvs._DisableFirstRenderCameraReset()

        # Create a new 'Render View'
        self.view = pvs.CreateView('RenderView')
        self.view.ViewSize = [600, 600]
        self.view.AxesGrid = 'GridAxes3DActor'
        self.view.CenterOfRotation = [0, 0, 0]
        self.view.StereoType = 'Crystal Eyes'
        self.view.CameraPosition = [-3, 3, 3]
        self.view.CameraFocalPoint = [0, 0, 0]
        self.view.CameraViewUp = [0, 0, 1]
        self.view.CameraFocalDisk = 1.0
        self.view.CameraParallelScale = 2
        self.view.BackEnd = 'OSPRay raycaster'
        self.view.OSPRayMaterialLibrary = pvs.GetMaterialLibrary()

        # Time string
        disp = pvs.Show(self.time_string, self.view,
                        'TextSourceRepresentation')

        # get color transfer function/color map for Bz initially
        bzLUT = pvs.GetColorTransferFunction('Bz')
        bzLUT.RGBPoints = [-10, 0.231373, 0.298039, 0.752941,
                           0, 0.865003, 0.865003, 0.865003,
                           10, 0.705882, 0.0156863, 0.14902]
        bzLUT.ScalarRangeInitialized = 1.0
        # get opacity transfer function/opacity map for 'Bz'
        bzPWF = pvs.GetOpacityTransferFunction('Bz')
        bzPWF.Points = [-10, 0.0, 0.5, 0.0,
                        10, 1.0, 0.5, 0.0]
        bzPWF.ScalarRangeInitialized = 1

        # CME Threshold
        disp = pvs.Show(self.cme, self.view,
                        'UniformGridRepresentation')
        self.displays[self.cme] = disp
        disp.Representation = 'Volume'
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.ScalarOpacityFunction = bzPWF
        disp.ScalarOpacityUnitDistance = 0.02
        disp.OpacityArrayName = [None, '']

        # CME Contours
        disp = pvs.Show(self.cme_contours, self.view,
                        'GeometryRepresentation')
        self.displays[self.cme_contours] = disp
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT

        disp = pvs.Show(self.threshold, self.view,
                        'UniformGridRepresentation')
        self.displays[self.threshold] = disp
        disp.Representation = 'Volume'
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.ScalarOpacityFunction = bzPWF
        disp.ScalarOpacityUnitDistance = 0.02
        disp.OpacityArrayName = [None, '']

        # Latitude
        disp = pvs.Show(self.lat_slice, self.view,
                        'GeometryRepresentation')
        self.displays[self.lat_slice] = disp
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['CELLS', 'Bz']
        disp.LookupTable = bzLUT

        # Longitude
        disp = pvs.Show(self.lon_slice, self.view, 'GeometryRepresentation')
        self.displays[self.lon_slice] = disp
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['CELLS', 'Bz']
        disp.LookupTable = bzLUT

        # Streamlines
        disp = pvs.Show(self.lon_streamlines, self.view,
                        'GeometryRepresentation')
        # Add in a magnetic polarity colormap (radial in or out)
        # with two values blue/red
        # separate=True makes sure it doesn't overwrite the Br of the
        # frontend choices
        bpLUT = pvs.GetColorTransferFunction('Br', disp, separate=True)
        bpLUT.RGBPoints = [-1e5, 0.5, 0.5, 0.5,
                           1e5, 0.9, 0.9, 0.9]
        bpLUT.ScalarRangeInitialized = 1.0
        bpLUT.NumberOfTableValues = 2
        self.displays[self.lon_streamlines] = disp
        disp.Representation = 'Surface'
        disp.ColorArrayName = ["POINTS", 'Br']
        disp.LookupTable = bpLUT

        # B-field vectors
        disp = pvs.Show(self.lon_arrows, self.view,
                        'GeometryRepresentation')
        self.displays[self.lon_arrows] = disp
        disp.Representation = 'Surface'
        disp.ColorArrayName = ["POINTS", 'Br']
        disp.LookupTable = bpLUT

        # Set colormaps
        for name in VARIABLE_MAP:
            self.set_colormap(name)

        # hide this data from the default initial view
        for x in [self.lon_slice, self.lon_arrows, self.lon_streamlines,
                  self.threshold, self.cme_contours]:
            pvs.Hide(x, self.view)

        # restore active source
        pvs.SetActiveSource(None)

    def _setup_satellites(self):
        """
        Initializes the satellites locations and plots them as spheres.
        """
        curr_time = self.get_current_time()

        for x in SATELLITE_COLORS:
            # Skip this satellite if it isn't in the data
            if x not in self.evolutions:
                continue

            # All satellites are represented as a box, except the Earth
            if x == "earth":
                sat = pvs.Sphere()
                sat.Radius = 0.025
            else:
                sat = pvs.Box()
                radius = 0.02
                sat.XLength = radius
                sat.YLength = radius
                sat.ZLength = radius
            setattr(self, x, sat)
            # TODO: What coordinate system do we want x/y/z to be in?
            #       The base model is rotated 180 degrees, should we
            #       automatically rotate it for the users?
            evo = self.evolutions[x]
            sat.Center = evo.get_position(curr_time)

            disp = pvs.Show(sat, self.view,
                            'GeometryRepresentation')
            # trace defaults for the display properties.
            disp.Representation = 'Surface'
            disp.AmbientColor = SATELLITE_COLORS[x]
            disp.ColorArrayName = [None, '']
            disp.DiffuseColor = SATELLITE_COLORS[x]

        # Sun representation
        self.sun = pvs.Sphere()
        self.sun.Center = [0.0, 0.0, 0.0]
        self.sun.Radius = 0.074
        self.sun.ThetaResolution = 50
        self.sun.PhiResolution = 50
        disp = pvs.Show(self.sun, self.view, 'GeometryRepresentation')

        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.AmbientColor = [0.8313725490196079, 0.8313725490196079, 0.0]
        disp.ColorArrayName = [None, '']
        disp.DiffuseColor = [0.8313725490196079, 0.8313725490196079, 0.0]

        # Apply an image to the Earth sphere
        self.apply_earth_texture()
        self.apply_solar_texture()

    def _add_streamlines(self, plane):
        """
        Add streamlines to the desired slice plane

        This requires creating several sources and filters.

        1. Ellipse source (Circle at 0.2 AU) in the proper plane
        2. Calculator filter for creating the vector components
        3. CellData -> PointData filter on the plane
        4. StreamTracer with custom source from (1)
        5. Tubes for better display of (4)
        6. Arrows to indicate direction of the arrows

        plane : str
            Name of the plane to add the streamlines on (lon, lat)
        """
        if plane not in ("lon", "lat"):
            raise ValueError("The plane must be one of ('lon', 'lat').")
        # Our stream tracer source needs to have the same plane
        # as our slice, and 0.2 for the radius
        curr_slice = getattr(self, f"{plane}_slice")
        stream_source = pvs.Ellipse(registrationName=f"{plane}-StreamSource")
        stream_source.Center = [0.0, 0.0, 0.0]
        stream_source.Normal = curr_slice.SliceType.Normal
        radius = [0.2, 0, 0] if plane == "lon" else [0, 0, 0.2]
        stream_source.MajorRadiusVector = radius
        # Controls how many streamlines we have
        stream_source.Resolution = 50

        # Create the magnetic field vectors through a PV Function
        bvec = pvs.Calculator(registrationName=f"{plane}-Bvec",
                              Input=curr_slice)
        bvec.AttributeType = 'Cell Data'
        bvec.ResultArrayName = 'Bvec'
        bvec.Function = 'Bx*iHat + By*jHat + Bz*kHat'
        # Make sure we are doing the CellDatatoPointData on the slice
        # directly and not up above, so that it is perfectly on this plane,
        # otherwise the streamlines will be out-of-plane potentially
        point_data = pvs.CellDatatoPointData(
            registrationName=f"{plane}-CellDatatoPointData",
            Input=bvec)
        # Limit the arrays to process to the vector components and
        # direction (Br) for coloring the polarity
        # NOTE: In the processing steps, we multiply
        #       Br * sign(BP), meaning the Br here contains the
        #       polarity implicitly.
        point_data.ProcessAllArrays = 0
        point_data.CellDataArraytoprocess = ['Br', 'Bvec']

        stream_input = pvs.StreamTracerWithCustomSource(
            registrationName=f"{plane}-StreamTracerWithCustomSource",
            Input=point_data,
            SeedSource=stream_source)
        stream_input.Vectors = ['POINTS', 'Bvec']
        stream_input.SurfaceStreamlines = 1
        stream_input.MaximumStreamlineLength = 3.4
        stream_input.ComputeVorticity = 0

        streamlines = pvs.Tube(registrationName=f"{plane}-Streamlines",
                               Input=stream_input)
        streamlines.Capping = 1
        streamlines.Radius = 0.005

        # create a new 'Glyph' in the slice (Arrow/vectors)
        arrows = pvs.Glyph(
            registrationName=f"{plane}-B-Arrows", Input=stream_input,
            GlyphType='Cone')
        arrows.OrientationArray = ['POINTS', 'Bvec']
        arrows.ScaleArray = ['POINTS', 'No scale array']
        arrows.ScaleFactor = 1
        arrows.GlyphType.Resolution = 60
        arrows.GlyphType.Radius = 0.02
        arrows.GlyphType.Height = 0.08
        arrows.GlyphTransform = 'Transform2'
        arrows.GlyphMode = 'Every Nth Point'
        arrows.GlyphMode = 'Uniform Spatial Distribution (Bounds Based)'
        arrows.MaximumNumberOfSamplePoints = 100
        # Store the objects we need for later as attributes on self
        # (self.lon_streamlines / self.lat_streamlines)
        setattr(self, f"{plane}_stream_source", stream_source)
        setattr(self, f"{plane}_streamlines", streamlines)
        setattr(self, f"{plane}_arrows", arrows)

    @exportRpc("pv.enlil.get_available_runs")
    def get_available_runs(self):
        """
        Get a list of available runs to choose from.

        Returns
        -------
        List of available runs
        """
        # We need to go up a directory from where we are currently.
        # This expects a flat list of available runs currently.
        #   /data/run1/pv-data-3d.nc
        #   /data/run2/pv-data-3d.nc
        # If we have loaded /data/run1/pv-data-3d.nc, then this
        # will return ["/data/run1", "/data/run2"] listing all directories
        # up one level from the current data file
        base_dir = os.path.abspath(os.path.join(self._data_dir, ".."))
        dirs = os.listdir(base_dir)
        # Add the base path into the listed directory
        dirs = [os.path.join(base_dir, x) for x in dirs]
        # Now search to see if there is a pv-data-3d.nc in that directory
        # and if not, ignore that entry
        dirs = [x for x in dirs
                if os.path.exists(os.path.join(x, "pv-data-3d.nc"))]
        return dirs

    @exportRpc("pv.enlil.get_variable_range")
    def get_variable_range(self, name):
        """
        Get the range of values for a variable at the current timestep.

        name : str
            Name of variable to colormap all of the surfaces by
        """
        variable = VARIABLE_MAP[name]
        return self.data.CellData.GetArray(variable).GetRange()

    @exportRpc("pv.enlil.directory")
    def update_dataset(self, dirname):
        """
        Change the dataset directory to the one specified by dirname

        dirname : str
            Path to the dataset file (/dirname/pv-data-3d.nc)
        """
        self._data_dir = dirname
        # Update the evolution files associated with the run
        # NOTE: We need to delete the evolutions first, there must be a
        #       dangling reference within the cpp that causes a segfault
        #       if we just update the object dictionary without removal
        del self.evolutions
        self.evolutions = {x.name: x for x in load_evolution_files(dirname)}
        # Update the primary data 3D data file
        self.data.FileName = os.path.join(dirname, "pv-data-3d.nc")
        # Force an update and re-render
        self.data.UpdatePipeline()
        pvs.Render(self.view)

    @exportRpc("pv.enlil.visibility")
    def change_visibility(self, obj, visibility):
        """
        Change the visibility of an object.

        obj : str
            Name of the object to update
        visibility : str ("on", "off")
            What to set the visibility to for the given object
        """
        if visibility == "on":
            pvs.Show(self.objs[obj], self.view)
            if obj == "lon_streamlines":
                # We also want to turn on the vectors
                pvs.Show(self.objs["lon_arrows"], self.view)
        elif visibility == "off":
            pvs.Hide(self.objs[obj], self.view)
            if obj == "lon_streamlines":
                # We also want to turn off the vectors
                pvs.Hide(self.objs["lon_arrows"], self.view)
        else:
            return ["Visibility can only be 'on' or 'off'"]

        # Keep track of whether the CME and Threshold variables are visible
        if obj == "cme":
            self._CME_VISIBLE = {"on": True, "off": False}[visibility]
        if obj == "threshold":
            self._THRESHOLD_VISIBLE = {"on": True, "off": False}[visibility]
        self.update(None, None)

    @exportRpc("pv.enlil.colorby")
    def change_color_variable(self, name):
        """
        Change the visibility of an object.

        name : str
            Name of variable to colormap all of the surfaces by
        """
        # Use a dictionary to map the variable received to the internal name
        variable = VARIABLE_MAP[name]
        label = VARIABLE_LABEL[name]

        # Update all displays to be colored by this variable
        for obj, disp in self.displays.items():
            if obj in (self.lon_arrows, self.lon_streamlines,
                       self.lat_arrows, self.lat_streamlines):
                # We don't want to update the longitude arrow colors
                continue
            pvs.ColorBy(disp, variable)

            # Also update the colorbar orientation information
            if disp.LookupTable is None:
                continue
            cbar = pvs.GetScalarBar(disp.LookupTable, self.view)
            cbar.AutoOrient = 0
            cbar.Orientation = 'Horizontal'
            cbar.TextPosition = 'Ticks left/bottom, annotations right/top'
            cbar.Title = label
            cbar.ComponentTitle = ""

        self.update_opacity(variable)
        self.update_lut(variable)

        # hides old scalarbars that aren't in the view and
        # shows the new variable we are using now
        pvs.UpdateScalarBars(self.view)
        # But we want to hide the streamlines colorbar
        disp = self.displays[self.lon_streamlines]
        disp.SetScalarBarVisibility(self.view, False)

        # restore active source
        pvs.SetActiveSource(None)
        # Render the view
        pvs.Render(self.view)

    def update_lut(self, variable):
        """
        Set the variable range of the lookup table.

        variable : str
            Name of variable to update
        """
        lut = pvs.GetColorTransferFunction(variable)
        lut.RescaleTransferFunction(LUT_RANGE[variable])
        lut.AutomaticRescaleRangeMode = 'Never'

    def update_opacity(self, variable):
        """
        Set the variable range of the opacity lookup table.

        variable : str
            Name of variable to update
        """
        opacity_map = pvs.GetOpacityTransferFunction(variable)
        # Create the control points
        points = OPACITY_VALUES[variable]
        data_range = LUT_RANGE[variable]
        # opacity_map.Points order of flattened list is
        # (data-value, opacity, mid-point, sharpness)
        if len(points) == 2:
            # only min/max values to map
            opacity_map.Points = [data_range[0], points[0], 0.5, 0,
                                  data_range[1], points[1], 0.5, 0]
        elif len(points) == 3:
            # min/mid/max
            opacity_map.Points = [data_range[0], points[0], 0.5, 0,
                                  (data_range[0] + data_range[1]) / 2,
                                  points[1], 0.5, 0,
                                  data_range[1], points[2], 0.5, 0]
        else:
            raise ValueError("Opacity needs 2 or 3 points to map")

    @exportRpc("pv.enlil.set_colormap")
    def set_colormap(self, name, cmap_name=None):
        """
        Set the colormap for the variable.

        name : str
            Name of the variable to set the colormap of
        cmap_name : str
            Name of the colormap to apply
        """
        # Use a dictionary to map the variable received to the internal name
        variable = VARIABLE_MAP[name]
        lut = pvs.GetColorTransferFunction(variable)
        # If cmap_name is None, use the default version
        lut.ApplyPreset(cmap_name or DEFAULT_CMAP[variable])
        lut.EnableOpacityMapping = 1

    @exportRpc("pv.enlil.set_range")
    def set_range(self, name, range):
        """
        Set the range of values used for colormapping.

        name : str
            Name of the variable to set the colormap of
        range : list[2]
            A list of the minimum and maximum values to colormap over
        """
        # Use a dictionary to map the variable received to the internal name
        variable = VARIABLE_MAP[name]
        LUT_RANGE[variable] = range
        self.update_lut(variable)

    @exportRpc("pv.enlil.set_opacity")
    def set_opacity(self, name, range):
        """
        Set the range of values used for opacity-mapping.

        Opacity values are between 0 and 1, and will be applied
        to the current data range of the variables.

        name : str
            Name of the variable to set the opacity of
        range : list[2 or 3]
            A list of the minimum and maximum opacity values if 2 elements
            are given, or the minimum, midpoint, and maximum if 3 elements
            are given.
        """
        # Use a dictionary to map the variable received to the internal name
        variable = VARIABLE_MAP[name]
        OPACITY_VALUES[variable] = range
        self.update_opacity(variable)

    @exportRpc("pv.enlil.set_threshold")
    def set_threshold(self, name, range):
        """
        Set the variable and range of values to be used for the threshold.

        name : str
            Name of the variable to use for the thresholding
        range : list[2]
            A list of the minimum and maximum values to threshold by
        """
        variable = VARIABLE_MAP[name]
        # The quantity of interest
        self.threshold_data.Scalars = ['CELLS', variable]
        self.threshold_data.ThresholdRange = range

    @exportRpc("pv.enlil.set_contours")
    def set_contours(self, name, values):
        """
        Set the variable and a list of values to be used for the contours.

        name : str
            Name of the variable to use for the contouring
        values : list
            A list of the values to contour by
        """
        variable = VARIABLE_MAP[name]
        # The quantity of interest
        self.cme_contours.ContourBy = ['POINTS', variable]
        self.cme_contours.Isosurfaces = values

    @exportRpc("pv.enlil.snap_solar_plane")
    def snap_solar_plane(self, clip):
        """Snap the solar plane to either the solar equator or Sun-Earth plane.

        clip : str
            Name of the plane to snap to. Either "ecliptic" (Earth) or
            "equator" (heliographic equator).
        """
        if clip == "ecliptic":
            if not hasattr(self, "earth"):
                # We don't have earth's location, so just return
                return
            # We have to have an Earth location for this to work
            # (-z, 0, x), y is frozen, so tilt is only in the xz plane
            # -z rotates by 90 degrees to get the normal vector to the plane
            # that contains Earth
            loc = [-self.earth.Center[2], 0, self.earth.Center[0]]
        elif clip == "equator":
            # Standard coordinates, so the plane is purely in xy and
            # z is perpendicular
            loc = [0, 0, 1]
        else:
            raise ValueError('The snapping clip plane must be either '
                             '"ecliptic" or "equator"')

        self.lon_slice.SliceType.Normal = loc
        # Also update the stream source so they stay in-sync
        self.lon_stream_source.Normal = self.lon_slice.SliceType.Normal

    @exportRpc("pv.enlil.rotate_plane")
    def rotate_plane(self, plane, angle):
        """
        Rotate the desired plane to the given angle.

        plane : str
            Plane to rotate (lon, lat).
        angle : float
            Angle (in degrees) for the tilt of the ecliptic/solar plane.
        """
        angle_rad = math.radians(angle)
        x = math.cos(angle_rad)
        y = math.sin(angle_rad)
        if plane == "lon":
            # We want the normal to the plane, so we are really rotating the
            # normal vector here, which swaps the x/z and adds a negative
            # tilt is only in the xz plane
            self.lon_slice.SliceType.Normal = [-y, 0, -x]
            self.lon_stream_source.Normal = self.lon_slice.SliceType.Normal
        elif plane == "lat":
            raise NotImplementedError("Updating the latitudinal plane angle "
                                      "is not implemented.")
            # TODO: Implement the latitudinal adjustment.
            #       Currently this fails with a segfault, which I
            #       assume has to do with some slicing of grid cells changing
            #       sizes when creating a new slice angle, but it is odd that
            #       it only happens for the latitudinal plane...
            # self.lat_slice.SliceType.Normal = [-y, x, 0]
        else:
            raise ValueError("You can only update the 'lon' or 'lat' plane.")

    @exportRpc("pv.enlil.snap_to_view")
    def snap_to_view(self, plane):
        """Snap to the given planar view

        plane : str
            Plane to snap to (ecliptic, meridional, initial)
        """
        if plane == "ecliptic":
            # Go up in Z
            self.view.CameraPosition = [0, 0, 7.5]
            # Point Earth to the right (Y down)
            self.view.CameraViewUp = [0, -1, 0]
        elif plane == "meridional":
            # Go out in Y
            self.view.CameraPosition = [0, 7.5, 0]
            # Point Earth to the right (Z up)
            self.view.CameraViewUp = [0, 0, 1]
        elif plane == "initial":
            # Reset to the initial values
            self.view.CameraPosition = [-3, 3, 3]
            self.view.CameraViewUp = [0, 0, 1]
        else:
            raise ValueError('Invalid string, only "ecliptic", "meridional", '
                             'and "initial" are allowed.')
        # Force the focal point to be the sun
        self.view.CameraFocalPoint = [0, 0, 0]

    @exportRpc("pv.enlil.get_satellite_times")
    def get_satellite_time(self, sat):
        """
        Returns a time-series of data for the given satellite and variable.

        sat : str
            Name of the satellite (earth, stereoa, stereob)
        variable : str
            Variable of interest

        Returns
        -------
        List of times from epoch
        """
        return self.evolutions[sat].get_times()

    @exportRpc("pv.enlil.get_satellite_data")
    def get_satellite_data(self, sat):
        """
        Returns a time-series of data for the given satellite and variable.

        sat : str
            Name of the satellite (earth, stereoa, stereob)

        Returns
        -------
        JSON formatted similarly to a LaTiS-response
        """
        return self.evolutions[sat].as_latis()

    def update(self, caller, event):
        """
        Update function to call every time the time variable has changed.

        The Threshold variables will ruin the view if they are shown and
        there is no data present. So, if they are clicked "on" by the
        frontend, but have no data we still want to force them to be hidden.
        When they have data again, we want to show that value without needing
        to click on/off by the user.
        """
        pv_time = pvs.GetAnimationScene().TimeKeeper.Time
        curr_time = self.get_current_time()
        self.time_string.Text = curr_time.strftime("%Y-%m-%d %H:00")

        for x in SATELLITE_COLORS:
            # Update the satellite positions based on the evolution data
            if hasattr(self, x):
                getattr(self, x).Center = self.evolutions[x].get_position(
                    curr_time)

        # We need to force an update of the filters to populate the data.
        # The CellData[variable] will be None if there is no data calculated
        # based on the thresholding. In that case, we want to hide the object
        # from view.
        pvs.UpdatePipeline(time=pv_time, proxy=self.threshold_cme)
        if (self.cme.Input.CellData['DP'] is None or not self._CME_VISIBLE):
            pvs.Hide(self.cme, self.view)
        else:
            pvs.Show(self.cme, self.view)

        pvs.UpdatePipeline(time=pv_time, proxy=self.threshold_data)
        # Get the variable associated with this threshold operation and see if
        # it is present within CellData
        var = self.threshold_data.Scalars[1]
        if (self.threshold.Input.CellData[var] is None or
                not self._THRESHOLD_VISIBLE):
            pvs.Hide(self.threshold, self.view)
        else:
            pvs.Show(self.threshold, self.view)

        # Update the rotation of the earth image
        self.rotate_earth()
        self.update_solar_image()

    def apply_earth_texture(self):
        """Applies a texture (image) to the Earth sphere.

        This will look for a local image asset to use, and if not found,
        try to go download it for the user.
        """
        if not hasattr(self, "earth"):
            # We don't have earth's location, so just return
            return
        import pathlib
        # Path to the Earth texture on our local system
        # cwd() is where paraview is launched from
        earth_path = pathlib.Path.cwd() / 'pvw' / 'server' / 'assets'
        earth_path /= 'land_shallow_topo_2048.jpg'
        # If we don't have the texture file, go download it.
        if not earth_path.exists():
            # Make the directories if they don't already exist
            earth_path.parent.mkdir(parents=True, exist_ok=True)
            import urllib.request
            url = ("https://eoimages.gsfc.nasa.gov/images/imagerecords/"
                   "57000/57752/land_shallow_topo_2048.jpg")
            # Make a request
            req = urllib.request.urlopen(url)
            # Write out the response to our local file
            with open(earth_path, "wb") as f:
                f.write(req.read())

        # We should have a local image file to use as a texture
        earth_image = pvs.CreateTexture(str(earth_path))
        # Set up the sphere source
        sphere = self.earth
        sphere.ThetaResolution = 50
        # We need to perturb the StartTheta a small amount to not have a
        # seam/mismatch in the texture at 0
        sphere.StartTheta = 1e-3
        sphere.PhiResolution = 50
        # create a new 'Texture Map to Sphere'
        texture_map = pvs.TextureMaptoSphere(registrationName='EarthImage',
                                             Input=sphere)
        texture_map.PreventSeam = 0

        # Move the Earth sphere center back to zero for a translation
        # after rotation later.
        self.earth.Center = [0, 0, 0]

        # We want to rotate the Earth image with the hour of the day
        # To do that, we need to translate our object to the opposite
        # location of what it truly is, then rotate about the origin's
        # z-axis, then translate to the final location after the rotation.
        # TODO: This could be turned into a rotation matrix eventually to
        #       only calculate this at one step, rather than three successive
        #       filters being applied.
        t = pvs.Transform(registrationName='EarthTranslation1',
                          Input=texture_map)
        t.Transform = 'Transform'
        t.Transform.Translate = [0.0, 0.0, 0.0]
        self.earth_translation1 = t
        t = pvs.Transform(registrationName='EarthRotation',
                          Input=self.earth_translation1)
        t.Transform = 'Transform'
        t.Transform.Rotate = [0.0, 0.0, 0.0]
        self.earth_rotation = t
        t = pvs.Transform(registrationName='EarthTranslation2',
                          Input=self.earth_rotation)
        t.Transform = 'Transform'
        t.Transform.Translate = [0.0, 0.0, -1]
        self.earth_translation2 = t

        # show data from the image and hide the plain sphere
        pvs.Hide(self.earth)
        texture_map_disp = pvs.Show(t, self.view,
                                    'GeometryRepresentation')

        # trace defaults for the display properties.
        texture_map_disp.Representation = 'Surface'
        texture_map_disp.ColorArrayName = [None, '']
        texture_map_disp.SelectTCoordArray = 'Texture Coordinates'
        texture_map_disp.SelectNormalArray = 'Normals'
        texture_map_disp.SelectTangentArray = 'None'
        texture_map_disp.Texture = earth_image
        # To get the proper orientation
        texture_map_disp.FlipTextures = 1

    def apply_solar_texture(self):
        """Applies a texture (image) to the Sun.

        This will look for a local image asset to use, and if not found,
        try to go download it for the user.
        """
        solar_dir = self._data_dir + '/solar_images'
        if not os.path.exists(solar_dir):
            return
        # Store a list of the solar images
        self._solar_images = [os.path.basename(x)
                              for x in glob.glob(solar_dir + '/*.jpg')]
        self._solar_images = sorted(self._solar_images)

        # Sun representation
        # Note we don't use the self.sun here because we want this
        # to be slightly larger than the sun sphere when representing it
        sun = pvs.Sphere()
        sun.Center = [0.0, 0.0, 0.0]
        r = 0.075
        sun.Radius = 0.075
        sun.ThetaResolution = 50
        sun.PhiResolution = 50

        # For the solar imagery we want texture map to plane because it is
        # a flat image instead of an unwrapped image.
        texture_map = pvs.TextureMaptoPlane(registrationName='SunImage',
                                            Input=sun)
        # Make the points form a square of radius r
        # Earth is in the -X direction, so we want our image plane to be
        # in the Y-Z direction, with the origin at (+Y, -Z) and the base
        # of the image extending out in the Y direction to (-Y, -Z).
        texture_map.Origin = [0, r, -r]
        texture_map.Point1 = [0, -r, -r]
        texture_map.Point2 = [0, r, r]

        # We also want to clip the sphere so we don't get any wrapping
        # into the back plane
        clip = pvs.Clip(registrationName='ClipSun', Input=texture_map)
        clip.ClipType = 'Plane'
        clip.HyperTreeGridClipper = 'Plane'
        clip.Scalars = ['POINTS', '']
        clip.Invert = 1

        # This is the plane to clip on. It doesn't cover the entire
        # half-sphere, so limit it a little bit in the X direction
        clip.ClipType.Origin = [-0.03, 0.0, 0.0]

        sun_display = pvs.Show(clip, self.view,
                               'GeometryRepresentation')

        # Create a texture from the first image
        sun_texture = pvs.CreateTexture(
            solar_dir + '/' + self._solar_images[0])
        self._previous_time = self.get_current_time()

        # trace defaults for the display properties.
        sun_display.Representation = 'Surface'
        sun_display.ColorArrayName = [None, '']
        sun_display.SelectTCoordArray = 'Texture Coordinates'
        sun_display.SelectNormalArray = 'Normals'
        sun_display.SelectTangentArray = 'None'
        sun_display.Texture = sun_texture
        # This hides the HMI image when looking from behind
        sun_display.BackfaceRepresentation = 'Cull Backface'
        self.sun_display = sun_display

    def get_current_time(self):
        """Retrieves the current time of the view.

        Returns
        -------
        datetime of the current timestep
        """
        pv_time = pvs.GetAnimationScene().TimeKeeper.Time
        # The internal time variable on the ViewTime attribute is stored as
        # seconds from 1970-01-01, so we use that epoch directly internally.
        return (datetime.datetime(1970, 1, 1) +
                datetime.timedelta(seconds=pv_time))

    def update_solar_image(self):
        if not hasattr(self, '_solar_images'):
            # No solar images
            return
        if self._previous_time == self.get_current_time():
            # This is the same timestep, so we don't need to
            # update anything now
            return
        # Iterate through the solar images, choosing
        # the one before this timestep.
        # filename looks like: 20170906_000000_M_color_4k.jpg
        i = len(self._solar_images) - 1
        image_name = self._solar_images[i]
        t = self.get_current_time().strftime("%Y%m%d_%H0000_M_color_4k.jpg")
        while image_name > t and i > 0:
            image_name = self._solar_images[i]
            i -= 1

        # We have our image_name now, so update the texture
        solar_dir = self._data_dir + '/solar_images'
        self.sun_display.Texture = pvs.CreateTexture(
            solar_dir + '/' + image_name)
        # Set the time for the next update
        self._previous_time = self.get_current_time()

    def rotate_earth(self):
        """Rotates the Earth image around with the hour of day."""
        if not hasattr(self, 'earth_rotation'):
            # There is no earth image to rotate
            return
        curr_time = self.get_current_time()
        # We want rotation to be from 0 -> 360
        rot = (curr_time.hour + curr_time.minute/60) / 24 * 360
        # Rotate around Z with the hours of the day
        earth_pos = self.evolutions['earth'].get_position(curr_time)
        # Move it negative first to apply the rotation
        self.earth_translation1.Transform.Translate = [-x for x in earth_pos]
        self.earth_rotation.Transform.Rotate = [0.0, 0.0, rot]
        # Then move it back positive to its actual location
        self.earth_translation2.Transform.Translate = earth_pos


def load_evolution_files(dirname):
    """
    Loads evolution files relative to the given file.

    dirname : str
        Directory path for the evolution files.

    Returns
    -------
    A list of Evolution objects.
    """
    # Find all json files in our current directory
    files = glob.glob(os.path.join(dirname, '*.json'))
    # Iterate over the files and create an Evolution object for each one
    return [Evolution(f) for f in files]
