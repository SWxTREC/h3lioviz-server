import datetime
import glob
import os

import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc

from evolution import Evolution

# Global definitions of variables
# Range for each lookup table
LUT_RANGE = {'Vr': [300, 900],
             'Density': [0, 30],
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
                  'T': [0.2, 0.9],
                  'Br': [0.9, 0.2, 0.9],
                  'Bx': [0.9, 0.2, 0.9],
                  'By': [0.9, 0.2, 0.9],
                  'Bz': [0.9, 0.2, 0.9]}

# Default colormaps to use for the variables
DEFAULT_CMAP = {'Vr': "Plasma (matplotlib)",
                'Density': "Viridis (matplotlib)",
                'T': "Inferno (matplotlib)",
                'Br': "Cool Warm",
                'Bx': "Cool Warm",
                'By': "Cool Warm",
                'Bz': "Cool Warm"}

# Name mapping from frontend to variables in dataset
VARIABLE_MAP = {'velocity': 'Vr',
                'density': 'Density',
                'temperature': 'T',
                'b': 'Br',
                'bx': 'Bx',
                'by': 'By',
                'bz': 'Bz'}

# List of satellite colors
SATELLITE_COLORS = {"earth": [0.0, 0.3333333333333333, 0.0],
                    "stereoa": [177/255, 138/255, 142/255],
                    "stereob": [94/255, 96/255, 185/255]}


class EnlilDataset(pv_protocols.ParaViewWebProtocol):
    def __init__(self, fname):
        """
        Enlil 4D dataset representation in Paraview.

        This class will read in the given NetCDF file that contains
        the Enlil output. It is designed to enable an easy storage
        and access layer to the data.
        """
        # Initialize the PV web protocols
        super().__init__()

        self.evolutions = {x.name: x for x in load_evolution_files(fname)}

        # create a new 'NetCDF Reader'
        self.data = pvs.NetCDFReader(
            registrationName='test_xarray.nc', FileName=[fname])
        self.data.Dimensions = '(longitude, latitude, radius)'

        self.time_string = pvs.Text(registrationName='Time')
        # Don't add in any text right now
        self.time_string.Text = ""
        # Need to keep track of whether the CME/Threshold should be visible
        # so that we can hide 0 value returns
        self._CME_VISIBLE = True
        self._THRESHOLD_VISIBLE = False

        # Create the magnetic field vectors through a PV Function
        self.bvec = pvs.Calculator(registrationName='Bvec', Input=self.data)
        self.bvec.AttributeType = 'Cell Data'
        self.bvec.ResultArrayName = 'Bvec'
        self.bvec.Function = 'Bx*iHat + By*jHat + Bz*kHat'

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
        # self.cme.SamplingBounds = [-1.5, 0, -1.5, 1.5, -1.5, 1.5]

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
            registrationName='Longitude', Input=self.bvec)
        self.lon_slice.SliceType = 'Plane'
        self.lon_slice.HyperTreeGridSlicer = 'Plane'
        self.lon_slice.SliceOffsetValues = [0.0]
        # Init the SliceType values
        self.lon_slice.SliceType.Origin = [0, 0, 0]
        self.lon_slice.SliceType.Normal = [0.0, 0.0, 1.0]
        # TODO: Not sure what this is...
        # init the 'Plane' selected for 'HyperTreeGridSlicer'
        # self.lon_slice.HyperTreeGridSlicer.Origin = [0, 0, 0]

        # create a new 'Stream Tracer' in longitude
        self.lon_streamlines = pvs.StreamTracer(
            registrationName='LonStreamlines', Input=self.lon_slice,
            SeedType='Point Cloud')
        self.lon_streamlines.Vectors = ['CELLS', 'Bvec']
        self.lon_streamlines.SurfaceStreamlines = 1
        self.lon_streamlines.MaximumStreamlineLength = 3.4
        # init the 'Point Cloud' selected for 'SeedType'
        self.lon_streamlines.SeedType.Center = [0, 0, 0]
        self.lon_streamlines.SeedType.NumberOfPoints = 500
        self.lon_streamlines.SeedType.Radius = 0.35

        # create a new 'Glyph' in longitude (Arrow/vectors)
        self.lon_arrows = pvs.Glyph(
            registrationName='Lon-B-Arrows', Input=self.lon_slice,
            GlyphType='Arrow')
        self.lon_arrows.OrientationArray = ['CELLS', 'Bvec']
        self.lon_arrows.ScaleArray = ['POINTS', 'No scale array']
        self.lon_arrows.ScaleFactor = 0.1
        self.lon_arrows.GlyphTransform = 'Transform2'
        self.lon_arrows.MaximumNumberOfSamplePoints = 500

        # Create a Latitude slice
        self.lat_slice = pvs.Slice(
            registrationName='Latitude', Input=self.bvec)
        self.lat_slice.SliceType = 'Plane'
        self.lat_slice.HyperTreeGridSlicer = 'Plane'
        self.lat_slice.SliceOffsetValues = [0.0]
        # init the 'Plane' selected for 'SliceType'
        self.lat_slice.SliceType.Origin = [0, 0, 0]
        self.lat_slice.SliceType.Normal = [0.0, 1.0, 0.0]
        # init the 'Plane' selected for 'HyperTreeGridSlicer'
        # self.lat_slice.HyperTreeGridSlicer.Origin = [0, 0, 0]

        # Dictionary mapping of string names to the object
        self.objs = {s: getattr(self, s) for s in (
            "lon_slice", "lat_slice", "bvec", "cme", "data",
            "lon_arrows", "lon_streamlines", "threshold")}
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
        self.view.CameraParallelScale = 2.1250001580766877
        self.view.BackEnd = 'OSPRay raycaster'
        self.view.OSPRayMaterialLibrary = pvs.GetMaterialLibrary()

        # Time string
        disp = pvs.Show(self.time_string, self.view,
                        'TextSourceRepresentation')

        # TODO: Show the base dataset?
        # pvs.Show(self.data, self.view, 'StructuredGridRepresentation')

        # get color transfer function/color map for 'Bz'
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
        # trace defaults for the display properties.
        disp.Representation = 'Volume'
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'None'
        disp.ScaleFactor = 0.09197479853610144
        disp.SelectScaleArray = 'None'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'None'
        disp.GaussianRadius = 0.004598739926805072
        disp.SetScaleArray = [None, '']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'
        disp.ScalarOpacityFunction = bzPWF
        disp.ScalarOpacityUnitDistance = 0.02090409368521722
        disp.OpacityArrayName = [None, '']

        disp = pvs.Show(self.threshold, self.view,
                        'UniformGridRepresentation')
        self.displays[self.threshold] = disp
        # trace defaults for the display properties.
        disp.Representation = 'Volume'
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'None'
        disp.ScaleFactor = 0.09197479853610144
        disp.SelectScaleArray = 'None'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'None'
        disp.GaussianRadius = 0.004598739926805072
        disp.SetScaleArray = [None, '']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'
        disp.ScalarOpacityFunction = bzPWF
        disp.ScalarOpacityUnitDistance = 0.02090409368521722
        disp.OpacityArrayName = [None, '']
        # TODO: show data from bvec?
        # pvs.Show(self.bvec, self.view, 'StructuredGridRepresentation')

        # Latitude
        disp = pvs.Show(self.lat_slice, self.view,
                        'GeometryRepresentation')
        self.displays[self.lat_slice] = disp

        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['CELLS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'Bvec'
        disp.ScaleFactor = 0.2944486496874232
        disp.SelectScaleArray = 'None'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'None'
        disp.GaussianRadius = 0.01472243248437116
        disp.SetScaleArray = [None, '']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'

        # Longitude
        disp = pvs.Show(self.lon_slice, self.view, 'GeometryRepresentation')
        self.displays[self.lon_slice] = disp
        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['CELLS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'Bvec'
        disp.ScaleFactor = 0.340000014164759
        disp.SelectScaleArray = 'None'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'None'
        disp.GaussianRadius = 0.017000000708237952
        disp.SetScaleArray = [None, '']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = [None, '']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'

        # Longitude streamlines
        disp = pvs.Show(self.lon_streamlines, self.view,
                        'GeometryRepresentation')
        self.displays[self.lon_streamlines] = disp
        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.ColorArrayName = [None, '']
        disp.OSPRayScaleArray = 'AngularVelocity'
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'Normals'
        disp.ScaleFactor = 0.3336487650871277
        disp.SelectScaleArray = 'AngularVelocity'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'AngularVelocity'
        disp.GaussianRadius = 0.016682438254356384
        disp.SetScaleArray = ['POINTS', 'AngularVelocity']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = ['POINTS', 'AngularVelocity']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        disp.ScaleTransferFunction.Points = [
            0.0, 0.0, 0.5, 0.0, 0, 1.0, 0.5, 0.0]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        disp.OpacityTransferFunction.Points = [
            0.0, 0.0, 0.5, 0.0, 0, 1.0, 0.5, 0.0]

        # Longitude B-field vectors
        disp = pvs.Show(self.lon_arrows, self.view,
                        'GeometryRepresentation')
        self.displays[self.lon_arrows] = disp
        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.ScaleFactor = 0.1
        disp.GlyphType = 'Arrow'
        disp.GaussianRadius = 0.005
        disp.AmbientColor = [1, 1, 1]
        disp.ColorArrayName = [None, '']
        disp.DiffuseColor = [1, 1, 1]

        # setup the color legend parameters for each legend in this view
        # get color legend/bar for bzLUT in view self.view
        bzLUTColorBar = pvs.GetScalarBar(bzLUT, self.view)
        bzLUTColorBar.Title = 'Bz'
        bzLUTColorBar.ComponentTitle = ''
        # set color bar visibility
        bzLUTColorBar.Visibility = 1

        # show color legend
        self.displays[self.lon_slice].SetScalarBarVisibility(self.view, True)

        # Set colormaps
        for name in VARIABLE_MAP:
            self.set_colormap(name)

        # hide this data from the default initial view
        for x in [self.lon_slice, self.lon_arrows, self.lon_streamlines,
                  self.threshold]:
            pvs.Hide(x, self.view)

        # restore active source
        pvs.SetActiveSource(None)

    def _setup_satellites(self):
        """
        Initializes the satellites locations and plots them as spheres.
        """
        pv_time = pvs.GetAnimationScene().TimeKeeper.Time
        # The internal time variable on the ViewTime attribute is stored as
        # seconds from 1970-01-01, so we use that epoch directly internally.
        curr_time = (datetime.datetime(1970, 1, 1) +
                     datetime.timedelta(seconds=pv_time))

        for x in SATELLITE_COLORS:
            # Skip this satellite if it isn't in the data
            if x not in self.evolutions:
                continue

            # All satellites are represented as a sphere
            sat = pvs.Sphere()
            setattr(self, x, sat)
            # TODO: What coordinate system do we want x/y/z to be in?
            #       The base model is rotated 180 degrees, should we
            #       automatically rotate it for the users?
            evo = self.evolutions[x]

            sat.Center = evo.get_position(curr_time)
            sat.Radius = 0.025

            disp = pvs.Show(sat, self.view,
                            'GeometryRepresentation')
            # trace defaults for the display properties.
            disp.Representation = 'Surface'
            disp.AmbientColor = SATELLITE_COLORS[x]
            disp.ColorArrayName = [None, '']
            disp.DiffuseColor = SATELLITE_COLORS[x]
            disp.OSPRayScaleArray = 'Normals'
            disp.OSPRayScaleFunction = 'PiecewiseFunction'
            disp.SelectOrientationVectors = 'None'
            disp.ScaleFactor = 0.005000000074505806
            disp.SelectScaleArray = 'None'
            disp.GlyphType = 'Arrow'
            disp.GlyphTableIndexArray = 'None'
            disp.GaussianRadius = 0.0002500000037252903
            disp.SetScaleArray = ['POINTS', 'Normals']
            disp.ScaleTransferFunction = 'PiecewiseFunction'
            disp.OpacityArray = ['POINTS', 'Normals']
            disp.OpacityTransferFunction = 'PiecewiseFunction'
            disp.DataAxesGrid = 'GridAxesRepresentation'
            disp.PolarAxes = 'PolarAxesRepresentation'

        # Sun representation
        self.sun = pvs.Sphere()
        self.sun.Center = [0.0, 0.0, 0.0]
        self.sun.Radius = 0.075
        disp = pvs.Show(self.sun, self.view, 'GeometryRepresentation')

        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.AmbientColor = [0.8313725490196079, 0.8313725490196079, 0.0]
        disp.ColorArrayName = [None, '']
        disp.DiffuseColor = [0.8313725490196079, 0.8313725490196079, 0.0]
        disp.OSPRayScaleArray = 'Normals'
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'None'
        disp.ScaleFactor = 0.020000000298023225
        disp.SelectScaleArray = 'None'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'None'
        disp.GaussianRadius = 0.0010000000149011613
        disp.SetScaleArray = ['POINTS', 'Normals']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = ['POINTS', 'Normals']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'

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
        elif visibility == "off":
            pvs.Hide(self.objs[obj], self.view)
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

        # Update all displays to be colored by this variable
        for obj, disp in self.displays.items():
            if obj == self.lon_arrows:
                # We don't want to update the longitude arrow colors
                continue
            pvs.ColorBy(disp, variable)

        self.update_opacity(variable)
        self.update_lut(variable)

        # hides old scalarbars that aren't in the view and
        # shows the new variable we are using now
        pvs.UpdateScalarBars(self.view)

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

    @exportRpc("pv.enlil.snap_to_view")
    def snap_to_view(self, plane):
        """Snap to the given planar view

        plane : str
            Plane to snap to (ecliptic, meridional)
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
        else:
            raise ValueError('Invalid string, only "ecliptic" or "meridional" '
                             ' are allowed.')
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
    def get_satellite_data(self, sat, variable):
        """
        Returns a time-series of data for the given satellite and variable.

        sat : str
            Name of the satellite (earth, stereoa, stereob)
        variable : str
            Variable of interest

        Returns
        -------
        List of data
        """
        # Transform the variable name from the frontend to the back
        var = VARIABLE_MAP[variable]
        return self.evolutions[sat].get_data(var)

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
        # The internal time variable on the ViewTime attribute is stored as
        # seconds from 1970-01-01, so we use that epoch directly internally.
        curr_time = (datetime.datetime(1970, 1, 1) +
                     datetime.timedelta(seconds=pv_time))
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


def load_evolution_files(fname):
    """
    Loads evolution files relative to the given file.

    fname : str
        File name of the full 4D file. The evolution files
        should all be in the same directory as this.

    Returns
    -------
    A list of Evolution objects.
    """
    # Find all json files in our current directory
    files = glob.glob(os.path.join(os.path.dirname(fname), '*.json'))
    # Iterate over the files and create an Evolution object for each one
    return [Evolution(f) for f in files]
