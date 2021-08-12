import subprocess

import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc


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

        # create a new 'NetCDF Reader'
        self.data = pvs.NetCDFReader(
            registrationName='test_xarray.nc', FileName=[fname])
        self.data.Dimensions = '(longitude, latitude, radius)'

        # TODO: Figure out how to get this information using the pvs API
        # We are making two system calls, one for the ncdump command and
        # the second to sift for the line of interest (time:units)
        # XXX
        # Using ncdump inside the container causes issues with incompatible
        # HDF header files. There is some issue with libraries conflicting
        # when including it, so fake the data for now.
        # x = subprocess.run(["ncdump", "-h", fname], capture_output=True)
        # x = subprocess.run(['grep', 'time:units'], input=x.stdout,
        #                    capture_output=True)
        # # split the line and grab the variable which is in quotes
        # # "seconds since 2017-09-07 12:00:10.351562"
        # self.start_time = x.stdout.decode('utf-8').split('"')[1]
        self.start_time = "seconds since 2017-09-07 12:00:10.351562"

        # Create the magnetic field vectors through a PV Function
        self.bvec = pvs.Calculator(registrationName='Bvec', Input=self.data)
        self.bvec.AttributeType = 'Cell Data'
        self.bvec.ResultArrayName = 'Bvec'
        self.bvec.Function = 'Bx*iHat + By*jHat + Bz*kHat'

        # create a new 'Threshold' to represent the CME
        self.threshold = pvs.Threshold(registrationName='CME', Input=self.data)
        # We really only want a minimum value, so just set the maximum high
        self.threshold.ThresholdRange = [1e-5, 1e5]
        # DP is the variable name in Enlil
        self.threshold.Scalars = ['CELLS', 'DP']
        # This resamples the CME to a uniform grid to make Volume rendering
        # work better and faster
        self.cme = pvs.ResampleToImage(registrationName='resampled_data',
                                       Input=self.threshold)
        # self.cme.SamplingBounds = [-1.5, 0, -1.5, 1.5, -1.5, 1.5]

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
        self.lon_arrows.ScaleFactor = 0.35
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
            "lon_arrows", "lon_streamlines")}
        # Initialize an empty dictionary to store the displays of the objects
        self.displays = {}
        self._setup_views()

    def _setup_views(self):
        """Setup the rendering view."""
        # disable automatic camera reset on 'Show'
        pvs._DisableFirstRenderCameraReset()

        # Create a new 'Render View'
        self.view = pvs.CreateView('RenderView')
        # self.view.Background = [38, 55, 90]
        self.view.ViewSize = [1008, 539]
        self.view.AxesGrid = 'GridAxes3DActor'
        self.view.CenterOfRotation = [0.5, 0.0, 0]
        self.view.StereoType = 'Crystal Eyes'
        self.view.CameraPosition = [-5, 5, 5]
        self.view.CameraFocalPoint = [0.5, 0, 0]
        self.view.CameraViewUp = [0.446016916276754,
                                  -0.28531285472026074, 0.8483309998616994]
        self.view.CameraFocalDisk = 1.0
        self.view.CameraParallelScale = 2.1250001580766877
        self.view.BackEnd = 'OSPRay raycaster'
        self.view.OSPRayMaterialLibrary = pvs.GetMaterialLibrary()

        pvs.SetActiveView(None)

        # create new layout object 'Layout #1'
        self.layout = pvs.CreateLayout(name='Layout #1')
        self.layout.AssignView(0, self.view)
        self.layout.SetSize(1008, 539)

        # restore active view
        pvs.SetActiveView(self.view)

        # Earth representation
        self.earth = pvs.Sphere()
        # TODO: What coordinate system do we want x/y/z to be in?
        #       The base model is rotated 180 degrees, should we
        #       automatically rotate it for the users?
        self.earth.Center = [-1.0, 0.0, 0.0]
        self.earth.Radius = 0.025
        disp = pvs.Show(self.earth, self.view,
                        'GeometryRepresentation')
        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.AmbientColor = [0.0, 0.3333333333333333, 0.0]
        disp.ColorArrayName = [None, '']
        disp.DiffuseColor = [0.0, 0.3333333333333333, 0.0]
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

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        disp.ScaleTransferFunction.Points = [-0.9749279022216797,
                                             0.0, 0.5, 0.0,
                                             0.9749279022216797,
                                             1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        disp.OpacityTransferFunction.Points = [-0.9749279022216797,
                                               0.0, 0.5, 0.0,
                                               0.9749279022216797,
                                               1.0, 0.5, 0.0]

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

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        disp.ScaleTransferFunction.Points = [-0.9749279022216797,
                                             0.0, 0.5, 0.0,
                                             0.9749279022216797,
                                             1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        disp.OpacityTransferFunction.Points = [-0.9749279022216797,
                                               0.0, 0.5, 0.0,
                                               0.9749279022216797,
                                               1.0, 0.5, 0.0]

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
        disp.ColorArrayName = ['POINTS', 'Bz']
        disp.LookupTable = bzLUT
        disp.OSPRayScaleArray = 'BP'
        disp.OSPRayScaleFunction = 'PiecewiseFunction'
        disp.SelectOrientationVectors = 'BP'
        disp.ScaleFactor = 0.3804943442344666
        disp.SelectScaleArray = 'BP'
        disp.GlyphType = 'Arrow'
        disp.GlyphTableIndexArray = 'BP'
        disp.GaussianRadius = 0.019024717211723326
        disp.SetScaleArray = ['POINTS', 'BP']
        disp.ScaleTransferFunction = 'PiecewiseFunction'
        disp.OpacityArray = ['POINTS', 'BP']
        disp.OpacityTransferFunction = 'PiecewiseFunction'
        disp.DataAxesGrid = 'GridAxesRepresentation'
        disp.PolarAxes = 'PolarAxesRepresentation'
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        disp.ScaleTransferFunction.Points = [
            -96.6900405883789, 0.0, 0.5, 0.0, 97.67322540283203, 1.0, 0.5, 0.0]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        disp.OpacityTransferFunction.Points = [
            -96.6900405883789, 0.0, 0.5, 0.0, 97.67322540283203, 1.0, 0.5, 0.0]

        # setup the color legend parameters for each legend in this view
        # get color legend/bar for bzLUT in view self.view
        bzLUTColorBar = pvs.GetScalarBar(bzLUT, self.view)
        bzLUTColorBar.Title = 'Bz'
        bzLUTColorBar.ComponentTitle = ''
        # set color bar visibility
        bzLUTColorBar.Visibility = 1

        # show color legend
        self.displays[self.cme].SetScalarBarVisibility(self.view, True)

        # Set colormaps
        for name in VARIABLE_MAP:
            self.set_colormap(name)

        # hide this data from the default initial view
        for x in [self.lon_slice, self.lon_arrows, self.lon_streamlines]:
            pvs.Hide(x, self.view)

        # restore active source
        pvs.SetActiveSource(None)

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
        for disp in self.displays.values():
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

    @exportRpc("pv.enlil.get_start_time")
    def get_start_time(self):
        return [self.start_time]
