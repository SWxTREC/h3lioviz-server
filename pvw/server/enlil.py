import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc


class EnlilDataset(pv_protocols.ParaViewWebProtocol):
    def __init__(self, fname):
        """Enlil 4D dataset representation in Paraview.

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

        # Create the magnetic field vectors through a PV Function
        self.bvec = pvs.Calculator(registrationName='Bvec', Input=self.data)
        self.bvec.AttributeType = 'Cell Data'
        self.bvec.ResultArrayName = 'Bvec'
        self.bvec.Function = 'Bx*iHat + By*jHat + Bz*kHat'

        # create a new 'Threshold' to represent the CME
        self.cme = pvs.Threshold(registrationName='CME', Input=self.data)
        self.cme.Scalars = ['CELLS', 'DP']  # DP is the variable name in Enlil
        # We really only want a minimum value, so just set it high
        self.cme.ThresholdRange = [1e-5, 1e5]

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
                        'UnstructuredGridRepresentation')
        self.displays[self.cme] = disp
        # trace defaults for the display properties.
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['CELLS', 'Bz']
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

        # hide data in view
        pvs.Hide(self.lon_slice, self.view)

        # restore active source
        pvs.SetActiveSource(None)

    @exportRpc("pv.enlil.visibility")
    def change_visibility(self, obj, visibility):
        """Change the visibility of an object.

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
    def change_color_variable(self, variable):
        """Change the visibility of an object.

        variable : str
            Name of variable to colormap all of the surfaces by
        """
        bzLUT = pvs.GetColorTransferFunction('Bz')
        bzLUT.RGBPoints = [-10, 0.231373, 0.298039, 0.752941,
                           0, 0.865003, 0.865003, 0.865003,
                           10, 0.705882, 0.0156863, 0.14902]
        bzLUT.ScalarRangeInitialized = 1.0

        # get color transfer function/color map for 'Density'
        densityLUT = pvs.GetColorTransferFunction('Density')
        densityLUT.AutomaticRescaleRangeMode = 'Never'
        densityLUT.EnableOpacityMapping = 0
        densityLUT.RGBPoints = [0.0, 0.001462, 0.000466, 0.013866, 0.19610000000000027, 0.002267, 0.00127, 0.01857, 0.39215000000000144, 0.003299, 0.002249, 0.024239, 0.5882500000000018, 0.004547, 0.003392, 0.030909, 0.7842999999999991, 0.006006, 0.004692, 0.038558, 0.9803999999999993, 0.007676, 0.006136, 0.046836, 1.1764500000000004, 0.009561, 0.007713, 0.055143, 1.3725500000000008, 0.011663, 0.009417, 0.06346, 1.5686499999999999, 0.013995, 0.011225, 0.071862, 1.7647000000000004, 0.016561, 0.013136, 0.080282, 1.9608000000000008, 0.019373, 0.015133, 0.088767, 2.15685, 0.022447, 0.017199, 0.097327, 2.352949999999999, 0.025793, 0.019331, 0.10593, 2.5490000000000004, 0.029432, 0.021503, 0.114621, 2.7451000000000016, 0.033385, 0.023702, 0.123397, 2.9412000000000007, 0.037668, 0.025921, 0.132232, 3.13725, 0.042253, 0.028139, 0.141141, 3.3333500000000016, 0.046915, 0.030324, 0.150164, 3.5294000000000008, 0.051644, 0.032474, 0.159254, 3.7255, 0.056449, 0.034569, 0.168414, 3.9215499999999994, 0.06134, 0.03659, 0.177642, 4.11765, 0.066331, 0.038504, 0.186962, 4.3137500000000015, 0.071429, 0.040294, 0.196354, 4.5097999999999985, 0.076637, 0.041905, 0.205799, 4.705899999999998, 0.081962, 0.043328, 0.215289, 4.901949999999999, 0.087411, 0.044556, 0.224813, 5.098050000000001, 0.09299, 0.045583, 0.234358, 5.2941, 0.098702, 0.046402, 0.243904, 5.490199999999999, 0.104551, 0.047008, 0.25343, 5.68625, 0.110536, 0.047399, 0.262912, 5.88235, 0.116656, 0.047574, 0.272321, 6.078449999999999, 0.122908, 0.047536, 0.281624, 6.2745, 0.129285, 0.047293, 0.290788, 6.470600000000001, 0.135778, 0.046856, 0.299776, 6.666650000000001, 0.142378, 0.046242, 0.308553, 6.862749999999998, 0.149073, 0.045468, 0.317085, 7.0588000000000015, 0.15585, 0.044559, 0.325338, 7.2549, 0.162689, 0.043554, 0.333277, 7.451, 0.169575, 0.042489, 0.340874, 7.647050000000001, 0.176493, 0.041402, 0.348111, 7.843150000000001, 0.183429, 0.040329, 0.354971, 8.0392, 0.190367, 0.039309, 0.361447, 8.2353, 0.197297, 0.0384, 0.367535, 8.43135, 0.204209, 0.037632, 0.373238, 8.62745, 0.211095, 0.03703, 0.378563, 8.823550000000001, 0.217949, 0.036615, 0.383522, 9.0196, 0.224763, 0.036405, 0.388129, 9.215699999999998, 0.231538, 0.036405, 0.3924, 9.411750000000001, 0.238273, 0.036621, 0.396353, 9.607849999999996, 0.244967, 0.037055, 0.400007, 9.8039, 0.25162, 0.037705, 0.403378, 10.0, 0.258234, 0.038571, 0.406485, 10.196099999999998, 0.26481, 0.039647, 0.409345, 10.39215, 0.271347, 0.040922, 0.411976, 10.58825, 0.27785, 0.042353, 0.414392, 10.7843, 0.284321, 0.043933, 0.416608, 10.9804, 0.290763, 0.045644, 0.418637, 11.17645, 0.297178, 0.04747, 0.420491, 11.372549999999999, 0.303568, 0.049396, 0.422182, 11.568650000000002, 0.309935, 0.051407, 0.423721, 11.7647, 0.316282, 0.05349, 0.425116, 11.9608, 0.32261, 0.055634, 0.426377, 12.15685, 0.328921, 0.057827, 0.427511, 12.352949999999998, 0.335217, 0.06006, 0.428524, 12.549, 0.3415, 0.062325, 0.429425, 12.7451, 0.347771, 0.064616, 0.430217, 12.941200000000002, 0.354032, 0.066925, 0.430906, 13.13725, 0.360284, 0.069247, 0.431497, 13.333349999999996, 0.366529, 0.071579, 0.431994, 13.529399999999999, 0.372768, 0.073915, 0.4324, 13.725499999999998, 0.379001, 0.076253, 0.432719, 13.921550000000002, 0.385228, 0.078591, 0.432955, 14.117650000000001, 0.391453, 0.080927, 0.433109, 14.31375, 0.397674, 0.083257, 0.433183, 14.5098, 0.403894, 0.08558, 0.433179, 14.7059, 0.410113, 0.087896, 0.433098, 14.90195, 0.416331, 0.090203, 0.432943, 15.098050000000002, 0.422549, 0.092501, 0.432714, 15.294100000000002, 0.428768, 0.09479, 0.432412, 15.490200000000003, 0.434987, 0.097069, 0.432039, 15.686249999999996, 0.441207, 0.099338, 0.431594, 15.88235, 0.447428, 0.101597, 0.43108, 16.07845, 0.453651, 0.103848, 0.430498, 16.2745, 0.459875, 0.106089, 0.429846, 16.4706, 0.4661, 0.108322, 0.429125, 16.66665, 0.472328, 0.110547, 0.428334, 16.862750000000005, 0.478558, 0.112764, 0.427475, 17.0588, 0.484789, 0.114974, 0.426548, 17.2549, 0.491022, 0.117179, 0.425552, 17.451, 0.497257, 0.119379, 0.424488, 17.64705, 0.503493, 0.121575, 0.423356, 17.843149999999998, 0.50973, 0.123769, 0.422156, 18.0392, 0.515967, 0.12596, 0.420887, 18.235299999999995, 0.522206, 0.12815, 0.419549, 18.43135, 0.528444, 0.130341, 0.418142, 18.62745, 0.534683, 0.132534, 0.416667, 18.823549999999997, 0.54092, 0.134729, 0.415123, 19.019599999999997, 0.547157, 0.136929, 0.413511, 19.21569999999999, 0.553392, 0.139134, 0.411829, 19.411749999999998, 0.559624, 0.141346, 0.410078, 19.60785, 0.565854, 0.143567, 0.408258, 19.803899999999995, 0.572081, 0.145797, 0.406369, 20.0, 0.578304, 0.148039, 0.404411, 20.196099999999998, 0.584521, 0.150294, 0.402385, 20.392149999999997, 0.590734, 0.152563, 0.40029, 20.588249999999995, 0.59694, 0.154848, 0.398125, 20.7843, 0.603139, 0.157151, 0.395891, 20.980399999999992, 0.60933, 0.159474, 0.393589, 21.176450000000003, 0.615513, 0.161817, 0.391219, 21.372550000000004, 0.621685, 0.164184, 0.388781, 21.568650000000005, 0.627847, 0.166575, 0.386276, 21.764700000000005, 0.633998, 0.168992, 0.383704, 21.9608, 0.640135, 0.171438, 0.381065, 22.156850000000002, 0.64626, 0.173914, 0.378359, 22.352950000000003, 0.652369, 0.176421, 0.375586, 22.549000000000007, 0.658463, 0.178962, 0.372748, 22.745099999999997, 0.66454, 0.181539, 0.369846, 22.9412, 0.670599, 0.184153, 0.366879, 23.13725, 0.676638, 0.186807, 0.363849, 23.333350000000003, 0.682656, 0.189501, 0.360757, 23.5294, 0.688653, 0.192239, 0.357603, 23.72549999999999, 0.694627, 0.195021, 0.354388, 23.92155, 0.700576, 0.197851, 0.351113, 24.117649999999998, 0.7065, 0.200728, 0.347777, 24.313749999999995, 0.712396, 0.203656, 0.344383, 24.509800000000002, 0.718264, 0.206636, 0.340931, 24.705899999999996, 0.724103, 0.20967, 0.337424, 24.90195, 0.729909, 0.212759, 0.333861, 25.09805, 0.735683, 0.215906, 0.330245, 25.294100000000004, 0.741423, 0.219112, 0.326576, 25.4902, 0.747127, 0.222378, 0.322856, 25.68625, 0.752794, 0.225706, 0.319085, 25.88235, 0.758422, 0.229097, 0.315266, 26.078449999999997, 0.76401, 0.232554, 0.311399, 26.2745, 0.769556, 0.236077, 0.307485, 26.4706, 0.775059, 0.239667, 0.303526, 26.666649999999997, 0.780517, 0.243327, 0.299523, 26.862750000000002, 0.785929, 0.247056, 0.295477, 27.058799999999998, 0.791293, 0.250856, 0.29139, 27.2549, 0.796607, 0.254728, 0.287264, 27.450999999999997, 0.801871, 0.258674, 0.283099, 27.64705, 0.807082, 0.262692, 0.278898, 27.84315, 0.812239, 0.266786, 0.274661, 28.039199999999997, 0.817341, 0.270954, 0.27039, 28.235300000000002, 0.822386, 0.275197, 0.266085, 28.43135, 0.827372, 0.279517, 0.26175, 28.62745, 0.832299, 0.283913, 0.257383, 28.823549999999997, 0.837165, 0.288385, 0.252988, 29.0196, 0.841969, 0.292933, 0.248564, 29.2157, 0.846709, 0.297559, 0.244113, 29.411749999999998, 0.851384, 0.30226, 0.239636, 29.607850000000003, 0.855992, 0.307038, 0.235133, 29.8039, 0.860533, 0.311892, 0.230606, 30.0, 0.865006, 0.316822, 0.226055, 30.196100000000005, 0.869409, 0.321827, 0.221482, 30.39215, 0.873741, 0.326906, 0.216886, 30.588250000000006, 0.878001, 0.33206, 0.212268, 30.784299999999998, 0.882188, 0.337287, 0.207628, 30.980400000000007, 0.886302, 0.342586, 0.202968, 31.17645, 0.890341, 0.347957, 0.198286, 31.37255, 0.894305, 0.353399, 0.193584, 31.56864999999999, 0.898192, 0.358911, 0.18886, 31.7647, 0.902003, 0.364492, 0.184116, 31.960799999999995, 0.905735, 0.37014, 0.17935, 32.15685, 0.90939, 0.375856, 0.174563, 32.35295, 0.912966, 0.381636, 0.169755, 32.549, 0.916462, 0.387481, 0.164924, 32.7451, 0.919879, 0.393389, 0.16007, 32.9412, 0.923215, 0.399359, 0.155193, 33.13725, 0.92647, 0.405389, 0.150292, 33.33335, 0.929644, 0.411479, 0.145367, 33.529399999999995, 0.932737, 0.417627, 0.140417, 33.72550000000001, 0.935747, 0.423831, 0.13544, 33.92155, 0.938675, 0.430091, 0.130438, 34.11764999999999, 0.941521, 0.436405, 0.125409, 34.31375, 0.944285, 0.442772, 0.120354, 34.5098, 0.946965, 0.449191, 0.115272, 34.7059, 0.949562, 0.45566, 0.110164, 34.90195, 0.952075, 0.462178, 0.105031, 35.09805, 0.954506, 0.468744, 0.099874, 35.2941, 0.956852, 0.475356, 0.094695, 35.4902, 0.959114, 0.482014, 0.089499, 35.68625, 0.961293, 0.488716, 0.084289, 35.88235000000001, 0.963387, 0.495462, 0.079073, 36.078450000000004, 0.965397, 0.502249, 0.073859, 36.274499999999996, 0.967322, 0.509078, 0.068659, 36.47059999999999, 0.969163, 0.515946, 0.063488, 36.666650000000004, 0.970919, 0.522853, 0.058367, 36.862750000000005, 0.97259, 0.529798, 0.053324, 37.058800000000005, 0.974176, 0.53678, 0.048392, 37.2549, 0.975677, 0.543798, 0.043618, 37.45100000000001, 0.977092, 0.55085, 0.03905, 37.64705, 0.978422, 0.557937, 0.034931, 37.84314999999998, 0.979666, 0.565057, 0.031409, 38.039199999999994, 0.980824, 0.572209, 0.028508, 38.2353, 0.981895, 0.579392, 0.02625, 38.43135000000001, 0.982881, 0.586606, 0.024661, 38.627449999999996, 0.983779, 0.593849, 0.02377, 38.823550000000004, 0.984591, 0.601122, 0.023606, 39.019600000000004, 0.985315, 0.608422, 0.024202, 39.2157, 0.985952, 0.61575, 0.025592, 39.41175, 0.986502, 0.623105, 0.027814, 39.607850000000006, 0.986964, 0.630485, 0.030908, 39.803900000000006, 0.987337, 0.63789, 0.034916, 40.0, 0.987622, 0.64532, 0.039886, 40.19610000000001, 0.987819, 0.652773, 0.045581, 40.39215, 0.987926, 0.66025, 0.05175, 40.58824999999998, 0.987945, 0.667748, 0.058329, 40.784299999999995, 0.987874, 0.675267, 0.065257, 40.9804, 0.987714, 0.682807, 0.072489, 41.176449999999996, 0.987464, 0.690366, 0.07999, 41.37255000000001, 0.987124, 0.697944, 0.087731, 41.56865, 0.986694, 0.70554, 0.095694, 41.76469999999999, 0.986175, 0.713153, 0.103863, 41.960799999999985, 0.985566, 0.720782, 0.112229, 42.156850000000006, 0.984865, 0.728427, 0.120785, 42.35295, 0.984075, 0.736087, 0.129527, 42.549, 0.983196, 0.743758, 0.138453, 42.74510000000001, 0.982228, 0.751442, 0.147565, 42.9412, 0.981173, 0.759135, 0.156863, 43.137249999999995, 0.980032, 0.766837, 0.166353, 43.333349999999996, 0.978806, 0.774545, 0.176037, 43.52940000000001, 0.977497, 0.782258, 0.185923, 43.725500000000004, 0.976108, 0.789974, 0.196018, 43.92155, 0.974638, 0.797692, 0.206332, 44.11765000000001, 0.973088, 0.805409, 0.216877, 44.31375, 0.971468, 0.813122, 0.227658, 44.5098, 0.969783, 0.820825, 0.238686, 44.70590000000001, 0.968041, 0.828515, 0.249972, 44.90195000000001, 0.966243, 0.836191, 0.261534, 45.09805, 0.964394, 0.843848, 0.273391, 45.2941, 0.962517, 0.851476, 0.285546, 45.490199999999994, 0.960626, 0.859069, 0.29801, 45.68625000000001, 0.95872, 0.866624, 0.31082, 45.88235000000002, 0.956834, 0.874129, 0.323974, 46.07845, 0.954997, 0.881569, 0.337475, 46.2745, 0.953215, 0.888942, 0.351369, 46.47060000000001, 0.951546, 0.896226, 0.365627, 46.66664999999999, 0.950018, 0.903409, 0.380271, 46.862749999999984, 0.948683, 0.910473, 0.395289, 47.0588, 0.947594, 0.917399, 0.410665, 47.2549, 0.946809, 0.924168, 0.426373, 47.45099999999998, 0.946392, 0.930761, 0.442367, 47.64705, 0.946403, 0.937159, 0.458592, 47.84315, 0.946903, 0.943348, 0.47497, 48.03920000000001, 0.947937, 0.949318, 0.491426, 48.235299999999995, 0.949545, 0.955063, 0.50786, 48.431349999999995, 0.95174, 0.960587, 0.524203, 48.62745, 0.954529, 0.965896, 0.540361, 48.82355, 0.957896, 0.971003, 0.556275, 49.019600000000004, 0.961812, 0.975924, 0.571925, 49.21570000000001, 0.966249, 0.980678, 0.587206, 49.41175, 0.971162, 0.985282, 0.602154, 49.60784999999999, 0.976511, 0.989753, 0.61676, 49.8039, 0.982257, 0.994109, 0.631017, 50.0, 0.988362, 0.998364, 0.644924]
        densityLUT.NanColor = [0.0, 1.0, 0.0]
        densityLUT.ScalarRangeInitialized = 1.0

        # get color transfer function/color map for 'T'
        tLUT = pvs.GetColorTransferFunction('T')
        tLUT.EnableOpacityMapping = 0
        tLUT.RGBPoints = [1838.11572265625, 0.0, 0.0, 0.0, 798113.6694335938, 0.901960784314, 0.0, 0.0, 1594389.2231445312, 0.901960784314, 0.901960784314, 0.0, 1992527.0, 1.0, 1.0, 1.0]
        tLUT.ColorSpace = 'RGB'
        tLUT.NanColor = [0.0, 0.498039215686, 1.0]
        tLUT.ScalarRangeInitialized = 1.0

        # get color transfer function/color map for 'Vr'
        vrLUT = pvs.GetColorTransferFunction('Vr')
        vrLUT.RGBPoints = [217.9088897705078, 0.001462, 0.000466, 0.013866, 222.26840839193724, 0.002267, 0.00127, 0.01857, 226.62681545838927, 0.003299, 0.002249, 0.024239, 230.98633407981873, 0.004547, 0.003392, 0.030909, 235.34474114627076, 0.006006, 0.004692, 0.038558, 239.7042597677002, 0.007676, 0.006136, 0.046836, 244.06266683415222, 0.009561, 0.007713, 0.055143, 248.42218545558165, 0.011663, 0.009417, 0.06346, 252.7817040770111, 0.013995, 0.011225, 0.071862, 257.1401111434632, 0.016561, 0.013136, 0.080282, 261.4996297648926, 0.019373, 0.015133, 0.088767, 265.8580368313446, 0.022447, 0.017199, 0.097327, 270.21755545277404, 0.025793, 0.019331, 0.10593, 274.57596251922604, 0.029432, 0.021503, 0.114621, 278.9354811406555, 0.033385, 0.023702, 0.123397, 283.294999762085, 0.037668, 0.025921, 0.132232, 287.653406828537, 0.042253, 0.028139, 0.141141, 292.0129254499664, 0.046915, 0.030324, 0.150164, 296.3713325164185, 0.051644, 0.032474, 0.159254, 300.7308511378479, 0.056449, 0.034569, 0.168414, 305.0892582042999, 0.06134, 0.03659, 0.177642, 309.44877682572934, 0.066331, 0.038504, 0.186962, 313.80829544715886, 0.071429, 0.040294, 0.196354, 318.16670251361086, 0.076637, 0.041905, 0.205799, 322.52622113504026, 0.081962, 0.043328, 0.215289, 326.8846282014923, 0.087411, 0.044556, 0.224813, 331.2441468229217, 0.09299, 0.045583, 0.234358, 335.6025538893738, 0.098702, 0.046402, 0.243904, 339.96207251080324, 0.104551, 0.047008, 0.25343, 344.32047957725524, 0.110536, 0.047399, 0.262912, 348.6799981986847, 0.116656, 0.047574, 0.272321, 353.03951682011405, 0.122908, 0.047536, 0.281624, 357.3979238865662, 0.129285, 0.047293, 0.290788, 361.7574425079957, 0.135778, 0.046856, 0.299776, 366.1158495744477, 0.142378, 0.046242, 0.308553, 370.475368195877, 0.149073, 0.045468, 0.317085, 374.83377526232914, 0.15585, 0.044559, 0.325338, 379.1932938837585, 0.162689, 0.043554, 0.333277, 383.55281250518794, 0.169575, 0.042489, 0.340874, 387.91121957164, 0.176493, 0.041402, 0.348111, 392.27073819306946, 0.183429, 0.040329, 0.354971, 396.6291452595215, 0.190367, 0.039309, 0.361447, 400.98866388095087, 0.197297, 0.0384, 0.367535, 405.347070947403, 0.204209, 0.037632, 0.373238, 409.70658956883244, 0.211095, 0.03703, 0.378563, 414.06610819026184, 0.217949, 0.036615, 0.383522, 418.4245152567139, 0.224763, 0.036405, 0.388129, 422.78403387814325, 0.231538, 0.036405, 0.3924, 427.14244094459536, 0.238273, 0.036621, 0.396353, 431.5019595660247, 0.244967, 0.037055, 0.400007, 435.8603666324768, 0.25162, 0.037705, 0.403378, 440.2198852539062, 0.258234, 0.038571, 0.406485, 444.5794038753357, 0.26481, 0.039647, 0.409345, 448.9378109417877, 0.271347, 0.040922, 0.411976, 453.29732956321715, 0.27785, 0.042353, 0.414392, 457.65573662966915, 0.284321, 0.043933, 0.416608, 462.01525525109867, 0.290763, 0.045644, 0.418637, 466.3736623175507, 0.297178, 0.04747, 0.420491, 470.73318093898007, 0.303568, 0.049396, 0.422182, 475.0926995604096, 0.309935, 0.051407, 0.423721, 479.4511066268616, 0.316282, 0.05349, 0.425116, 483.810625248291, 0.32261, 0.055634, 0.426377, 488.169032314743, 0.328921, 0.057827, 0.427511, 492.52855093617245, 0.335217, 0.06006, 0.428524, 496.8869580026245, 0.3415, 0.062325, 0.429425, 501.2464766240539, 0.347771, 0.064616, 0.430217, 505.6059952454835, 0.354032, 0.066925, 0.430906, 509.9644023119354, 0.360284, 0.069247, 0.431497, 514.3239209333649, 0.366529, 0.071579, 0.431994, 518.6823279998168, 0.372768, 0.073915, 0.4324, 523.0418466212462, 0.379001, 0.076253, 0.432719, 527.4002536876983, 0.385228, 0.078591, 0.432955, 531.7597723091277, 0.391453, 0.080927, 0.433109, 536.1192909305573, 0.397674, 0.083257, 0.433183, 540.4776979970092, 0.403894, 0.08558, 0.433179, 544.8372166184388, 0.410113, 0.087896, 0.433098, 549.1956236848907, 0.416331, 0.090203, 0.432943, 553.5551423063202, 0.422549, 0.092501, 0.432714, 557.9135493727722, 0.428768, 0.09479, 0.432412, 562.2730679942016, 0.434987, 0.097069, 0.432039, 566.6314750606537, 0.441207, 0.099338, 0.431594, 570.9909936820831, 0.447428, 0.101597, 0.43108, 575.3505123035126, 0.453651, 0.103848, 0.430498, 579.7089193699646, 0.459875, 0.106089, 0.429846, 584.068437991394, 0.4661, 0.108322, 0.429125, 588.4268450578461, 0.472328, 0.110547, 0.428334, 592.7863636792756, 0.478558, 0.112764, 0.427475, 597.1447707457276, 0.484789, 0.114974, 0.426548, 601.5042893671571, 0.491022, 0.117179, 0.425552, 605.8638079885865, 0.497257, 0.119379, 0.424488, 610.2222150550385, 0.503493, 0.121575, 0.423356, 614.5817336764679, 0.50973, 0.123769, 0.422156, 618.94014074292, 0.515967, 0.12596, 0.420887, 623.2996593643493, 0.522206, 0.12815, 0.419549, 627.6580664308015, 0.528444, 0.130341, 0.418142, 632.017585052231, 0.534683, 0.132534, 0.416667, 636.3771036736601, 0.54092, 0.134729, 0.415123, 640.7355107401124, 0.547157, 0.136929, 0.413511, 645.0950293615416, 0.553392, 0.139134, 0.411829, 649.4534364279938, 0.559624, 0.141346, 0.410078, 653.8129550494232, 0.565854, 0.143567, 0.408258, 658.1713621158752, 0.572081, 0.145797, 0.406369, 662.5308807373046, 0.578304, 0.148039, 0.404411, 666.8903993587342, 0.584521, 0.150294, 0.402385, 671.2488064251861, 0.590734, 0.152563, 0.40029, 675.6083250466156, 0.59694, 0.154848, 0.398125, 679.9667321130676, 0.603139, 0.157151, 0.395891, 684.3262507344971, 0.60933, 0.159474, 0.393589, 688.684657800949, 0.615513, 0.161817, 0.391219, 693.0441764223785, 0.621685, 0.164184, 0.388781, 697.403695043808, 0.627847, 0.166575, 0.386276, 701.76210211026, 0.633998, 0.168992, 0.383704, 706.1216207316895, 0.640135, 0.171438, 0.381065, 710.4800277981415, 0.64626, 0.173914, 0.378359, 714.8395464195709, 0.652369, 0.176421, 0.375586, 719.1979534860229, 0.658463, 0.178962, 0.372748, 723.5574721074523, 0.66454, 0.181539, 0.369846, 727.916990728882, 0.670599, 0.184153, 0.366879, 732.2753977953339, 0.676638, 0.186807, 0.363849, 736.6349164167635, 0.682656, 0.189501, 0.360757, 740.9933234832154, 0.688653, 0.192239, 0.357603, 745.3528421046449, 0.694627, 0.195021, 0.354388, 749.7112491710968, 0.700576, 0.197851, 0.351113, 754.0707677925263, 0.7065, 0.200728, 0.347777, 758.4302864139557, 0.712396, 0.203656, 0.344383, 762.7886934804078, 0.718264, 0.206636, 0.340931, 767.1482121018371, 0.724103, 0.20967, 0.337424, 771.5066191682893, 0.729909, 0.212759, 0.333861, 775.8661377897187, 0.735683, 0.215906, 0.330245, 780.2245448561707, 0.741423, 0.219112, 0.326576, 784.5840634776, 0.747127, 0.222378, 0.322856, 788.9424705440521, 0.752794, 0.225706, 0.319085, 793.3019891654816, 0.758422, 0.229097, 0.315266, 797.6615077869109, 0.76401, 0.232554, 0.311399, 802.0199148533629, 0.769556, 0.236077, 0.307485, 806.3794334747924, 0.775059, 0.239667, 0.303526, 810.7378405412444, 0.780517, 0.243327, 0.299523, 815.0973591626741, 0.785929, 0.247056, 0.295477, 819.4557662291259, 0.791293, 0.250856, 0.29139, 823.8152848505554, 0.796607, 0.254728, 0.287264, 828.1748034719848, 0.801871, 0.258674, 0.283099, 832.533210538437, 0.807082, 0.262692, 0.278898, 836.8927291598663, 0.812239, 0.266786, 0.274661, 841.2511362263183, 0.817341, 0.270954, 0.27039, 845.6106548477477, 0.822386, 0.275197, 0.266085, 849.9690619141999, 0.827372, 0.279517, 0.26175, 854.3285805356293, 0.832299, 0.283913, 0.257383, 858.6880991570586, 0.837165, 0.288385, 0.252988, 863.0465062235106, 0.841969, 0.292933, 0.248564, 867.4060248449401, 0.846709, 0.297559, 0.244113, 871.7644319113922, 0.851384, 0.30226, 0.239636, 876.1239505328218, 0.855992, 0.307038, 0.235133, 880.4823575992735, 0.860533, 0.311892, 0.230606, 884.841876220703, 0.865006, 0.316822, 0.226055, 889.2013948421327, 0.869409, 0.321827, 0.221482, 893.5598019085847, 0.873741, 0.326906, 0.216886, 897.919320530014, 0.878001, 0.33206, 0.212268, 902.277727596466, 0.882188, 0.337287, 0.207628, 906.6372462178956, 0.886302, 0.342586, 0.202968, 910.9956532843477, 0.890341, 0.347957, 0.198286, 915.355171905777, 0.894305, 0.353399, 0.193584, 919.7146905272062, 0.898192, 0.358911, 0.18886, 924.0730975936585, 0.902003, 0.364492, 0.184116, 928.4326162150879, 0.905735, 0.37014, 0.17935, 932.7910232815399, 0.90939, 0.375856, 0.174563, 937.1505419029694, 0.912966, 0.381636, 0.169755, 941.5089489694212, 0.916462, 0.387481, 0.164924, 945.8684675908509, 0.919879, 0.393389, 0.16007, 950.2279862122803, 0.923215, 0.399359, 0.155193, 954.5863932787323, 0.92647, 0.405389, 0.150292, 958.9459119001617, 0.929644, 0.411479, 0.145367, 963.3043189666138, 0.932737, 0.417627, 0.140417, 967.6638375880433, 0.935747, 0.423831, 0.13544, 972.0222446544952, 0.938675, 0.430091, 0.130438, 976.3817632759246, 0.941521, 0.436405, 0.125409, 980.741281897354, 0.944285, 0.442772, 0.120354, 985.0996889638062, 0.946965, 0.449191, 0.115272, 989.4592075852356, 0.949562, 0.45566, 0.110164, 993.8176146516876, 0.952075, 0.462178, 0.105031, 998.177133273117, 0.954506, 0.468744, 0.099874, 1002.5355403395691, 0.956852, 0.475356, 0.094695, 1006.8950589609985, 0.959114, 0.482014, 0.089499, 1011.2534660274506, 0.961293, 0.488716, 0.084289, 1015.61298464888, 0.963387, 0.495462, 0.079073, 1019.9725032703095, 0.965397, 0.502249, 0.073859, 1024.3309103367615, 0.967322, 0.509078, 0.068659, 1028.6904289581908, 0.969163, 0.515946, 0.063488, 1033.0488360246427, 0.970919, 0.522853, 0.058367, 1037.4083546460724, 0.97259, 0.529798, 0.053324, 1041.7667617125244, 0.974176, 0.53678, 0.048392, 1046.1262803339541, 0.975677, 0.543798, 0.043618, 1050.4857989553834, 0.977092, 0.55085, 0.03905, 1054.8442060218354, 0.978422, 0.557937, 0.034931, 1059.2037246432646, 0.979666, 0.565057, 0.031409, 1063.562131709717, 0.980824, 0.572209, 0.028508, 1067.9216503311463, 0.981895, 0.579392, 0.02625, 1072.2800573975983, 0.982881, 0.586606, 0.024661, 1076.6395760190278, 0.983779, 0.593849, 0.02377, 1080.9990946404569, 0.984591, 0.601122, 0.023606, 1085.3575017069093, 0.985315, 0.608422, 0.024202, 1089.7170203283385, 0.985952, 0.61575, 0.025592, 1094.0754273947907, 0.986502, 0.623105, 0.027814, 1098.4349460162202, 0.986964, 0.630485, 0.030908, 1102.7933530826722, 0.987337, 0.63789, 0.034916, 1107.1528717041015, 0.987622, 0.64532, 0.039886, 1111.512390325531, 0.987819, 0.652773, 0.045581, 1115.8707973919832, 0.987926, 0.66025, 0.05175, 1120.2303160134124, 0.987945, 0.667748, 0.058329, 1124.5887230798644, 0.987874, 0.675267, 0.065257, 1128.9482417012937, 0.987714, 0.682807, 0.072489, 1133.306648767746, 0.987464, 0.690366, 0.07999, 1137.6661673891754, 0.987124, 0.697944, 0.087731, 1142.0256860106047, 0.986694, 0.70554, 0.095694, 1146.3840930770566, 0.986175, 0.713153, 0.103863, 1150.7436116984863, 0.985566, 0.720782, 0.112229, 1155.1020187649385, 0.984865, 0.728427, 0.120785, 1159.4615373863676, 0.984075, 0.736087, 0.129527, 1163.81994445282, 0.983196, 0.743758, 0.138453, 1168.1794630742493, 0.982228, 0.751442, 0.147565, 1172.538981695679, 0.981173, 0.759135, 0.156863, 1176.8973887621307, 0.980032, 0.766837, 0.166353, 1181.2569073835602, 0.978806, 0.774545, 0.176037, 1185.6153144500122, 0.977497, 0.782258, 0.185923, 1189.9748330714417, 0.976108, 0.789974, 0.196018, 1194.3332401378937, 0.974638, 0.797692, 0.206332, 1198.6927587593232, 0.973088, 0.805409, 0.216877, 1203.0522773807525, 0.971468, 0.813122, 0.227658, 1207.4106844472046, 0.969783, 0.820825, 0.238686, 1211.7702030686341, 0.968041, 0.828515, 0.249972, 1216.128610135086, 0.966243, 0.836191, 0.261534, 1220.4881287565154, 0.964394, 0.843848, 0.273391, 1224.8465358229673, 0.962517, 0.851476, 0.285546, 1229.206054444397, 0.960626, 0.859069, 0.29801, 1233.564461510849, 0.95872, 0.866624, 0.31082, 1237.9239801322785, 0.956834, 0.874129, 0.323974, 1242.2834987537078, 0.954997, 0.881569, 0.337475, 1246.64190582016, 0.953215, 0.888942, 0.351369, 1251.0014244415895, 0.951546, 0.896226, 0.365627, 1255.3598315080415, 0.950018, 0.903409, 0.380271, 1259.7193501294707, 0.948683, 0.910473, 0.395289, 1264.077757195923, 0.947594, 0.917399, 0.410665, 1268.4372758173522, 0.946809, 0.924168, 0.426373, 1272.7967944387817, 0.946392, 0.930761, 0.442367, 1277.155201505234, 0.946403, 0.937159, 0.458592, 1281.5147201266632, 0.946903, 0.943348, 0.47497, 1285.8731271931151, 0.947937, 0.949318, 0.491426, 1290.2326458145446, 0.949545, 0.955063, 0.50786, 1294.5910528809968, 0.95174, 0.960587, 0.524203, 1298.950571502426, 0.954529, 0.965896, 0.540361, 1303.3100901238558, 0.957896, 0.971003, 0.556275, 1307.6684971903076, 0.961812, 0.975924, 0.571925, 1312.028015811737, 0.966249, 0.980678, 0.587206, 1316.386422878189, 0.971162, 0.985282, 0.602154, 1320.7459414996185, 0.976511, 0.989753, 0.61676, 1325.1043485660707, 0.982257, 0.994109, 0.631017, 1329.4638671875, 0.988362, 0.998364, 0.644924]
        vrLUT.NanColor = [0.0, 1.0, 0.0]
        vrLUT.ScalarRangeInitialized = 1.0

        # (Lookup table, paraview variable name)
        # Use a dictionary to index into the quantities
        lut, name = {'velocity': (vrLUT, 'Vr'),
                     'density': (densityLUT, 'Density'),
                     'temperature': (tLUT, 'T'),
                     'b': (bzLUT, 'Br'),
                     'bx': (bzLUT, 'Bx'),
                     'by': (bzLUT, 'By'),
                     'bz': (bzLUT, 'Bz')}[variable]

        # Hide the colorbar
        self.displays[self.cme].SetScalarBarVisibility(self.view, False)

        # Update all displays to be colored by this variable
        for disp in self.displays.values():
            # disp.LookupTable = lut
            pvs.ColorBy(disp, name)

        # Reshow the colorbar
        self.displays[self.cme].SetScalarBarVisibility(self.view, True)

        # restore active source
        pvs.SetActiveSource(None)
        # Render the view
        pvs.Render(self.view)
