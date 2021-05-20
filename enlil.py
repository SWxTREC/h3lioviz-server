import paraview.simple as pvs


class EnlilDataset:

    def __init__(self, fname):
        """Enlil 4D dataset representation in Paraview.

        This class will read in the given NetCDF file that contains
        the Enlil output. It is designed to enable an easy storage
        and access layer to the data.
        """
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

        # TODO: Show the base dataset?
        # pvs.Show(self.data, self.view, 'StructuredGridRepresentation')

        # get color transfer function/color map for 'Bz'
        bzLUT = pvs.GetColorTransferFunction('Bz')
        bzLUT.RGBPoints = [-10.103191375732422, 0.231373, 0.298039, 0.752941,
                           1.3882222175598145, 0.865003, 0.865003, 0.865003,
                           12.87963581085205, 0.705882, 0.0156863, 0.14902]
        bzLUT.ScalarRangeInitialized = 1.0

        # get opacity transfer function/opacity map for 'Bz'
        bzPWF = pvs.GetOpacityTransferFunction('Bz')
        bzPWF.Points = [-10.103191375732422, 0.0, 0.5, 0.0, 12.87963581085205,
                        1.0, 0.5, 0.0]
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
