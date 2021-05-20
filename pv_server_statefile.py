
# import to process args
import os

# import paraview modules.
from paraview.web import pv_wslink
from paraview.web import protocols as pv_protocols

import paraview
from paraview import simple
from wslink import server

import argparse

# =============================================================================
# Create custom PVServerProtocol class to handle clients requests
# =============================================================================


class _DemoServer(pv_wslink.PVServerProtocol):
    authKey = "wslink-secret"

    def initialize(self):
        # Bring used components
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebMouseHandler())
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebViewPort())
        self.registerVtkWebProtocol(
            pv_protocols.ParaViewWebPublishImageDelivery(decode=False))
        self.updateSecret(_DemoServer.authKey)

        # tell the C++ web app to use no encoding.
        # ParaViewWebPublishImageDelivery must be set to decode=False to match.
        self.getApplication().SetImageEncoding(0)

        # Disable interactor-based render calls
        simple.GetRenderView().EnableRenderOnInteraction = 0
        simple.GetRenderView().Background = [38, 55, 90]

        # The NetCDF file with the data
        fname = 'data/test_xarray.nc'
        read_statefile(fname)


def read_statefile(fname):
    # state file generated using paraview version 5.9.1

    #### disable automatic camera reset on 'Show'
    simple._DisableFirstRenderCameraReset()

    # ----------------------------------------------------------------
    # setup views used in the visualization
    # ----------------------------------------------------------------

    # get the material library
    materialLibrary1 = simple.GetMaterialLibrary()

    # Create a new 'Render View'
    renderView1 = simple.CreateView('RenderView')
    # renderView1.Background = [38, 55, 90]
    renderView1.ViewSize = [1008, 539]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.CenterOfRotation = [0.4249999523162842, 0.0, 5.960464477539063e-08]
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraPosition = [-8.183350055920059, 5.466998450462119, 6.364584945029508]
    renderView1.CameraFocalPoint = [0.4249999523162841, 1.0256390200333957e-15, 5.960464525220887e-08]
    renderView1.CameraViewUp = [0.446016916276754, -0.28531285472026074, 0.8483309998616994]
    renderView1.CameraFocalDisk = 1.0
    renderView1.CameraParallelScale = 2.1250001580766877
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.OSPRayMaterialLibrary = materialLibrary1

    simple.SetActiveView(None)

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = simple.CreateLayout(name='Layout #1')
    layout1.AssignView(0, renderView1)
    layout1.SetSize(1008, 539)

    # ----------------------------------------------------------------
    # restore active view
    simple.SetActiveView(renderView1)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'NetCDF Reader'
    test_xarraync = simple.NetCDFReader(registrationName='test_xarray.nc', FileName=[fname])
    test_xarraync.Dimensions = '(longitude, latitude, radius)'

    # create a new 'Calculator'
    bvec = simple.Calculator(registrationName='Bvec', Input=test_xarraync)
    bvec.AttributeType = 'Cell Data'
    bvec.ResultArrayName = 'Bvec'
    bvec.Function = 'Bx*iHat + By*jHat + Bz*kHat'

    # create a new 'Slice'
    longitude = simple.Slice(registrationName='Longitude', Input=bvec)
    longitude.SliceType = 'Plane'
    longitude.HyperTreeGridSlicer = 'Plane'
    longitude.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    longitude.SliceType.Origin = [3.008704396734174e-14, 2.8199664825478976e-14, 5.659217572340225e-08]
    longitude.SliceType.Normal = [0.0, 0.0, 1.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    longitude.HyperTreeGridSlicer.Origin = [3.008704396734174e-14, 2.8199664825478976e-14, 5.659217572340225e-08]

    # create a new 'Stream Tracer'
    streamlines = simple.StreamTracer(registrationName='Streamlines', Input=longitude,
        SeedType='Point Cloud')
    streamlines.Vectors = ['CELLS', 'Bvec']
    streamlines.SurfaceStreamlines = 1
    streamlines.MaximumStreamlineLength = 3.40000014164759

    # init the 'Point Cloud' selected for 'SeedType'
    streamlines.SeedType.Center = [3.008704396734174e-14, 2.8199664825478976e-14, 5.659217537125256e-08]
    streamlines.SeedType.NumberOfPoints = 500
    streamlines.SeedType.Radius = 0.340000014164759

    # create a new 'Threshold'
    cME = simple.Threshold(registrationName='CME', Input=test_xarraync)
    cME.Scalars = ['CELLS', 'DP']
    cME.ThresholdRange = [0.001, 15.765473365783691]

    # create a new 'Slice'
    latitude = simple.Slice(registrationName='Latitude', Input=bvec)
    latitude.SliceType = 'Plane'
    latitude.HyperTreeGridSlicer = 'Plane'
    latitude.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    latitude.SliceType.Origin = [3.008704396734174e-14, 2.8199664825478976e-14, 5.659217572340225e-08]
    latitude.SliceType.Normal = [0.0, 1.0, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    latitude.HyperTreeGridSlicer.Origin = [3.008704396734174e-14, 2.8199664825478976e-14, 5.659217572340225e-08]

    # create a new 'Glyph'
    arrows = simple.Glyph(registrationName='Arrows', Input=longitude,
        GlyphType='Arrow')
    arrows.OrientationArray = ['CELLS', 'Bvec']
    arrows.ScaleArray = ['POINTS', 'No scale array']
    arrows.ScaleFactor = 0.340000014164759
    arrows.GlyphTransform = 'Transform2'
    arrows.MaximumNumberOfSamplePoints = 500

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from test_xarraync
    test_xarrayncDisplay = simple.Show(test_xarraync, renderView1, 'StructuredGridRepresentation')

    # get color transfer function/color map for 'Bz'
    bzLUT = simple.GetColorTransferFunction('Bz')
    bzLUT.RGBPoints = [-10.103191375732422, 0.231373, 0.298039, 0.752941, 1.3882222175598145, 0.865003, 0.865003, 0.865003, 12.87963581085205, 0.705882, 0.0156863, 0.14902]
    bzLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'Bz'
    bzPWF = simple.GetOpacityTransferFunction('Bz')
    bzPWF.Points = [-10.103191375732422, 0.0, 0.5, 0.0, 12.87963581085205, 1.0, 0.5, 0.0]
    bzPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    test_xarrayncDisplay.Representation = 'Surface'
    test_xarrayncDisplay.ColorArrayName = ['CELLS', 'Bz']
    test_xarrayncDisplay.LookupTable = bzLUT
    test_xarrayncDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    test_xarrayncDisplay.SelectOrientationVectors = 'None'
    test_xarrayncDisplay.ScaleFactor = 0.3400000143191053
    test_xarrayncDisplay.SelectScaleArray = 'None'
    test_xarrayncDisplay.GlyphType = 'Arrow'
    test_xarrayncDisplay.GlyphTableIndexArray = 'None'
    test_xarrayncDisplay.GaussianRadius = 0.017000000715955265
    test_xarrayncDisplay.SetScaleArray = [None, '']
    test_xarrayncDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    test_xarrayncDisplay.OpacityArray = [None, '']
    test_xarrayncDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    test_xarrayncDisplay.DataAxesGrid = 'GridAxesRepresentation'
    test_xarrayncDisplay.PolarAxes = 'PolarAxesRepresentation'
    test_xarrayncDisplay.ScalarOpacityFunction = bzPWF
    test_xarrayncDisplay.ScalarOpacityUnitDistance = 0.03188458069611

    # show data from cME
    cMEDisplay = simple.Show(cME, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    cMEDisplay.Representation = 'Surface'
    cMEDisplay.ColorArrayName = ['CELLS', 'Bz']
    cMEDisplay.LookupTable = bzLUT
    cMEDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    cMEDisplay.SelectOrientationVectors = 'None'
    cMEDisplay.ScaleFactor = 0.09197479853610144
    cMEDisplay.SelectScaleArray = 'None'
    cMEDisplay.GlyphType = 'Arrow'
    cMEDisplay.GlyphTableIndexArray = 'None'
    cMEDisplay.GaussianRadius = 0.004598739926805072
    cMEDisplay.SetScaleArray = [None, '']
    cMEDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    cMEDisplay.OpacityArray = [None, '']
    cMEDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    cMEDisplay.DataAxesGrid = 'GridAxesRepresentation'
    cMEDisplay.PolarAxes = 'PolarAxesRepresentation'
    cMEDisplay.ScalarOpacityFunction = bzPWF
    cMEDisplay.ScalarOpacityUnitDistance = 0.02090409368521722
    cMEDisplay.OpacityArrayName = [None, '']

    # show data from bvec
    bvecDisplay = simple.Show(bvec, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    bvecDisplay.Representation = 'Outline'
    bvecDisplay.ColorArrayName = ['CELLS', 'Bz']
    bvecDisplay.LookupTable = bzLUT
    bvecDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    bvecDisplay.SelectOrientationVectors = 'Bvec'
    bvecDisplay.ScaleFactor = 0.3400000143191053
    bvecDisplay.SelectScaleArray = 'None'
    bvecDisplay.GlyphType = 'Arrow'
    bvecDisplay.GlyphTableIndexArray = 'None'
    bvecDisplay.GaussianRadius = 0.017000000715955265
    bvecDisplay.SetScaleArray = [None, '']
    bvecDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    bvecDisplay.OpacityArray = [None, '']
    bvecDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    bvecDisplay.DataAxesGrid = 'GridAxesRepresentation'
    bvecDisplay.PolarAxes = 'PolarAxesRepresentation'
    bvecDisplay.ScalarOpacityFunction = bzPWF
    bvecDisplay.ScalarOpacityUnitDistance = 0.03188458069611

    # show data from latitude
    latitudeDisplay = simple.Show(latitude, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    latitudeDisplay.Representation = 'Surface'
    latitudeDisplay.ColorArrayName = ['CELLS', 'Bz']
    latitudeDisplay.LookupTable = bzLUT
    latitudeDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    latitudeDisplay.SelectOrientationVectors = 'Bvec'
    latitudeDisplay.ScaleFactor = 0.2944486496874232
    latitudeDisplay.SelectScaleArray = 'None'
    latitudeDisplay.GlyphType = 'Arrow'
    latitudeDisplay.GlyphTableIndexArray = 'None'
    latitudeDisplay.GaussianRadius = 0.01472243248437116
    latitudeDisplay.SetScaleArray = [None, '']
    latitudeDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    latitudeDisplay.OpacityArray = [None, '']
    latitudeDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    latitudeDisplay.DataAxesGrid = 'GridAxesRepresentation'
    latitudeDisplay.PolarAxes = 'PolarAxesRepresentation'

    # show data from longitude
    longitudeDisplay = simple.Show(longitude, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    longitudeDisplay.Representation = 'Surface'
    longitudeDisplay.ColorArrayName = ['CELLS', 'Bz']
    longitudeDisplay.LookupTable = bzLUT
    longitudeDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    longitudeDisplay.SelectOrientationVectors = 'Bvec'
    longitudeDisplay.ScaleFactor = 0.340000014164759
    longitudeDisplay.SelectScaleArray = 'None'
    longitudeDisplay.GlyphType = 'Arrow'
    longitudeDisplay.GlyphTableIndexArray = 'None'
    longitudeDisplay.GaussianRadius = 0.017000000708237952
    longitudeDisplay.SetScaleArray = [None, '']
    longitudeDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    longitudeDisplay.OpacityArray = [None, '']
    longitudeDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    longitudeDisplay.DataAxesGrid = 'GridAxesRepresentation'
    longitudeDisplay.PolarAxes = 'PolarAxesRepresentation'

    # show data from streamlines
    streamlinesDisplay = simple.Show(streamlines, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    streamlinesDisplay.Representation = 'Surface'
    streamlinesDisplay.ColorArrayName = [None, '']
    streamlinesDisplay.OSPRayScaleArray = 'AngularVelocity'
    streamlinesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    streamlinesDisplay.SelectOrientationVectors = 'Normals'
    streamlinesDisplay.ScaleFactor = 0.3336487650871277
    streamlinesDisplay.SelectScaleArray = 'AngularVelocity'
    streamlinesDisplay.GlyphType = 'Arrow'
    streamlinesDisplay.GlyphTableIndexArray = 'AngularVelocity'
    streamlinesDisplay.GaussianRadius = 0.016682438254356384
    streamlinesDisplay.SetScaleArray = ['POINTS', 'AngularVelocity']
    streamlinesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    streamlinesDisplay.OpacityArray = ['POINTS', 'AngularVelocity']
    streamlinesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    streamlinesDisplay.DataAxesGrid = 'GridAxesRepresentation'
    streamlinesDisplay.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    streamlinesDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    streamlinesDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # show data from arrows
    arrowsDisplay = simple.Show(arrows, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    arrowsDisplay.Representation = 'Surface'
    arrowsDisplay.ColorArrayName = ['POINTS', 'Bz']
    arrowsDisplay.LookupTable = bzLUT
    arrowsDisplay.OSPRayScaleArray = 'BP'
    arrowsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    arrowsDisplay.SelectOrientationVectors = 'BP'
    arrowsDisplay.ScaleFactor = 0.3804943442344666
    arrowsDisplay.SelectScaleArray = 'BP'
    arrowsDisplay.GlyphType = 'Arrow'
    arrowsDisplay.GlyphTableIndexArray = 'BP'
    arrowsDisplay.GaussianRadius = 0.019024717211723326
    arrowsDisplay.SetScaleArray = ['POINTS', 'BP']
    arrowsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    arrowsDisplay.OpacityArray = ['POINTS', 'BP']
    arrowsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    arrowsDisplay.DataAxesGrid = 'GridAxesRepresentation'
    arrowsDisplay.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    arrowsDisplay.ScaleTransferFunction.Points = [-96.6900405883789, 0.0, 0.5, 0.0, 97.67322540283203, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    arrowsDisplay.OpacityTransferFunction.Points = [-96.6900405883789, 0.0, 0.5, 0.0, 97.67322540283203, 1.0, 0.5, 0.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for bzLUT in view renderView1
    bzLUTColorBar = simple.GetScalarBar(bzLUT, renderView1)
    bzLUTColorBar.Title = 'Bz'
    bzLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    bzLUTColorBar.Visibility = 1

    # show color legend
    test_xarrayncDisplay.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    simple.Hide(test_xarraync, renderView1)

    # show color legend
    cMEDisplay.SetScalarBarVisibility(renderView1, True)

    # show color legend
    bvecDisplay.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    simple.Hide(bvec, renderView1)

    # show color legend
    latitudeDisplay.SetScalarBarVisibility(renderView1, True)

    # show color legend
    longitudeDisplay.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    simple.Hide(longitude, renderView1)

    # show color legend
    arrowsDisplay.SetScalarBarVisibility(renderView1, True)

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # restore active source
    simple.SetActiveSource(bvec)
    # ----------------------------------------------------------------


# =============================================================================
# Main: Parse args and start server
# =============================================================================


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="ParaViewWeb Demo")

    # Add default arguments
    server.add_arguments(parser)

    # Extract arguments
    args = parser.parse_args()

    # Start server
    server.start_webserver(options=args, protocol=_DemoServer)
