import paraview.simple as pvs

from display import LazyDisplay, color_by


class Slice:
    """
    A slice through the 3D domain

    Parameters
    ----------
    celldata : Cell Data
        The full 3D Paraview dataset
    slice_type : str
        The type of slice ("Plane", "Sphere")
    normal : 3-tuple
        The normal to the plane (only used for "Plane" slice_type)
    radius : Number
        The radius of the sphere (only used for "Sphere" slice_type)
    name : str
        S
    """

    def __init__(
        self,
        celldata,
        slice_type="Plane",
        normal=(0, 0, 1),
        radius=1,
        name="",
        view=None,
    ):
        self.celldata = celldata
        self.name = name
        self.view = view
        # Create the Paraview Slice
        self.slice_data = pvs.Slice(
            registrationName=f"{name}-Slice", Input=self.celldata
        )
        self.slice_data.SliceType = slice_type
        self.slice_data.HyperTreeGridSlicer = "Plane"
        self.slice_data.SliceOffsetValues = [0.0]
        if slice_type == "Plane":
            self.slice_data.SliceType.Origin = [0, 0, 0]
            self.slice_data.SliceType.Normal = normal
        elif slice_type == "Sphere":
            self.slice_data.SliceType.Center = [0, 0, 0]
            self.slice_data.SliceType.Radius = radius
        else:
            raise ValueError("Can only use a Plane or Sphere slice type")
        # Now make point data on that slice
        self.slice = pvs.CellDatatoPointData(
            registrationName=f"{name}-Slice-CellDatatoPointData", Input=self.slice_data
        )
        self.slice.ProcessAllArrays = 1

        # Set up the display, created when the slice is first shown
        self._variable = "Bz"
        self.slice_display = LazyDisplay(
            self.slice, self.view, setup=self._setup_slice_display
        )
        self.streamlines_display = None
        self.arrows_display = None

        # TODO: Handle spherical sources as well?
        #       Trying out spherical sources didn't look great at first glance
        if slice_type == "Sphere":
            return

        # Set up additional filters for streamlines
        # 1. Ellipse source (Circle at 0.2 AU) in the proper plane
        # 2. Calculator filter for creating the vector components
        # 3. CellData -> PointData filter on the plane
        # 4. StreamTracer with custom source from (1)
        # 5. Tubes for better display of (4)
        # 6. Arrows to indicate direction of the arrows

        # Our stream tracer source needs to have the same plane
        # as our slice, and 0.2 for the radius
        stream_source = pvs.Ellipse(registrationName=f"{self.name}-StreamSource")
        stream_source.Center = [0.0, 0.0, 0.0]
        stream_source.Normal = self.slice_data.SliceType.Normal
        # TODO: Look into whether this needs to "tilt" with the normal
        radius = [0.2, 0, 0] if name == "Longitude" else [0, 0, 0.2]
        stream_source.MajorRadiusVector = radius
        # Controls how many streamlines we have
        stream_source.Resolution = 50

        # Create the magnetic field vectors through a PV Function, which
        # we want to be in Point Data, not Cell Data
        bvec = pvs.Calculator(registrationName=f"{self.name}-Bvec", Input=self.slice)
        bvec.AttributeType = "Point Data"
        bvec.ResultArrayName = "Bvec"
        bvec.Function = "Bx*iHat + By*jHat + Bz*kHat"

        stream_input = pvs.StreamTracerWithCustomSource(
            registrationName=f"{self.name}-StreamTracerWithCustomSource",
            Input=bvec,
            SeedSource=stream_source,
        )
        stream_input.Vectors = ["POINTS", "Bvec"]
        stream_input.SurfaceStreamlines = 1
        stream_input.MaximumStreamlineLength = 3.4
        stream_input.ComputeVorticity = 0

        streamlines = pvs.Tube(
            registrationName=f"{self.name}-Streamlines", Input=stream_input
        )
        streamlines.Capping = 1
        streamlines.Radius = 0.005

        # create a new 'Glyph' in the slice (Arrow/vectors)
        arrows = pvs.Glyph(
            registrationName=f"{self.name}-B-Arrows",
            Input=stream_input,
            GlyphType="Cone",
        )
        arrows.OrientationArray = ["POINTS", "Bvec"]
        arrows.ScaleArray = ["POINTS", "No scale array"]
        arrows.ScaleFactor = 1
        arrows.GlyphType.Resolution = 60
        arrows.GlyphType.Radius = 0.02
        arrows.GlyphType.Height = 0.08
        arrows.GlyphTransform = "Transform2"
        arrows.GlyphMode = "Every Nth Point"
        arrows.GlyphMode = "Uniform Spatial Distribution (Bounds Based)"
        arrows.MaximumNumberOfSamplePoints = 100
        # Store the objects we need for later as attributes on self
        self.stream_source = stream_source
        self.streamlines = streamlines
        self.arrows = arrows

        # Streamlines and B-field vectors, created when first shown
        self._bpLUT = None
        self.streamlines_display = LazyDisplay(
            self.streamlines, self.view, setup=self._setup_streamlines_display
        )
        self.arrows_display = LazyDisplay(
            self.arrows, self.view, setup=self._setup_arrows_display
        )

    def _setup_slice_display(self, disp):
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = pvs.GetColorTransferFunction("Bz", disp)
        if self._variable != "Bz":
            # The coloring a display shown from the start would have now
            color_by(disp, self._variable)

    def _setup_streamlines_display(self, disp):
        # Add in a magnetic polarity colormap (radial in or out)
        # with two values blue/red
        # separate=True makes sure it doesn't overwrite the Br of the
        # frontend choices
        bpLUT = pvs.GetColorTransferFunction("Br", disp, separate=True)
        bpLUT.RGBPoints = [-1e5, 0.5, 0.5, 0.5, 1e5, 0.9, 0.9, 0.9]
        bpLUT.ScalarRangeInitialized = 1.0
        bpLUT.NumberOfTableValues = 2
        self._bpLUT = bpLUT
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = bpLUT

    def _setup_arrows_display(self, disp):
        # Uses the streamlines' colormap, so they are always shown first
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = self._bpLUT

    @property
    def variable(self):
        """The variable to use for colormapping this slice"""
        return self._variable

    @variable.setter
    def variable(self, var):
        self._variable = var
        if self.slice_display.created:
            pvs.ColorBy(self.slice_display.display, var)

    def hide(self):
        """Hide the plane"""
        self.slice_display.hide()

    def hide_streamlines(self):
        """Hide the streamlines on the plane"""
        if self.streamlines_display is None:
            # A spherical slice has no streamlines
            return
        self.streamlines_display.hide()
        self.arrows_display.hide()

    def show(self):
        """Show the plane"""
        self.slice_display.show()

    def show_streamlines(self):
        """Show the streamlines on the plane"""
        if self.streamlines_display is None:
            # A spherical slice has no streamlines
            return
        # The streamlines first: their setup makes the colormap the arrows use
        self.streamlines_display.show()
        self.arrows_display.show()
