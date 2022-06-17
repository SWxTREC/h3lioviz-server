import paraview.simple as pvs


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
        self.slice_data.SliceType.Origin = [0, 0, 0]
        if slice_type == "Plane":
            self.slice_data.SliceType.Normal = normal
        elif slice_type == "Sphere":
            self.slice_data.SliceType.Radius = radius
        else:
            raise ValueError("Can only use a Plane or Sphere slice type")
        # Now make point data on that slice
        self.slice = pvs.CellDatatoPointData(
            registrationName=f"{name}-Slice-CellDatatoPointData", Input=self.slice_data
        )
        self.slice.ProcessAllArrays = 1

        # Set up additional filters for streamlines
        # 1. Ellipse source (Circle at 0.2 AU) in the proper plane
        # 2. Calculator filter for creating the vector components
        # 3. CellData -> PointData filter on the plane
        # 4. StreamTracer with custom source from (1)
        # 5. Tubes for better display of (4)
        # 6. Arrows to indicate direction of the arrows

        # Our stream tracer source needs to have the same plane
        # as our slice, and 0.2 for the radius
        # TODO: Handle spherical sources as well?
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

        # Set up the displays
        self.slice_disp = pvs.Show(self.slice, self.view, "GeometryRepresentation")
        self.slice_disp.Representation = "Surface"
        self._variable = "Bz"
        self.slice_disp.ColorArrayName = ["POINTS", self._variable]
        self.slice_disp.LookupTable = pvs.GetColorTransferFunction(
            "Bz", self.slice_disp
        )

        # Streamlines
        self.streamlines_disp = pvs.Show(
            self.streamlines, self.view, "GeometryRepresentation"
        )
        # Add in a magnetic polarity colormap (radial in or out)
        # with two values blue/red
        # separate=True makes sure it doesn't overwrite the Br of the
        # frontend choices
        bpLUT = pvs.GetColorTransferFunction("Br", self.streamlines_disp, separate=True)
        bpLUT.RGBPoints = [-1e5, 0.5, 0.5, 0.5, 1e5, 0.9, 0.9, 0.9]
        bpLUT.ScalarRangeInitialized = 1.0
        bpLUT.NumberOfTableValues = 2
        self.streamlines_disp.Representation = "Surface"
        self.streamlines_disp.ColorArrayName = ["POINTS", "Br"]
        self.streamlines_disp.LookupTable = bpLUT

        # B-field vectors
        self.arrows_disp = pvs.Show(self.arrows, self.view, "GeometryRepresentation")
        self.arrows_disp.Representation = "Surface"
        self.arrows_disp.ColorArrayName = ["POINTS", "Br"]
        self.arrows_disp.LookupTable = bpLUT

    @property
    def variable(self):
        """The variable to use for colormapping this slice"""
        return self._variable

    @variable.setter
    def variable(self, var):
        self._variable = var
        pvs.ColorBy(self.slice_disp, var)

    def hide(self):
        """Hide the plane"""
        pvs.Hide(self.slice)

    def hide_streamlines(self):
        """Hide the streamlines on the plane"""
        pvs.Hide(self.streamlines)
        pvs.Hide(self.arrows)

    def show(self):
        """Show the plane"""
        pvs.Show(self.slice)

    def show_streamlines(self):
        """Show the streamlines on the plane"""
        pvs.Show(self.streamlines)
        pvs.Show(self.arrows)
