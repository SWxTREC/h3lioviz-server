import datetime
import math
import pathlib

import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc

import satellite

# Global definitions of variables
# Range for each lookup table
LUT_RANGE = {
    "Vr": [300, 900],
    "Density": [0, 30],
    "Pressure": [1e5, 1e7],
    "T": [1e4, 1e6],
    "Br": [-10, 10],
    "Bx": [-10, 10],
    "By": [-10, 10],
    "Bz": [-10, 10],
    "DP": [0, 1],
}

# Control points for the opacity mapping
# Can be either 2 or 3 values
# 2: Min/max opacity corresponding to the min/max data
# 3: Min, middle, max opacity corresponding to min/center/max data
OPACITY_VALUES = {
    "Vr": [0.2, 0.9],
    "Density": [0.2, 0.9],
    "Pressure": [0.2, 0.9],
    "T": [0.2, 0.9],
    "Br": [0.9, 0.2, 0.9],
    "Bx": [0.9, 0.2, 0.9],
    "By": [0.9, 0.2, 0.9],
    "Bz": [0.9, 0.2, 0.9],
    "DP": [0.2, 0.9]
}

# Default colormaps to use for the variables
DEFAULT_CMAP = {
    "Vr": "Plasma (matplotlib)",
    "Density": "Viridis (matplotlib)",
    "Pressure": "Viridis (matplotlib)",
    "T": "Inferno (matplotlib)",
    "Br": "Cool to Warm",
    "Bx": "Cool to Warm",
    "By": "Cool to Warm",
    "Bz": "Cool to Warm",
    "DP": "Plasma (matplotlib)",
}

# Name mapping from frontend to variables in dataset
VARIABLE_MAP = {
    "velocity": "Vr",
    "density": "Density",
    "pressure": "Pressure",
    "temperature": "T",
    "b": "Br",
    "bx": "Bx",
    "by": "By",
    "bz": "Bz",
    "dp": "DP",
}

VARIABLE_LABEL = {
    "velocity": "Velocity (km/s)",
    "density": "Density (r$^2$N/cm$^3$)",
    "pressure": "Ram pressure (r$^2$N/cm$^3$ * km$^2$/s$^2$)",
    "temperature": "Temperature (K)",
    "b": "Br (nT)",
    "bx": "Bx (nT)",
    "by": "By (nT)",
    "bz": "Bz (nT)",
    "dp": "Cloud tracer (-)",
}


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
        self._data_dir = pathlib.Path(dirname)

        # create a new 'NetCDF Reader' from the full data path
        fname = str(self._data_dir / "pv-data-3d.nc")
        self.celldata = pvs.NetCDFReader(
            registrationName="enlil-data", FileName=[fname]
        )
        self.celldata.Dimensions = "(longitude, latitude, radius)"

        # Force all cell data to point data in the volume
        self.data = pvs.CellDatatoPointData(
            registrationName=f"3D-CellDatatoPointData", Input=self.celldata
        )
        self.data.ProcessAllArrays = 1
        self.data.PassCellData = 1

        self.time_string = pvs.Text(registrationName="Time")
        # Don't add in any text right now
        self.time_string.Text = ""
        self._previous_time = None

        # create a new 'Threshold' to represent the CME
        self.threshold_cme = pvs.Threshold(registrationName="CME", Input=self.data)
        # We really only want a minimum value, so just set the maximum high
        self.threshold_cme.ThresholdRange = [1e-5, 1e5]
        # DP is the variable name in Enlil
        self.threshold_cme.Scalars = ["CELLS", "DP"]
        self.cme = pvs.Contour(registrationName="contoured_cme", Input=self.data)
        self.cme.ContourBy = ["POINTS", "DP"]
        self.cme.ComputeNormals = 0
        self.cme.Isosurfaces = [0.2]
        self.cme.PointMergeMethod = 'Uniform Binning'

        self.cme_contours = pvs.Contour(
            registrationName="CME-contour", Input=self.threshold_cme
        )
        self.cme_contours.ContourBy = ["POINTS", "Density"]
        self.cme_contours.Isosurfaces = []
        self.cme_contours.PointMergeMethod = "Uniform Binning"

        # Create a threshold that can be modified by the user, we use
        # two contours here instead because it looks a bit nicer.
        self.threshold = pvs.Contour(registrationName="Threshold", Input=self.data)
        self.threshold.ContourBy = ["POINTS", "Density"]
        self.threshold.Isosurfaces = [10, 50]
        self.threshold.PointMergeMethod = "Uniform Binning"

        # Create a Longitude slice
        self.lon_slice_data = pvs.Slice(
            registrationName="Longitude", Input=self.celldata
        )
        self.lon_slice_data.SliceType = "Plane"
        self.lon_slice_data.HyperTreeGridSlicer = "Plane"
        self.lon_slice_data.SliceOffsetValues = [0.0]
        self.lon_slice_data.SliceType.Origin = [0, 0, 0]
        self.lon_slice_data.SliceType.Normal = [0.0, 0.0, 1.0]
        # Now make point data on that slice
        self.lon_slice = pvs.CellDatatoPointData(
            registrationName=f"lon-slice-CellDatatoPointData", Input=self.lon_slice_data
        )
        self.lon_slice.ProcessAllArrays = 1
        self._add_streamlines("lon")

        # Create a Latitude slice
        self.lat_slice_data = pvs.Slice(
            registrationName="Latitude", Input=self.celldata
        )
        self.lat_slice_data.SliceType = "Plane"
        self.lat_slice_data.HyperTreeGridSlicer = "Plane"
        self.lat_slice_data.SliceOffsetValues = [0.0]
        self.lat_slice_data.SliceType.Origin = [0, 0, 0]
        self.lat_slice_data.SliceType.Normal = [0.0, 1.0, 0.0]
        # Now make point data on that slice
        self.lat_slice = pvs.CellDatatoPointData(
            registrationName=f"lat-slice-CellDatatoPointData", Input=self.lat_slice_data
        )
        self.lat_slice.ProcessAllArrays = 1
        self._add_streamlines("lat")

        # Dictionary mapping of string names to the object
        self.objs = {
            s: getattr(self, s)
            for s in (
                "lon_slice",
                "lat_slice",
                "cme",
                "data",
                "lon_arrows",
                "lon_streamlines",
                "lat_arrows",
                "lat_streamlines",
                "threshold",
                "cme_contours",
            )
        }
        # Initialize an empty dictionary to store the displays of the objects
        self.displays = {}
        self._setup_views()
        # After we have a view, we can start adding the satellites
        self._setup_satellites()
        self.update(None, None)

        # Call the update function every time the TimeKeeper gets modified
        # The final 1.0 is optional, but sets it as high priority to first
        # do this before other rendering.
        pvs.GetAnimationScene().TimeKeeper.AddObserver(
            "PropertyModifiedEvent", self.update, 1.0
        )

        # Set the default to have the slice go through Earth
        self.snap_solar_plane("ecliptic")

    def _setup_views(self):
        """Setup the rendering view."""
        # disable automatic camera reset on 'Show'
        pvs._DisableFirstRenderCameraReset()

        # Get the initial 'Render View'
        self.view = pvs.GetActiveView()
        self.view.ViewSize = [600, 600]
        self.view.AxesGrid = "GridAxes3DActor"
        self.view.CenterOfRotation = [0, 0, 0]
        self.view.StereoType = "Crystal Eyes"
        self.view.CameraPosition = [-3, 3, 3]
        self.view.CameraFocalPoint = [0, 0, 0]
        self.view.CameraViewUp = [0, 0, 1]
        self.view.CameraFocalDisk = 1.0
        self.view.CameraParallelScale = 2
        self.view.BackEnd = "OSPRay raycaster"
        self.view.OSPRayMaterialLibrary = pvs.GetMaterialLibrary()

        # Time string
        disp = pvs.Show(self.time_string, self.view, "TextSourceRepresentation")

        # get color transfer function/color map for Bz initially
        bzLUT = pvs.GetColorTransferFunction("Bz")
        bzLUT.RGBPoints = [
            -10,
            0.231373,
            0.298039,
            0.752941,
            0,
            0.865003,
            0.865003,
            0.865003,
            10,
            0.705882,
            0.0156863,
            0.14902,
        ]
        bzLUT.ScalarRangeInitialized = 1.0
        # get opacity transfer function/opacity map for 'Bz'
        bzPWF = pvs.GetOpacityTransferFunction("Bz")
        bzPWF.Points = [-10, 0.0, 0.5, 0.0, 10, 1.0, 0.5, 0.0]
        bzPWF.ScalarRangeInitialized = 1

        # CME Threshold
        disp = pvs.Show(self.cme, self.view, "GeometryRepresentation")
        self.displays[self.cme] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = [None, ""]
        disp.Opacity = 0.25

        # CME Contours
        disp = pvs.Show(self.cme_contours, self.view, "GeometryRepresentation")
        self.displays[self.cme_contours] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = bzLUT

        disp = pvs.Show(self.threshold, self.view, "GeometryRepresentation")
        self.displays[self.threshold] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = bzLUT

        # Latitude
        disp = pvs.Show(self.lat_slice, self.view, "GeometryRepresentation")
        self.displays[self.lat_slice] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = bzLUT

        # Longitude
        disp = pvs.Show(self.lon_slice, self.view, "GeometryRepresentation")
        self.displays[self.lon_slice] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = bzLUT

        # Streamlines
        disp = pvs.Show(self.lon_streamlines, self.view, "GeometryRepresentation")
        # Add in a magnetic polarity colormap (radial in or out)
        # with two values blue/red
        # separate=True makes sure it doesn't overwrite the Br of the
        # frontend choices
        bpLUT = pvs.GetColorTransferFunction("Br", disp, separate=True)
        bpLUT.RGBPoints = [-1e5, 0.5, 0.5, 0.5, 1e5, 0.9, 0.9, 0.9]
        bpLUT.ScalarRangeInitialized = 1.0
        bpLUT.NumberOfTableValues = 2
        self.displays[self.lon_streamlines] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = bpLUT

        # B-field vectors
        disp = pvs.Show(self.lon_arrows, self.view, "GeometryRepresentation")
        self.displays[self.lon_arrows] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = bpLUT

        # latitudinal plane
        disp = pvs.Show(self.lat_slice, self.view, "GeometryRepresentation")
        self.displays[self.lat_slice] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Bz"]
        disp.LookupTable = bzLUT

        # Streamlines
        disp = pvs.Show(self.lat_streamlines, self.view, "GeometryRepresentation")
        # Add in a magnetic polarity colormap (radial in or out)
        # with two values blue/red
        # separate=True makes sure it doesn't overwrite the Br of the
        # frontend choices
        bpLUT = pvs.GetColorTransferFunction("Br", disp, separate=True)
        bpLUT.RGBPoints = [-1e5, 0.5, 0.5, 0.5, 1e5, 0.9, 0.9, 0.9]
        bpLUT.ScalarRangeInitialized = 1.0
        bpLUT.NumberOfTableValues = 2
        self.displays[self.lat_streamlines] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = bpLUT

        # B-field vectors
        disp = pvs.Show(self.lat_arrows, self.view, "GeometryRepresentation")
        self.displays[self.lat_arrows] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", "Br"]
        disp.LookupTable = bpLUT

        # Set colormaps
        for name in VARIABLE_MAP:
            self.set_colormap(name)

        # hide this data from the default initial view
        for x in [
            self.lon_slice,
            self.lon_arrows,
            self.lon_streamlines,
            self.lat_arrows,
            self.lat_streamlines,
            self.threshold,
            self.cme_contours,
        ]:
            pvs.Hide(x, self.view)

        # restore active source
        pvs.SetActiveSource(None)

    def _setup_satellites(self):
        """
        Initializes the satellites locations and plots them as spheres.
        """
        sats = list(self._data_dir.glob("*.json"))
        # strip evo.name.json to only keep "name"
        # TODO: Use the other satellites eventually
        #       For now, we are only using the stereo spacecraft
        # self.satellites = {x.name[4:-5]: satellite.Satellite(x.name[4:-5], x, view=self.view) for x in sats if "earth" not in x.name}
        self.satellites = {
            x.name[4:-5]: satellite.Satellite(x.name[4:-5], x, view=self.view)
            for x in sats
            if "stereo" in x.name
        }
        # Filter for Earth in the name
        self.earth = satellite.Earth(
            list(filter(lambda x: "earth" in x.name, sats))[0], view=self.view
        )
        self.sun = satellite.Sun(self._data_dir / "solar_images", view=self.view)

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
        curr_slice_data = getattr(self, f"{plane}_slice_data")
        stream_source = pvs.Ellipse(registrationName=f"{plane}-StreamSource")
        stream_source.Center = [0.0, 0.0, 0.0]
        stream_source.Normal = curr_slice_data.SliceType.Normal
        radius = [0.2, 0, 0] if plane == "lon" else [0, 0, 0.2]
        stream_source.MajorRadiusVector = radius
        # Controls how many streamlines we have
        stream_source.Resolution = 50

        # Create the magnetic field vectors through a PV Function
        curr_slice = getattr(self, f"{plane}_slice")
        bvec = pvs.Calculator(registrationName=f"{plane}-Bvec", Input=curr_slice)
        bvec.AttributeType = "Point Data"
        bvec.ResultArrayName = "Bvec"
        bvec.Function = "Bx*iHat + By*jHat + Bz*kHat"

        stream_input = pvs.StreamTracerWithCustomSource(
            registrationName=f"{plane}-StreamTracerWithCustomSource",
            Input=bvec,
            SeedSource=stream_source,
        )
        stream_input.Vectors = ["POINTS", "Bvec"]
        stream_input.SurfaceStreamlines = 1
        stream_input.MaximumStreamlineLength = 3.4
        stream_input.ComputeVorticity = 0

        streamlines = pvs.Tube(
            registrationName=f"{plane}-Streamlines", Input=stream_input
        )
        streamlines.Capping = 1
        streamlines.Radius = 0.005

        # create a new 'Glyph' in the slice (Arrow/vectors)
        arrows = pvs.Glyph(
            registrationName=f"{plane}-B-Arrows", Input=stream_input, GlyphType="Cone"
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
        base_dir = self._data_dir / ".."
        # Now search to see if there is a pv-data-3d.nc in the directories
        # and if not, ignore that entry
        dirs = [x for x in base_dir.iterdir() if (x / "pv-data-3d.nc").exists()]
        return dirs

    @exportRpc("pv.enlil.get_variable_range")
    def get_variable_range(self, name):
        """
        Get the range of values for a variable at the current timestep.

        name : str
            Name of variable to colormap all of the surfaces by
        """
        variable = VARIABLE_MAP[name]
        return self.celldata.CellData.GetArray(variable).GetRange()

    @exportRpc("pv.enlil.directory")
    def update_dataset(self, dirname):
        """
        Change the dataset directory to the one specified by dirname

        dirname : str
            Path to the dataset file (/dirname/pv-data-3d.nc)
        """
        self._data_dir = pathlib.Path(dirname)
        # Update the evolution files associated with the run
        # NOTE: We need to delete the evolutions first, there must be a
        #       dangling reference within the cpp that causes a segfault
        #       if we just update the object dictionary without removal
        del self.satellites
        del self.earth
        del self.sun
        self._setup_satellites()
        # Update the primary data 3D data file
        self.data.FileName = str(dirname / "pv-data-3d.nc")
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
            elif obj == "lat_streamlines":
                # We also want to turn on the vectors
                pvs.Show(self.objs["lat_arrows"], self.view)
        elif visibility == "off":
            pvs.Hide(self.objs[obj], self.view)
            if obj == "lon_streamlines":
                # We also want to turn off the vectors
                pvs.Hide(self.objs["lon_arrows"], self.view)
            elif obj == "lat_streamlines":
                # We also want to turn off the vectors
                pvs.Hide(self.objs["lat_arrows"], self.view)
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
        label = VARIABLE_LABEL[name]

        # Update all displays to be colored by this variable
        for obj, disp in self.displays.items():
            if obj in (
                self.lon_arrows,
                self.lon_streamlines,
                self.lat_arrows,
                self.lat_streamlines,
                self.cme,
            ):
                # We don't want to update the longitude arrow colors
                continue
            pvs.ColorBy(disp, variable)

            # Also update the colorbar orientation information
            if disp.LookupTable is None:
                continue
            cbar = pvs.GetScalarBar(disp.LookupTable, self.view)
            cbar.AutoOrient = 0
            cbar.Orientation = "Horizontal"
            cbar.TextPosition = "Ticks left/bottom, annotations right/top"
            cbar.Title = label
            cbar.ComponentTitle = ""
            # Disables the endpoints which can be formatted differently
            cbar.AddRangeLabels = 0

        self.update_opacity(variable)
        self.update_lut(variable)

        # hides old scalarbars that aren't in the view and
        # shows the new variable we are using now
        pvs.UpdateScalarBars(self.view)
        # But we want to hide the streamlines colorbar
        disp = self.displays[self.lon_streamlines]
        disp.SetScalarBarVisibility(self.view, False)
        disp = self.displays[self.lat_streamlines]
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
        lut.AutomaticRescaleRangeMode = "Never"

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
            opacity_map.Points = [
                data_range[0],
                points[0],
                0.5,
                0,
                data_range[1],
                points[1],
                0.5,
                0,
            ]
        elif len(points) == 3:
            # min/mid/max
            opacity_map.Points = [
                data_range[0],
                points[0],
                0.5,
                0,
                (data_range[0] + data_range[1]) / 2,
                points[1],
                0.5,
                0,
                data_range[1],
                points[2],
                0.5,
                0,
            ]
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
        self.threshold.ContourBy = ["POINTS", variable]
        self.threshold.Isosurfaces = range

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
        self.cme_contours.ContourBy = ["POINTS", variable]
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
            loc = [-self.earth.sat.Center[2], 0, self.earth.sat.Center[0]]
        elif clip == "equator":
            # Standard coordinates, so the plane is purely in xy and
            # z is perpendicular
            loc = [0, 0, 1]
        else:
            raise ValueError(
                "The snapping clip plane must be either " '"ecliptic" or "equator"'
            )

        self.lon_slice_data.SliceType.Normal = loc
        # Also update the stream source so they stay in-sync
        self.lon_stream_source.Normal = loc

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
            loc = [-y, 0, -x]
            self.lon_slice_data.SliceType.Normal = loc
            self.lon_stream_source.Normal = loc
        elif plane == "lat":
            loc = [-y, x, 0]
            self.lat_slice_data.SliceType.Normal = loc
            self.lat_stream_source.Normal = loc
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
            raise ValueError(
                'Invalid string, only "ecliptic", "meridional", '
                'and "initial" are allowed.'
            )
        # Force the focal point to be the sun
        self.view.CameraFocalPoint = [0, 0, 0]

    @exportRpc("pv.enlil.toggle_satellites")
    def toggle_satellites(self, visibility):
        """
        Toggles the visibility of the satellites on/off from the view

        visibility : str ("on", "off")
            What to set the visibility to
        """
        if visibility == "on":
            hide_show = "show"
        elif visibility == "off":
            hide_show = "hide"
        else:
            return ["Visibility can only be 'on' or 'off'"]
        for sat in self.satellites:
            # Call the hide() or show() method
            getattr(self.satellites[sat], hide_show)()

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
        if sat == "earth":
            return self.earth.get_times()
        return self.satellites[sat].get_times()

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
        if sat == "earth":
            return self.earth.evolution.as_latis()
        return self.satellites[sat].evolution.as_latis()

    def update(self, caller, event):
        """
        Update function to call every time the time variable has changed.
        """
        curr_time = self.get_current_time()
        if self._previous_time == curr_time:
            # This is the same timestep, so we don't need to
            # update anything now
            return
        self.time_string.Text = curr_time.strftime("%Y-%m-%d %H:00")

        for x in self.satellites:
            # Update the satellite positions based on the evolution data
            self.satellites[x].update(curr_time)

        # Update the solar image and rotate the Earth image
        self.sun.update(curr_time)
        self.earth.update(curr_time)
        self._previous_time == curr_time

    def get_current_time(self):
        """Retrieves the current time of the view.

        Returns
        -------
        datetime of the current timestep
        """
        pv_time = pvs.GetAnimationScene().TimeKeeper.Time
        # The internal time variable on the ViewTime attribute is stored as
        # seconds from 1970-01-01, so we use that epoch directly internally.
        return datetime.datetime(1970, 1, 1) + datetime.timedelta(seconds=pv_time)
