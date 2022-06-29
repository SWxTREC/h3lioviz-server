import datetime
import math
import pathlib

import paraview.simple as pvs
from paraview.web import protocols as pv_protocols
from wslink import register as exportRpc

import models
import satellite
import slice

# Global definitions of variables
# Range for each lookup table
LUT_RANGE = {
    "velocity": [300, 900],
    "density": [0, 30],
    "pressure": [1e5, 1e7],
    "temperature": [1e4, 1e6],
    "b": [-10, 10],
    "bx": [-10, 10],
    "by": [-10, 10],
    "bz": [-10, 10],
    "dp": [0, 1],
}

# Control points for the opacity mapping
# Can be either 2 or 3 values
# 2: Min/max opacity corresponding to the min/max data
# 3: Min, middle, max opacity corresponding to min/center/max data
OPACITY_VALUES = {
    "velocity": [0.2, 0.9],
    "density": [0.2, 0.9],
    "pressure": [0.2, 0.9],
    "temperature": [0.2, 0.9],
    "b": [0.9, 0.2, 0.9],
    "bx": [0.9, 0.2, 0.9],
    "by": [0.9, 0.2, 0.9],
    "bz": [0.9, 0.2, 0.9],
    "dp": [0.2, 0.9],
}

# Default colormaps to use for the variables
DEFAULT_CMAP = {
    "velocity": "Plasma (matplotlib)",
    "density": "Viridis (matplotlib)",
    "pressure": "Viridis (matplotlib)",
    "temperature": "Inferno (matplotlib)",
    "b": "Cool to Warm",
    "bx": "Cool to Warm",
    "by": "Cool to Warm",
    "bz": "Cool to Warm",
    "dp": "Plasma (matplotlib)",
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


class App(pv_protocols.ParaViewWebProtocol):
    def __init__(self, dirname):
        """
        Heliosphere 4D dataset representation in Paraview.

        This class is the base app designed to enable an easy storage
        and access layer to the data, which is stored in the models
        module.
        """
        # Initialize the PV web protocols
        super().__init__()
        # Save the data directory
        self._data_dir = pathlib.Path(dirname)

        self.model = models.Enlil(self._data_dir)
        # TODO: Implement a check for whether we are given Enlil or EUHFORIA data
        # self.model = models.Euhforia(self._data_dir / "euhforia_cone_cme_example")

        # Store the data locally as well
        self.celldata = self.model.data

        # Force all cell data to point data in the volume
        self.data = pvs.CellDatatoPointData(
            registrationName=f"3D-CellDatatoPointData", Input=self.celldata
        )
        self.data.ProcessAllArrays = 1
        self.data.PassCellData = 1

        self.bvec = pvs.Calculator(registrationName="3D-Bvec", Input=self.data)
        self.bvec.AttributeType = "Point Data"
        self.bvec.ResultArrayName = "Bvec"
        self.bvec.Function = (
            f"{self.model.get_variable('bx')}*iHat"
            f" + {self.model.get_variable('by')}*jHat"
            f" + {self.model.get_variable('bz')}*kHat"
        )

        self.time_string = pvs.Text(registrationName="Time")
        # Don't add in any text right now
        self.time_string.Text = ""
        self._previous_time = None

        # create a new 'Threshold' to represent the CME
        self.threshold_cme = pvs.Threshold(registrationName="CME", Input=self.data)
        # We really only want a minimum value, so just set the maximum high
        self.threshold_cme.ThresholdRange = [1e-5, 1e5]
        # DP is the variable name in Enlil
        self.threshold_cme.Scalars = ["CELLS", self.model.get_variable("dp")]
        self.cme = pvs.Contour(registrationName="contoured_cme", Input=self.data)
        self.cme.ContourBy = ["POINTS", self.model.get_variable("dp")]
        self.cme.ComputeNormals = 0
        self.cme.Isosurfaces = [0.2]
        self.cme.PointMergeMethod = "Uniform Binning"

        self.cme_contours = pvs.Contour(
            registrationName="CME-contour", Input=self.threshold_cme
        )
        self.cme_contours.ContourBy = ["POINTS", self.model.get_variable("density")]
        self.cme_contours.Isosurfaces = []
        self.cme_contours.PointMergeMethod = "Uniform Binning"

        # Create a threshold that can be modified by the user, we use
        # two contours here instead because it looks a bit nicer.
        self.threshold = pvs.Contour(registrationName="Threshold", Input=self.data)
        self.threshold.ContourBy = ["POINTS", self.model.get_variable("density")]
        self.threshold.Isosurfaces = [10, 50]
        self.threshold.PointMergeMethod = "Uniform Binning"

        # Create the slices
        self.lon_slice = slice.Slice(
            self.celldata, slice_type="Plane", normal=(0, 0, 1), name="Longitude"
        )
        self.lat_slice = slice.Slice(
            self.celldata, slice_type="Plane", normal=(0, 1, 0), name="Latitude"
        )
        self.radial_slice = slice.Slice(
            self.celldata, slice_type="Sphere", radius=1, name="Radial"
        )

        # Dictionary mapping of string names to the object
        self.objs = {
            s: getattr(self, s)
            for s in (
                "cme",
                "data",
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
        bzLUT = pvs.GetColorTransferFunction(self.model.get_variable("bz"))
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
        bzPWF = pvs.GetOpacityTransferFunction(self.model.get_variable("bz"))
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
        disp.ColorArrayName = ["POINTS", self.model.get_variable("bz")]
        disp.LookupTable = bzLUT

        disp = pvs.Show(self.threshold, self.view, "GeometryRepresentation")
        self.displays[self.threshold] = disp
        disp.Representation = "Surface"
        disp.ColorArrayName = ["POINTS", self.model.get_variable("bz")]
        disp.LookupTable = bzLUT

        # Set colormaps
        for name in DEFAULT_CMAP:
            self.set_colormap(name)

        # hide this data from the default initial view
        for x in [
            self.threshold,
            self.cme_contours,
        ]:
            pvs.Hide(x, self.view)

        self.lon_slice.hide()
        self.radial_slice.hide()
        # We don't want any lat streamlines
        self.lat_slice.hide_streamlines()

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
        # Add the fieldlines to the satellites + Earth
        for sat in self.satellites:
            self.satellites[sat].add_fieldline(self.bvec)
        self.earth.add_fieldline(self.bvec)

    @exportRpc("pv.h3lioviz.get_available_runs")
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

    @exportRpc("pv.h3lioviz.get_variable_range")
    def get_variable_range(self, name):
        """
        Get the range of values for a variable at the current timestep.

        name : str
            Name of variable to colormap all of the surfaces by
        """
        variable = self.model.get_variable(name)
        return self.celldata.CellData.GetArray(variable).GetRange()

    @exportRpc("pv.h3lioviz.directory")
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

    @exportRpc("pv.h3lioviz.visibility")
    def change_visibility(self, obj, visibility):
        """
        Change the visibility of an object.

        obj : str
            Name of the object to update
        visibility : str ("on", "off")
            What to set the visibility to for the given object
        """
        if visibility == "on":
            if obj == "lon_slice":
                self.lon_slice.show()
            elif obj == "lat_slice":
                self.lat_slice.show()
            elif obj == "radial_slice":
                self.radial_slice.show()
            elif obj == "lon_streamlines":
                self.lon_slice.show_streamlines()
            elif obj == "lat_streamlines":
                self.lat_slice.show_streamlines()
            elif obj == "sat_fieldlines":
                for sat in self.satellites:
                    self.satellites[sat].show_fieldline()
                self.earth.show_fieldline()
            else:
                pvs.Show(self.objs[obj], self.view)

        elif visibility == "off":
            if obj == "lon_slice":
                self.lon_slice.hide()
            elif obj == "lat_slice":
                self.lat_slice.hide()
            elif obj == "radial_slice":
                self.radial_slice.hide()
            elif obj == "lon_streamlines":
                self.lon_slice.hide_streamlines()
            elif obj == "lat_streamlines":
                self.lat_slice.hide_streamlines()
            elif obj == "sat_fieldlines":
                for sat in self.satellites:
                    self.satellites[sat].hide_fieldline()
                self.earth.hide_fieldline()
            else:
                pvs.Hide(self.objs[obj], self.view)
        else:
            return ["Visibility can only be 'on' or 'off'"]

    @exportRpc("pv.h3lioviz.colorby")
    def change_color_variable(self, name):
        """
        Change the visibility of an object.

        name : str
            Name of variable to colormap all of the surfaces by
        """
        # Use a dictionary to map the variable received to the internal name
        variable = self.model.get_variable(name)
        label = VARIABLE_LABEL[name]

        # Update all displays to be colored by this variable
        for obj, disp in self.displays.items():
            if obj in (self.cme,):
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

        self.lon_slice.variable = variable
        self.lat_slice.variable = variable
        self.radial_slice.variable = variable
        self.update_opacity(name)
        self.update_lut(name)

        # hides old scalarbars that aren't in the view and
        # shows the new variable we are using now
        pvs.UpdateScalarBars(self.view)
        # But we want to hide the streamlines colorbar
        disp = self.lon_slice.streamlines_disp.SetScalarBarVisibility(self.view, False)
        disp = self.lat_slice.streamlines_disp.SetScalarBarVisibility(self.view, False)

        # restore active source
        pvs.SetActiveSource(None)
        # Render the view
        pvs.Render(self.view)

    def update_lut(self, name):
        """
        Set the variable range of the lookup table.

        name : str
            Name of variable to update
        """
        variable = self.model.get_variable(name)
        lut = pvs.GetColorTransferFunction(variable)
        lut.RescaleTransferFunction(LUT_RANGE[name])
        lut.AutomaticRescaleRangeMode = "Never"

    def update_opacity(self, name):
        """
        Set the variable range of the opacity lookup table.

        name : str
            Name of variable to update
        """
        variable = self.model.get_variable(name)
        opacity_map = pvs.GetOpacityTransferFunction(variable)
        # Create the control points
        points = OPACITY_VALUES[name]
        data_range = LUT_RANGE[name]
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

    @exportRpc("pv.h3lioviz.set_colormap")
    def set_colormap(self, name, cmap_name=None):
        """
        Set the colormap for the variable.

        name : str
            Name of the variable to set the colormap of
        cmap_name : str
            Name of the colormap to apply
        """
        # Use a dictionary to map the variable received to the internal name
        variable = self.model.get_variable(name)
        lut = pvs.GetColorTransferFunction(variable)
        # If cmap_name is None, use the default version
        lut.ApplyPreset(cmap_name or DEFAULT_CMAP[name])
        lut.EnableOpacityMapping = 1

    @exportRpc("pv.h3lioviz.set_range")
    def set_range(self, name, range):
        """
        Set the range of values used for colormapping.

        name : str
            Name of the variable to set the colormap of
        range : list[2]
            A list of the minimum and maximum values to colormap over
        """
        LUT_RANGE[name] = range
        self.update_lut(name)

    @exportRpc("pv.h3lioviz.set_opacity")
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
        OPACITY_VALUES[name] = range
        self.update_opacity(name)

    @exportRpc("pv.h3lioviz.set_threshold")
    def set_threshold(self, name, range):
        """
        Set the variable and range of values to be used for the threshold.

        name : str
            Name of the variable to use for the thresholding
        range : list[2]
            A list of the minimum and maximum values to threshold by
        """
        variable = self.model.get_variable(name)
        # The quantity of interest
        self.threshold.ContourBy = ["POINTS", variable]
        self.threshold.Isosurfaces = range

    @exportRpc("pv.h3lioviz.set_contours")
    def set_contours(self, name, values):
        """
        Set the variable and a list of values to be used for the contours.

        name : str
            Name of the variable to use for the contouring
        values : list
            A list of the values to contour by
        """
        variable = self.model.get_variable(name)
        # The quantity of interest
        self.cme_contours.ContourBy = ["POINTS", variable]
        self.cme_contours.Isosurfaces = values

    @exportRpc("pv.h3lioviz.snap_solar_plane")
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

        self.lon_slice.slice_data.SliceType.Normal = loc
        # Also update the stream source so they stay in-sync
        self.lon_slice.stream_source.Normal = loc

    @exportRpc("pv.h3lioviz.rotate_plane")
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
            self.lon_slice.slice_data.SliceType.Normal = loc
            self.lon_slice.stream_source.Normal = loc
        elif plane == "lat":
            loc = [-y, x, 0]
            self.lat_slice.slice_data.SliceType.Normal = loc
            self.lat_slice.stream_source.Normal = loc
        else:
            raise ValueError("You can only update the 'lon' or 'lat' plane.")

    @exportRpc("pv.h3lioviz.snap_to_view")
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

    @exportRpc("pv.h3lioviz.toggle_satellites")
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

    @exportRpc("pv.h3lioviz.get_satellite_times")
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

    @exportRpc("pv.h3lioviz.get_satellite_data")
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
