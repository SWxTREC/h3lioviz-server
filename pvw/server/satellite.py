import pathlib
import paraview.simple as pvs
from evolution import Evolution

# List of satellite colors
SATELLITE_COLORS = {
    "earth": [0.0, 0.3333333333333333, 0.0],
    "stereoa": [177 / 255, 138 / 255, 142 / 255],
    "stereob": [94 / 255, 96 / 255, 185 / 255],
}

# Keep track of the name mapping that we want to show to users
SATELLITE_NAMES = {"earth": "Earth", "stereoa": "STEREO-A", "stereob": "STEREO-B"}


class Satellite:
    """
    A satellite or observation point in space.

    Parameters
    ----------
    name : str
        Name of the satellite
    fname : Path
        Path of the evolution data file for this satellite
    representation : str
        What type of object to display the satellite as, one of
        ("box", "sphere")
    size : Number
        Edge-length of the box, or radius of the sphere
    center : 3-tuple
        (x, y, z) position of the satellite
    view : paraview.Simple.View
        A view to render the satellite into
    """

    def __init__(self, name, fname, representation="box", size=0.02, view=None):
        self.name = name
        self.size = size
        self.view = view
        if representation == "box":
            self.sat = pvs.Box()
            self.sat.XLength = size
            self.sat.YLength = size
            self.sat.ZLength = size
        elif representation == "sphere":
            self.sat = pvs.Sphere()
            self.sat.Radius = size
        else:
            raise ValueError("Satellite representation can only be 'box' or 'sphere'")

        self.evolution = Evolution(fname)
        self.color = (
            SATELLITE_COLORS[name] if name in SATELLITE_COLORS else [0.5, 0.5, 0.5]
        )
        self.sat_disp = pvs.Show(self.sat, self.view, "GeometryRepresentation")
        self.sat_disp.Representation = "Surface"
        self.sat_disp.AmbientColor = self.color
        self.sat_disp.ColorArrayName = [None, ""]
        self.sat_disp.DiffuseColor = self.color

        self.label = pvs.Text()
        self.label.Text = SATELLITE_NAMES[name] if name in SATELLITE_NAMES else name
        self.label_disp = pvs.Show(self.label, self.view, "TextSourceRepresentation")
        self.label_disp.TextPropMode = "Billboard 3D Text"
        self.label_disp.FontSize = 14
        self.label_disp.WindowLocation = "AnyLocation"
        # Offset the center by the width to give some separation
        self.label_disp.BillboardPosition = [x + size for x in self.sat.Center]
        self.label_disp.Color = [0, 0, 0]  # Black text

    def add_fieldline(self, data):
        """
        Add a fieldline trace to this satellite.

        The seed point is at the satellite location and the tube
        uses the same color as the satellite.

        Parameters
        ----------
        data : Paraview point dataset
            3D Point dataset used for fieldline integration
        """
        stream_input = pvs.StreamTracer(
            registrationName=f"{self.name}-StreamTracer",
            Input=data,
            SeedType="Point Cloud",
        )
        stream_input.Vectors = ["POINTS", "Bvec"]
        # One point, at the satellite's center
        stream_input.SeedType.NumberOfPoints = 1
        stream_input.SeedType.Center = self.sat.Center
        stream_input.SeedType.Radius = 0.0
        stream_input.MaximumStreamlineLength = 3.4
        stream_input.ComputeVorticity = 0

        streamlines = pvs.Tube(
            registrationName=f"{self.name}-Streamlines", Input=stream_input
        )
        streamlines.Capping = 1
        streamlines.Radius = 0.0025

        # Store the objects we need for later as attributes on self
        self.stream_input = stream_input
        self.streamlines = streamlines

        # Streamlines
        self.streamlines_disp = pvs.Show(
            self.streamlines, self.view, "GeometryRepresentation"
        )
        self.streamlines_disp.Representation = "Surface"
        self.streamlines_disp.ColorArrayName = [None, ""]
        self.streamlines_disp.AmbientColor = self.color
        self.streamlines_disp.DiffuseColor = self.color

    def show(self):
        """Make the satellite visible"""
        pvs.Show(self.sat, self.view)
        pvs.Show(self.label, self.view)

    def hide(self):
        """Hide the satellite"""
        pvs.Hide(self.sat, self.view)
        pvs.Hide(self.label, self.view)

    def show_fieldline(self):
        pvs.Show(self.streamlines, self.view)

    def hide_fieldline(self):
        pvs.Hide(self.streamlines, self.view)

    def update(self, time):
        """
        Update the satellite position to this time

        Parameters
        ----------
        time : datetime
            The time of the visualization
        """
        self.sat.Center = self.evolution.get_position(time=time)
        self.label_disp.BillboardPosition = [x + self.size for x in self.sat.Center]
        if hasattr(self, "stream_input"):
            self.stream_input.SeedType.Center = self.sat.Center


class Sun:
    """
    Representation of the sun

    Parameters
    ----------
    texture_dir : Path
        Path of the textures for the sun
    size : Number
        radius of the sun
    view : paraview.Simple.View
        A view to render the sun into
    """

    def __init__(self, texture_dir=None, size=0.075, view=None):
        self.texture_dir = texture_dir
        self.size = size
        self.view = view

        # Sun representation
        self.sphere = pvs.Sphere()
        self.sphere.Center = [0.0, 0.0, 0.0]
        self.sphere.Radius = self.size - 0.001
        self.sphere.ThetaResolution = 50
        self.sphere.PhiResolution = 50
        disp = pvs.Show(self.sphere, self.view, "GeometryRepresentation")

        # Yellow color
        self.color = [0.8313725490196079, 0.8313725490196079, 0.0]
        disp.Representation = "Surface"
        disp.AmbientColor = self.color
        disp.ColorArrayName = [None, ""]
        disp.DiffuseColor = self.color

        # Now try to create the solar texture
        self._create_solar_texture()

    def _create_solar_texture(self):
        if self.texture_dir is None or not self.texture_dir.exists():
            return
        # Store a list of the solar images
        self._solar_images = [x.name for x in self.texture_dir.glob("*.jpg")]
        self._solar_images = sorted(self._solar_images)

        # Sun representation
        sun = pvs.Sphere()
        sun.Center = self.sphere.Center
        sun.Radius = self.size
        sun.ThetaResolution = 50
        sun.PhiResolution = 50

        # For the solar imagery we want texture map to plane because it is
        # a flat image instead of an unwrapped image.
        texture_map = pvs.TextureMaptoPlane(registrationName="SunImage", Input=sun)
        # Make the points form a square of radius r
        # Earth is in the -X direction, so we want our image plane to be
        # in the Y-Z direction, with the origin at (+Y, -Z) and the base
        # of the image extending out in the Y direction to (-Y, -Z).
        r = self.size
        texture_map.Origin = [0, r, -r]
        texture_map.Point1 = [0, -r, -r]
        texture_map.Point2 = [0, r, r]

        # We also want to clip the sphere so we don't get any wrapping
        # into the back plane
        clip = pvs.Clip(registrationName="ClipSun", Input=texture_map)
        clip.ClipType = "Plane"
        clip.HyperTreeGridClipper = "Plane"
        clip.Scalars = ["POINTS", ""]
        clip.Invert = 1

        # This is the plane to clip on. It doesn't cover the entire
        # half-sphere, so limit it a little bit in the X direction
        clip.ClipType.Origin = [-0.03, 0.0, 0.0]

        sun_display = pvs.Show(clip, self.view, "GeometryRepresentation")

        # Create a texture from the first image
        sun_texture = pvs.CreateTexture(str(self.texture_dir / self._solar_images[0]))

        # trace defaults for the display properties.
        sun_display.Representation = "Surface"
        sun_display.ColorArrayName = [None, ""]
        sun_display.SelectTCoordArray = "Texture Coordinates"
        sun_display.SelectNormalArray = "Normals"
        sun_display.SelectTangentArray = "None"
        sun_display.Texture = sun_texture
        # This hides the HMI image when looking from behind
        sun_display.BackfaceRepresentation = "Cull Backface"
        self.texture_disp = sun_display

    def update(self, time):
        """
        Update the sun image

        Parameters
        ----------
        time : datetime
            The time of the visualization
        """
        if not hasattr(self, "_solar_images"):
            # No solar images
            return

        # Iterate through the solar images, choosing
        # the one before this timestep.
        # filename looks like: 20170906_000000_M_color_4k.jpg
        i = len(self._solar_images) - 1
        image_name = self._solar_images[i]
        t = time.strftime("%Y%m%d_%H0000_M_color_4k.jpg")
        while image_name > t and i > 0:
            image_name = self._solar_images[i]
            i -= 1

        # We have our image_name now, so update the texture
        self.texture_disp.Texture = pvs.CreateTexture(
            str(self.texture_dir / image_name)
        )


class Earth(Satellite):
    """
    Representation of the Earth

    Parameters
    ----------
    fname : Path
        Path of the evolution data file for the Earth
    size : Number
        radius of the Earth
    view : paraview.Simple.View
        A view to render the sun into
    """

    def __init__(self, fname, size=0.025, view=None):
        # inherit from Satellites for evolution data file parsing
        # and setting up the initial sphere
        super().__init__("earth", fname, representation="sphere", size=size, view=view)

        # Path to the Earth texture on our local system
        # cwd() is where paraview is launched from
        earth_path = pathlib.Path.cwd() / "pvw" / "server" / "assets"
        earth_path /= "land_shallow_topo_2048.jpg"
        # If we don't have the texture file, go download it.
        if not earth_path.exists():
            # Make the directories if they don't already exist
            earth_path.parent.mkdir(parents=True, exist_ok=True)
            import urllib.request

            url = (
                "https://eoimages.gsfc.nasa.gov/images/imagerecords/"
                "57000/57752/land_shallow_topo_2048.jpg"
            )
            # Make a request
            req = urllib.request.urlopen(url)
            # Write out the response to our local file
            with open(earth_path, "wb") as f:
                f.write(req.read())

        # We should have a local image file to use as a texture
        earth_image = pvs.CreateTexture(str(earth_path))
        # Set up the sphere source
        sphere = self.sat
        sphere.ThetaResolution = 50
        # We need to perturb the StartTheta a small amount to not have a
        # seam/mismatch in the texture at 0
        sphere.StartTheta = 1e-3
        sphere.PhiResolution = 50
        # create a new 'Texture Map to Sphere'
        texture_map = pvs.TextureMaptoSphere(
            registrationName="EarthImage", Input=sphere
        )
        texture_map.PreventSeam = 0

        # Move the Earth sphere center back to zero for a translation
        # after rotation later.
        self.sat.Center = [0, 0, 0]

        # We want to rotate the Earth image with the hour of the day
        # To do that, we need to translate our object to the opposite
        # location of what it truly is, then rotate about the origin's
        # z-axis, then translate to the final location after the rotation.
        # TODO: This could be turned into a rotation matrix eventually to
        #       only calculate this at one step, rather than three successive
        #       filters being applied.
        t = pvs.Transform(registrationName="EarthTranslation1", Input=texture_map)
        t.Transform = "Transform"
        t.Transform.Translate = [0.0, 0.0, 0.0]
        self.translation1 = t
        t = pvs.Transform(registrationName="EarthRotation", Input=self.translation1)
        t.Transform = "Transform"
        t.Transform.Rotate = [0.0, 0.0, 0.0]
        self.rotation = t
        t = pvs.Transform(registrationName="EarthTranslation2", Input=self.rotation)
        t.Transform = "Transform"
        t.Transform.Translate = [0.0, 0.0, -1]
        self.translation2 = t

        # show data from the image and hide the plain sphere and label
        pvs.Hide(self.sat)
        pvs.Hide(self.label)
        texture_map_disp = pvs.Show(
            self.translation2, self.view, "GeometryRepresentation"
        )

        # trace defaults for the display properties.
        texture_map_disp.Representation = "Surface"
        texture_map_disp.ColorArrayName = [None, ""]
        texture_map_disp.SelectTCoordArray = "Texture Coordinates"
        texture_map_disp.SelectNormalArray = "Normals"
        texture_map_disp.SelectTangentArray = "None"
        texture_map_disp.Texture = earth_image
        # To get the proper orientation
        texture_map_disp.FlipTextures = 1

    def show(self):
        pvs.Show(self.translation2, self.view)

    def hide(self):
        pvs.Hide(self.translation2, self.view)

    def update(self, time):
        """
        Update the Earth image to this time

        Parameters
        ----------
        time : datetime
            The time of the visualization
        """
        # Update the position of the Earth first
        super().update(time)
        if not hasattr(self, "rotation"):
            # There is no earth image to rotate
            return
        # We want rotation to be from 0 -> 360
        rot = (time.hour + time.minute / 60) / 24 * 360
        # Rotate around Z with the hours of the day
        earth_pos = self.sat.Center
        # Move it negative first to apply the rotation
        self.translation1.Transform.Translate = [-x for x in earth_pos]
        self.rotation.Transform.Rotate = [0.0, 0.0, rot]
        # Then move it back positive to its actual location
        self.translation2.Transform.Translate = earth_pos
