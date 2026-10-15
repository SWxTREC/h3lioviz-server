import paraview.simple as pvs
from paraview import servermanager


def color_by(display, variable):
    """
    Color a display by a point array as pvs.ColorBy does, but keep the color
    map's range. ColorBy also rescales the color map to this display's data,
    which the app always overrides with its own range right after
    (App.update_lut); a display created later must not undo that.
    """
    display.SetScalarColoring(variable, servermanager.GetAssociationFromString("POINTS"))


class LazyDisplay:
    """
    A display (representation) that is only created the first time it is shown.

    Under ParaView 6.2 creating a geometry representation costs about 0.2 s and
    75 MB whatever the data (measured in the 6.2 image), most of it for the
    data axes grid each one builds.
    Most of what h3lioviz builds is hidden until the user asks for it, so only
    the representations that are shown get built.

    Parameters
    ----------
    source : paraview proxy
        The source or filter to display
    view : paraview view
        The view to show it in
    representation : str
        The representation type passed to Show
    setup : callable, optional
        Called with the new display once, right after it is created, to set
        its properties
    """

    def __init__(self, source, view, representation="GeometryRepresentation", setup=None):
        self.source = source
        self.view = view
        self.representation = representation
        self.setup = setup
        self.display = None

    @property
    def created(self):
        """Whether the display exists yet"""
        return self.display is not None

    def show(self):
        """Show the source, creating its display the first time"""
        if self.display is None:
            self.display = pvs.Show(self.source, self.view, self.representation)
            if self.setup is not None:
                self.setup(self.display)
        else:
            pvs.Show(self.source, self.view)

    def hide(self):
        """Hide the source, if it was ever shown"""
        if self.display is not None:
            pvs.Hide(self.source, self.view)
