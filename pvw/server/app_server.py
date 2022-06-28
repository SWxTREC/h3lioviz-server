import argparse

# import paraview modules.
from paraview.web import pv_wslink
from paraview.web import protocols as pv_protocols

from paraview import simple
from wslink import server

from app import App

# =============================================================================
# Create custom PVServerProtocol class to handle clients requests
# =============================================================================


class _AppServer(pv_wslink.PVServerProtocol):
    authKey = "wslink-secret"
    viewportScale = 1.0
    viewportMaxWidth = 2560
    viewportMaxHeight = 1440
    settingsLODThreshold = 102400

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--dir",
            default="/data",
            help=("Path to the data directory to load"),
            dest="data_dir",
        )
        parser.add_argument(
            "--viewport-scale",
            default=1.0,
            type=float,
            help="Viewport scaling factor",
            dest="viewportScale",
        )
        parser.add_argument(
            "--viewport-max-width",
            default=2560,
            type=int,
            help="Viewport maximum size in width",
            dest="viewportMaxWidth",
        )
        parser.add_argument(
            "--viewport-max-height",
            default=1440,
            type=int,
            help="Viewport maximum size in height",
            dest="viewportMaxHeight",
        )
        parser.add_argument(
            "--settings-lod-threshold",
            default=102400,
            type=int,
            help="LOD Threshold in Megabytes",
            dest="settingsLODThreshold",
        )

    @staticmethod
    def configure(args):
        # Update this server based on the passed in arguments
        _AppServer.authKey = args.authKey
        _AppServer.data_dir = args.data_dir
        _AppServer.viewportScale = args.viewportScale
        _AppServer.viewportMaxWidth = args.viewportMaxWidth
        _AppServer.viewportMaxHeight = args.viewportMaxHeight
        _AppServer.settingsLODThreshold = args.settingsLODThreshold

    def initialize(self):
        # Bring used components
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebMouseHandler())
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebViewPort())
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebTimeHandler())
        self.registerVtkWebProtocol(
            pv_protocols.ParaViewWebPublishImageDelivery(decode=False)
        )
        self.updateSecret(_AppServer.authKey)

        # tell the C++ web app to use no encoding.
        # ParaViewWebPublishImageDelivery must be set to decode=False to match.
        self.getApplication().SetImageEncoding(0)

        # Disable interactor-based render calls
        simple.GetRenderView().EnableRenderOnInteraction = 0

        # The directory containing the NetCDF file with the data
        self.app = App(self.data_dir)
        # Register the Paraview protocols for dispatching methods
        self.registerVtkWebProtocol(self.app)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="H3LIO Viewer")

    # Add default arguments
    server.add_arguments(parser)
    _AppServer.add_arguments(parser)

    # Extract arguments
    args = parser.parse_args()
    _AppServer.configure(args)

    # Start server
    server.start_webserver(options=args, protocol=_AppServer)
