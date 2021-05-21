import argparse

# import paraview modules.
from paraview.web import pv_wslink
from paraview.web import protocols as pv_protocols

from paraview import simple
from wslink import server

from enlil import EnlilDataset
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
        self.enlil = EnlilDataset(fname)
        # Register the Paraview protocols for dispatching methods
        self.registerVtkWebProtocol(self.enlil)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="ParaViewWeb Demo")

    # Add default arguments
    server.add_arguments(parser)

    # Extract arguments
    args = parser.parse_args()

    # Start server
    server.start_webserver(options=args, protocol=_DemoServer)
