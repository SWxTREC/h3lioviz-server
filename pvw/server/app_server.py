import argparse
import os
import sys
import time

# import paraview modules.
from paraview.web import pv_wslink
from paraview.web import protocols as pv_protocols

from wslink import schedule_callback, server

from app import App


class _ImageDelivery(pv_protocols.ParaViewWebPublishImageDelivery):
    """
    ParaView's image delivery, with its retries for stale images fixed.

    When a render's image is not encoded yet, the reply is marked stale and
    the image is asked for again later. ParaView (5.10.1, 6.1.1 and 6.2.0)
    means to retry up to staleCountLimit times, retry n once (n + 1) times
    deltaStaleTimeBeforeRender have passed since the last stale reply. But
    each retry's own stale reply resets that time, so only the first retry
    ever happens, and an image that takes longer than one delay to encode is
    not sent until something renders again. Here each retry is rescheduled
    until its time has come.
    """

    def renderStaleImage(self, vId, staleCount=0):
        if self.staleHandlerCount.get(vId, 0) <= 0:
            return
        self.staleHandlerCount[vId] -= 1
        if self.lastStaleTime[vId] == 0:
            # An image that was not stale has been sent since
            return
        wait = self.deltaStaleTimeBeforeRender * (staleCount + 1) - (
            time.time() - self.lastStaleTime[vId]
        )
        if wait > 0:
            self.staleHandlerCount[vId] += 1
            schedule_callback(wait + 0.001, lambda: self.renderStaleImage(vId, staleCount))
        elif staleCount < self.staleCountLimit:
            self.pushRender(vId, False, staleCount + 1)


# =============================================================================
# Create custom PVServerProtocol class to handle clients requests
# =============================================================================


class _AppServer(pv_wslink.PVServerProtocol):
    authKey = "wslink-secret"
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
        _AppServer.settingsLODThreshold = args.settingsLODThreshold
        _AppServer.host = args.host

    def initialize(self):
        # Bring used components
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebMouseHandler())
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebViewPort())
        self.registerVtkWebProtocol(pv_protocols.ParaViewWebTimeHandler())
        image_delivery = _ImageDelivery(decode=False)
        # Ask for a stale image again after 0.02 s rather than ParaView's
        # 0.1 s, so it reaches the client sooner. With the retries fixed, an
        # image that takes longer to encode is still sent: the retries wait
        # 0.02, 0.04, ... 0.2 s, 1.1 s in all.
        image_delivery.deltaStaleTimeBeforeRender = 0.02
        self.registerVtkWebProtocol(image_delivery)
        self.updateSecret(_AppServer.authKey)

        # tell the C++ web app to use no encoding.
        # ParaViewWebPublishImageDelivery must be set to decode=False to match.
        self.getApplication().SetImageEncoding(0)

        # The directory containing the NetCDF file with the data
        self.app = App(self.data_dir)
        # Register the Paraview protocols for dispatching methods
        self.registerVtkWebProtocol(self.app)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="H3lioViz")

    # Add default arguments
    server.add_arguments(parser)
    _AppServer.add_arguments(parser)

    # Extract arguments
    args = parser.parse_args()
    _AppServer.configure(args)

    # Start server
    args.host = "0.0.0.0"  # aiohttp doesn't like "localhost" here
    server.start_webserver(options=args, protocol=_AppServer)

    # The session is over. Exit without running ParaView's C++ teardown:
    # under ParaView 6.x it can segfault (in vtkXRenderWindowInteractor's
    # Finalize, though this process is headless) or hang, and a hung
    # process keeps its memory until the instance stops.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)
