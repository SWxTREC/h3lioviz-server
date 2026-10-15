import argparse
import asyncio
import os
import signal
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
    # Started as the launcher's spare process (pvw/launcher/spare_launcher.py)
    spare = False
    # The idle timeout once a spare process belongs to a client
    client_timeout = 0
    # The wslink web server, once wslink has created it
    web_server = None
    # Connected clients
    clients = 0

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
        parser.add_argument(
            "--spare",
            action="store_true",
            help=(
                "Started by the launcher before any client asked for it: no idle "
                "timeout until the launcher hands it to a client (SIGUSR1), then "
                "--timeout applies"
            ),
            dest="spare",
        )

    @staticmethod
    def configure(args):
        # Update this server based on the passed in arguments
        _AppServer.authKey = args.authKey
        _AppServer.data_dir = args.data_dir
        _AppServer.settingsLODThreshold = args.settingsLODThreshold
        _AppServer.host = args.host
        _AppServer.spare = args.spare

    def set_server(self, web_server):
        """Called by wslink with its web server, before it starts serving"""
        _AppServer.web_server = web_server

    def onConnect(self, request, client_id):
        _AppServer.clients += 1

    def onClose(self, client_id):
        # wslink also calls this when a connection failed before onConnect
        _AppServer.clients = max(0, _AppServer.clients - 1)

    @staticmethod
    def handed_to_client():
        """
        The launcher has handed this spare process to a client: from now on
        it times out when idle like a process started for that client,
        including when the client never connects.
        """
        if not _AppServer.spare:
            return
        _AppServer.spare = False
        # wslink 1.12 also keeps its web server in server.ws_server
        web_server = _AppServer.web_server or getattr(server, "ws_server", None)
        if web_server is None:
            # Not reached with wslink 1.12, which sets both. Keep serving the
            # client rather than drop it, but say the timeout is not applied
            print("Handed to a client, but no web server to set the idle timeout on")
            return
        web_server.config["timeout"] = _AppServer.client_timeout
        # wslink registers a connection before it calls onConnect, so ask its
        # handlers (WebAppServer._ws_handlers in wslink 1.12) whether one is
        # open; the onConnect count is the fallback
        handlers = getattr(web_server, "_ws_handlers", None)
        if handlers is not None:
            connected = any(handler.connections for handler in handlers)
        else:
            connected = _AppServer.clients > 0
        if not connected:
            web_server.shutdown_schedule()

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
    if args.spare:
        try:
            from wslink.protocol import AbstractWebApp

            spare_supported = hasattr(AbstractWebApp, "shutdown_schedule")
        except ImportError:
            spare_supported = False
        if spare_supported:
            # No idle timeout until the launcher hands this process to a
            # client, which it signals with SIGUSR1
            _AppServer.client_timeout = args.timeout
            args.timeout = 0
            signal.signal(
                signal.SIGUSR1,
                lambda signum, frame: asyncio.get_event_loop().call_soon_threadsafe(
                    _AppServer.handed_to_client
                ),
            )
        else:
            # This wslink cannot change the timeout later: run as an ordinary
            # session, and ignore the launcher's signal rather than die of it
            print("--spare needs a newer wslink (1.12.4 works); running as an ordinary session")
            args.spare = False
            signal.signal(signal.SIGUSR1, signal.SIG_IGN)
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
