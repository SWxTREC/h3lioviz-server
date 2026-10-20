"""
wslink's launcher, with one session process started ahead of time.

Starting a session process (pvpython, the imports and the App) takes seconds,
and with the stock launcher the user waits for it after every launch request
(POST). This launcher keeps one "spare" session process running for one
application and hands it to the next request for that application, then
starts another spare. Otherwise it is wslink's own launcher, with the same
configuration file, the same reply to the POST, and the same proxy-mapping
file that Apache reads to route each session's websocket, plus three
changes that apply with or without a spare: each session process runs in
its own process group and is signalled as a group, the launcher ignores
SIGUSR1, and a port that something still listens on is not handed out.

The launcher configuration gets one more top-level entry:

    "spare": {"application": "visualizer", "delay": 10}

"application" is the app to keep a spare of; "delay" is how many seconds
after handing one out to start the next, so that it does not compete with
the new session's first load. Without the entry no spare is kept.

The spare is started with "--spare" added to its command line, and gets
SIGUSR1 when it is handed out; see pvw/server/app_server.py.

Usage, as for wslink.launcher:

    pvpython spare_launcher.py config.json
"""

import asyncio
import io
import logging
import os
import signal
import socket
import string
import subprocess

from wslink import launcher
from wslink.backends.aiohttp import launcher as aiohttp_launcher

logger = logging.getLogger("wslink.spare_launcher")

SPARE_ARGUMENT = "--spare"

# Session fields the launcher sets itself, whatever the request says. (The
# stock launcher also overwrites the app's sessionData keys, such as updir,
# after the request's; on the spare path the request's win. Neither is in the
# command line or the reply.)
LAUNCHER_FIELDS = {"id", "host", "port", "sessionURL", "cmd"}


def template_fields(templates):
    """The ${name} and $name fields used in a list of template strings"""
    fields = set()
    for template in templates:
        for match in string.Template.pattern.finditer(template):
            name = match.group("named") or match.group("braced")
            if name:
                fields.add(name)
    return fields


def port_is_free(port):
    """Whether a session process could listen on this port (app_server.py
    listens on 0.0.0.0)"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        # As the session's own listening socket does (asyncio sets
        # SO_REUSEADDR on Unix), so a port held only by a closed session's
        # TIME_WAIT connections counts as free
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind(("0.0.0.0", port))
            return True
        except OSError:
            return False


class FreePortResourceManager(launcher.ResourceManager):
    """
    wslink's port allocation, skipping ports that something still listens
    on: a session process the launcher has lost track of would make every
    session given that port fail to start
    """

    def getNextResource(self):
        busy = []
        host, port = super().getNextResource()
        while host is not None and not port_is_free(port):
            logger.warning("Port %s is still in use; not handing it out", port)
            busy.append((host, port))
            host, port = super().getNextResource()
        for busy_host, busy_port in busy:
            # Free again, but at the end of the queue, to be tried last
            self.resources[busy_host]["used"].remove(busy_port)
            self.resources[busy_host]["available"].insert(0, busy_port)
        return host, port


class SessionProcessManager(launcher.ProcessManager):
    """
    wslink's process manager, with each session started in a process group
    of its own and signalled as a group.

    On Linux, ParaView's bin/pvpython is a small launcher that runs
    bin/pvpython-real as a child process and waits for it. A signal sent to
    that launcher alone never reaches pvpython-real, and kills the launcher if
    it does not handle the signal. The session then counts as ended and its
    port is freed while pvpython-real still serves on it.
    """

    def startProcess(self, session):
        # As wslink's, plus start_new_session
        log_path = self._getLogFilePath(session["id"])
        with io.open(log_path, mode="a+", buffering=1, encoding="utf-8") as log_file:
            try:
                proc = subprocess.Popen(
                    session["cmd"],
                    stdout=log_file,
                    stderr=log_file,
                    start_new_session=True,
                )
            except OSError:
                logger.error("The command line failed")
                logger.error(" ".join(map(str, session["cmd"])))
                return None
        self.processes[session["id"]] = proc
        return proc

    def signal(self, session_id, signum):
        """Send a signal to every process of a session; False if none is left"""
        try:
            os.killpg(self.processes[session_id].pid, signum)
            return True
        except (KeyError, ProcessLookupError):
            return False

    def stopProcess(self, id):
        # The whole group: a pvpython-real left running keeps the port
        self.signal(id, signal.SIGTERM)
        del self.processes[id]

    def __del__(self):
        for session_id in list(self.processes):
            self.signal(session_id, signal.SIGTERM)


class SpareLauncherResource(aiohttp_launcher.LauncherResource):
    def __init__(self, options, config):
        super().__init__(options, config)
        self.process_manager = SessionProcessManager(config)
        self.session_manager.resources = FreePortResourceManager(config["resources"])
        spare_config = config.get("spare", {})
        self.spare_application = spare_config.get("application")
        self.spare_delay = float(spare_config.get("delay", 10))
        self.spare = None
        self.spare_pending = False
        if self.spare_application:
            app = config["apps"][self.spare_application]
            # Request fields that would change the spare's command line or
            # session URL: a request setting any of them gets its own process
            self.spare_fields = template_fields(
                app["cmd"]
                + [config["configuration"]["sessionURL"]]
                + list(config.get("sessionData", {}).values())
            ) - {"id", "host", "port"}
            self.start_spare()

    def start_spare(self):
        """Start a spare session process, if there is none"""
        self.spare_pending = False
        if self.spare is not None:
            if self.session_running(self.spare):
                return
            logger.warning("Spare session %s has exited", self.spare["id"])
            self.release(self.spare)
            self.spare = None
        session = self.session_manager.createSession(
            {"application": self.spare_application}
        )
        if session is None:
            logger.warning("No free port for a spare session process")
            return
        session["cmd"] = session["cmd"] + [SPARE_ARGUMENT]
        if not self.process_manager.startProcess(session):
            self.session_manager.deleteSession(session["id"])
            return
        logger.info("Started spare session %s on port %s", session["id"], session["port"])
        self.spare = session

    def schedule_spare(self):
        """Start the next spare after the configured delay"""
        if self.spare_application and self.spare is None and not self.spare_pending:
            self.spare_pending = True
            asyncio.get_running_loop().call_later(self.spare_delay, self.start_spare)

    def session_running(self, session):
        """Whether a session's process is still known and running"""
        session_id = session["id"]
        return (
            session_id in self.session_manager.sessions
            and session_id in self.process_manager.processes
            and self.process_manager.isRunning(session_id)
        )

    def release(self, session):
        """Free a session's port and process, whatever state they are in"""
        if session["id"] in self.session_manager.sessions:
            self.session_manager.deleteSession(session["id"])
        if session["id"] in self.process_manager.processes:
            self.process_manager.stopProcess(session["id"])

    def take_spare(self, payload):
        """Remove and return the spare session, if it can serve this request"""
        if (
            self.spare is None
            or not isinstance(payload, dict)
            or payload.get("application") != self.spare_application
            or self.spare_fields & set(payload)
        ):
            return None
        spare, self.spare = self.spare, None
        if not self.session_running(spare):
            logger.warning("Spare session %s has exited", spare["id"])
            self.release(spare)
            return None
        return spare

    def release_ended(self):
        """Free the ports and mapping lines of sessions whose process has
        ended, as the stock launcher does on every request"""
        for session_id in self.process_manager.listEndedProcess():
            self.session_manager.deleteSession(session_id)
            self.process_manager.stopProcess(session_id)

    async def handle_post(self, request):
        try:
            payload = await request.json()
        except ValueError:
            payload = None
        self.release_ended()
        spare = self.take_spare(payload)
        if spare is None:
            # The stock launcher starts a process for this request
            response = await super().handle_post(request)
            self.schedule_spare()
            return response

        # Make the session what the stock launcher would have made for this
        # request; the fields that matter to the command line are excluded
        # by take_spare
        for key, value in payload.items():
            if key not in LAUNCHER_FIELDS:
                spare[key] = value
        response = await self._waitForReady(spare, request)
        if response.status == launcher.STATUS_OK:
            # The process belongs to this client now: start its idle
            # timeout, as for a process started for this request
            if not self.process_manager.signal(spare["id"], signal.SIGUSR1):
                # It exited after all: start one for this request instead
                self.release(spare)
                response = await super().handle_post(request)
        self.schedule_spare()
        return response


def start():
    # Session processes inherit this. ParaView's pvpython launcher process
    # (see SessionProcessManager) then ignores the SIGUSR1 meant for the
    # app_server.py it runs, which sets its own handler for it (--spare).
    signal.signal(signal.SIGUSR1, signal.SIG_IGN)
    # wslink's launcher builds its LauncherResource from this module global
    aiohttp_launcher.LauncherResource = SpareLauncherResource
    launcher.start(description="wslink launcher with a spare session process")


if __name__ == "__main__":
    start()
