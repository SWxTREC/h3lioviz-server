import logging

from flask import Flask, abort, jsonify

import get_timeseries as timeseries
import h3lioviz_list_runs as runs
import sync_metadata as metadata

app = Flask(__name__)

logger = logging.getLogger(__name__)


@app.route("/metadata/health", methods=["GET"])
def health():
    return jsonify(status="ok"), 200


@app.route("/metadata/getTimeSeries/<run_id>/<satellite>", methods=["GET"])
def get_timeseries(run_id: str, satellite: str):
    try:
        payload = timeseries.get_timeseries(run_id, satellite)
    except FileNotFoundError:
        abort(404, description="Time series not found")
    except RuntimeError as exc:
        abort(500, description=str(exc))
    except Exception as e:
        logger.exception("Failed to fetch time series for run %s", run_id)
        abort(500, description=f"Failed to fetch time series: {str(e)}")
    return jsonify(payload), 200


@app.route("/metadata/availableRuns", methods=["GET"])
def available_runs():
    try:
        data = runs.list_runs()
    except RuntimeError as exc:
        abort(500, description=str(exc))
    except Exception as e:
        logger.exception("Failed to list available runs")
        abort(500, description=f"Failed to list runs: {e}")
    return jsonify(data), 200


@app.route("/metadata/syncMetadata", methods=["GET"])
def sync_metadata():
    try:
        report = metadata.sync_metadata()
    except RuntimeError as exc:
        abort(500, description=str(exc))
    except Exception:
        logger.exception("Failed to sync metadata")
        abort(500, description="Failed to sync metadata")
    return jsonify(report), 200

if __name__ == '__main__':
    # Listen on all interfaces so Apache can reverse-proxy to it
    app.run(host='0.0.0.0', port=5000)
