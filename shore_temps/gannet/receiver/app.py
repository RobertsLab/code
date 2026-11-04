"""Receiver for SHORE tank temperature hubs.

Each hub POSTs its readings to /ingest with "Authorization: Bearer <token>".
The body is one reading or a list of readings (a backlog after a network outage):

    {"hub": "shore-temps-a", "uptime_s": 61, "epoch": 1793642085,
     "probes": [{"id": "28FF4A1B6E180301", "temp_c": 20.81}, ...]}

"epoch" (Unix time the reading was taken, from NTP) is optional.
Every probe reading becomes one row in DATA_DIR/raw/<hub>/<YYYY-MM-DD>.csv,
dated by when the server received it (UTC). Data is stored raw: no tank
names or offsets are applied here.
"""

import csv
import hmac
import json
import os
import re
import threading
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

RAW_DIR = Path(os.environ.get("DATA_DIR", "/data")) / "raw"
PORT = int(os.environ.get("PORT", "8080"))
MAX_BODY = 256 * 1024                       # bytes; a large backlog is still well under this
FIELDS = ["received_utc", "reading_utc", "hub", "uptime_s", "probe_id", "temp_c"]
PROBE_ID = re.compile(r"^[0-9A-F]{16}$")
MIN_EPOCH = 1704067200                      # 2024-01-01; anything earlier means the hub's clock isn't set

# HUB_TOKENS="shore-temps-a:<token>,shore-temps-b:<token>"
HUB_TOKENS = {}
for pair in os.environ["HUB_TOKENS"].split(","):
    hub, token = pair.strip().split(":", 1)
    if len(token) < 20:
        raise SystemExit(f"Token for {hub} is too short; generate a long random one.")
    HUB_TOKENS[hub] = token

write_lock = threading.Lock()


def hub_for_token(token):
    for hub, expected in HUB_TOKENS.items():
        if hmac.compare_digest(token.encode(), expected.encode()):
            return hub
    return None


def iso_utc(ts):
    return datetime.fromtimestamp(ts, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def to_rows(reading, hub, received):
    """Validate one reading and turn it into CSV rows. Raises ValueError if malformed."""
    if not isinstance(reading, dict) or not isinstance(reading.get("probes"), list):
        raise ValueError("reading must be an object with a 'probes' list")
    if reading.get("hub") != hub:
        raise ValueError(f"hub name {reading.get('hub')!r} does not match token for {hub!r}")
    epoch = reading.get("epoch")
    reading_utc = iso_utc(epoch) if isinstance(epoch, (int, float)) and epoch >= MIN_EPOCH else ""
    uptime = reading.get("uptime_s", "")
    rows = []
    for p in reading["probes"]:
        probe_id = p.get("id") if isinstance(p, dict) else None
        if not isinstance(probe_id, str) or not PROBE_ID.match(probe_id):
            raise ValueError(f"bad probe id {probe_id!r}")
        temp = p.get("temp_c")
        if temp is not None and not isinstance(temp, (int, float)):
            raise ValueError(f"bad temp_c {temp!r}")
        rows.append([received, reading_utc, hub, uptime, probe_id, "" if temp is None else temp])
    return rows


def append_rows(hub, rows):
    hub_dir = RAW_DIR / hub
    hub_dir.mkdir(parents=True, exist_ok=True)
    path = hub_dir / f"{datetime.now(timezone.utc):%Y-%m-%d}.csv"
    with write_lock:
        new_file = not path.exists()
        with path.open("a", newline="") as f:
            w = csv.writer(f)
            if new_file:
                w.writerow(FIELDS)
            w.writerows(rows)


class Handler(BaseHTTPRequestHandler):
    timeout = 15                            # drop clients that stall mid-request

    def reply(self, status, body):
        data = (json.dumps(body) + "\n").encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path == "/health":
            self.reply(200, {"status": "ok"})
        else:
            self.reply(404, {"error": "not found"})

    def do_POST(self):
        if self.path != "/ingest":
            return self.reply(404, {"error": "not found"})

        auth = self.headers.get("Authorization", "")
        hub = hub_for_token(auth[7:]) if auth.startswith("Bearer ") else None
        if hub is None:
            return self.reply(401, {"error": "bad or missing token"})

        length = int(self.headers.get("Content-Length") or 0)
        if not 0 < length <= MAX_BODY:
            return self.reply(413, {"error": f"body must be 1-{MAX_BODY} bytes"})

        received = iso_utc(datetime.now(timezone.utc).timestamp())
        try:
            body = json.loads(self.rfile.read(length))
            readings = body if isinstance(body, list) else [body]
            rows = [row for r in readings for row in to_rows(r, hub, received)]
        except ValueError as e:             # json.JSONDecodeError is a ValueError
            self.log_message("rejected from %s: %s", hub, e)
            return self.reply(400, {"error": str(e)})

        append_rows(hub, rows)
        self.reply(200, {"stored_rows": len(rows)})

    def log_request(self, code="-", size="-"):
        if not str(code).startswith("2"):  # keep the log to problems only
            super().log_request(code, size)


if __name__ == "__main__":
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Listening on port {PORT}; hubs: {', '.join(HUB_TOKENS)}; writing to {RAW_DIR}", flush=True)
    ThreadingHTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
