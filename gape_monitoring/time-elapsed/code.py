"""
oyster_gape_logger.py  ->  save this file as code.py on CIRCUITPY

Logs three-axis magnetometer readings from an Adafruit TLV493D
(STEMMA QT, connected via I2C) to a CSV file on the SD card of an
Adafruit Feather RP2040 Adalogger, for oyster valve-gaping detection.

SD card mounting: NOT done in this script. The Feather RP2040 Adalogger
is one of the boards CircuitPython automatically mounts a present SD
card for at boot, to /sd, before code.py ever runs (no storage.mount()
call needed here). This was confirmed on the actual hardware: a bare
open("/sd/test.txt", "w") succeeded with no mounting code at all.

Timestamping: NO real-time clock is installed yet. Each row records
seconds elapsed since this script started (time.monotonic()), not a
wall-clock date/time. Write down the actual start date/time by hand
when you deploy the logger -- the companion R Markdown file adds it
back in to reconstruct real timestamps. An RTC module (e.g. PCF8523 /
DS3231) is planned for a future revision; once installed, replace the
elapsed-time column with rtc.datetime directly.

A new log file is created each time the board (re)boots, named
gapelog_001.csv, gapelog_002.csv, etc., so re-running the script never
overwrites previous deployment data already on the card.

Event logging: every boot and every error condition (SD card, sensor,
or write failures) is appended to a single persistent /sd/eventlog.csv
that accumulates across all deployments -- it is never recreated or
truncated. Each row is tagged with a boot number (from /sd/boot_count.txt,
incremented once per boot) and the seconds elapsed since that boot's
script start, so events from different sessions can be told apart even
without an RTC. Note there is no way to detect clean power-OFF in
CircuitPython (power just disappears) -- only power-ON (boot) events
are observable, and that's what gets logged.
"""

import time
import board
import neopixel
import adafruit_tlv493d

# ----------------------------------------------------------------------
# CONFIGURATION -- tune these without touching the rest of the script
# ----------------------------------------------------------------------
SAMPLE_INTERVAL = 1.0     # seconds between magnetometer readings
SD_MOUNT_POINT = "/sd"
LOG_PREFIX = "gapelog_"
LOG_SUFFIX = ".csv"
FLUSH_EVERY_N_ROWS = 1    # flush to disk after every N rows (1 = safest)
EVENT_LOG_PATH = SD_MOUNT_POINT + "/eventlog.csv"
BOOT_COUNT_PATH = SD_MOUNT_POINT + "/boot_count.txt"

boot_start_time = time.monotonic()

# ----------------------------------------------------------------------
# Onboard NeoPixel status LED (optional heartbeat; safe no-op if board
# has none). Requires neopixel.mpy in /lib on CIRCUITPY.
#
# Green = a sample was read and written to the log (data is being
# collected). Red = an error condition (SD card, sensor, or write
# failure) -- see call sites below.
# ----------------------------------------------------------------------
COLOR_OFF = (0, 0, 0)
COLOR_LOGGING = (0, 255, 0)
COLOR_ERROR = (255, 0, 0)

try:
    pixel = neopixel.NeoPixel(board.NEOPIXEL, 1, brightness=0.2, auto_write=True)
    pixel.fill(COLOR_OFF)
except AttributeError:
    pixel = None


def blink(times=1, duration=0.05, color=COLOR_LOGGING):
    if pixel is None:
        return
    for _ in range(times):
        pixel.fill(color)
        time.sleep(duration)
        pixel.fill(COLOR_OFF)
        time.sleep(duration)


# ----------------------------------------------------------------------
# Confirm the SD card is present and writable (it should already be
# auto-mounted at /sd by CircuitPython itself before this script runs
# -- see module docstring). This just verifies the mount succeeded
# rather than performing it.
# ----------------------------------------------------------------------
def check_sd_card():
    import os
    try:
        os.listdir(SD_MOUNT_POINT)
    except OSError as e:
        raise OSError(
            "SD card not found at {}. Check card is seated/formatted "
            "FAT32, or power-cycle the board. ({})".format(SD_MOUNT_POINT, e)
        )


try:
    check_sd_card()
    print("SD card confirmed at {}".format(SD_MOUNT_POINT))
except OSError as e:
    # Can't persist this one -- the SD card itself is the thing that failed.
    print("FATAL: {}".format(e))
    while True:
        blink(3, 0.1, color=COLOR_ERROR)
        time.sleep(1)

# ----------------------------------------------------------------------
# Persistent event log (boot events + errors), separate from the
# per-boot sensor-data CSV. Unlike gapelog_NNN.csv, this file is never
# recreated -- it accumulates across every deployment so boot/error
# history survives reboots.
# ----------------------------------------------------------------------
def next_boot_number():
    n = 0
    try:
        with open(BOOT_COUNT_PATH, "r") as f:
            n = int(f.read().strip())
    except (OSError, ValueError):
        n = 0
    n += 1
    try:
        with open(BOOT_COUNT_PATH, "w") as f:
            f.write(str(n))
    except OSError:
        pass
    return n


boot_number = next_boot_number()

try:
    with open(EVENT_LOG_PATH, "r"):
        pass
except OSError:
    with open(EVENT_LOG_PATH, "w") as f:
        f.write("boot_number,elapsed_s,event,detail\n")


def log_event(event, detail=""):
    elapsed = time.monotonic() - boot_start_time
    row = "{},{:.3f},{},{}\n".format(boot_number, elapsed, event, detail)
    try:
        with open(EVENT_LOG_PATH, "a") as f:
            f.write(row)
    except OSError as e:
        print("Event log write failed: {}".format(e))


log_event("BOOT", "SD card confirmed at {}".format(SD_MOUNT_POINT))

# ----------------------------------------------------------------------
# Initialize the TLV493D magnetometer (I2C / STEMMA QT)
# ----------------------------------------------------------------------
try:
    i2c = board.I2C()
    tlv = adafruit_tlv493d.TLV493D(i2c)
    print("TLV493D magnetometer initialized")
    log_event("SENSOR_INIT_OK", "TLV493D magnetometer initialized")
except (OSError, ValueError, RuntimeError) as e:
    print("FATAL: could not initialize TLV493D: {}".format(e))
    print("Check wiring on the I2C / STEMMA QT bus.")
    log_event("SENSOR_INIT_ERROR", str(e))
    while True:
        blink(5, 0.1, color=COLOR_ERROR)
        time.sleep(1)

# ----------------------------------------------------------------------
# Pick a log filename that doesn't collide with existing files
# ----------------------------------------------------------------------
def next_log_path():
    existing = set(os_listdir_safe(SD_MOUNT_POINT))
    n = 1
    while True:
        name = "{}{:03d}{}".format(LOG_PREFIX, n, LOG_SUFFIX)
        if name not in existing:
            return SD_MOUNT_POINT + "/" + name
        n += 1


def os_listdir_safe(path):
    import os
    try:
        return os.listdir(path)
    except OSError:
        return []


log_path = next_log_path()
print("Logging to {}".format(log_path))
log_event("LOG_START", "Logging sensor data to {}".format(log_path))

with open(log_path, "w") as f:
    f.write("elapsed_s,x_uT,y_uT,z_uT\n")

# ----------------------------------------------------------------------
# Main logging loop
# ----------------------------------------------------------------------
start_time = time.monotonic()
row_count = 0

print("Starting logging loop. Sample interval: {} s".format(SAMPLE_INTERVAL))
print("NOTE: elapsed_s is seconds since this script started, "
      "NOT a wall-clock timestamp (no RTC installed).")

while True:
    loop_start = time.monotonic()
    elapsed = loop_start - start_time

    try:
        x, y, z = tlv.magnetic
    except (OSError, RuntimeError) as e:
        print("Sensor read failed: {}".format(e))
        log_event("SENSOR_READ_ERROR", str(e))
        time.sleep(SAMPLE_INTERVAL)
        continue

    row = "{:.3f},{:.3f},{:.3f},{:.3f}\n".format(elapsed, x, y, z)

    # Echo each reading to the serial console as it's captured -- useful
    # for testing over USB/REPL. Remove or comment out this line for
    # unattended field deployment if you want a quieter console.
    print("elapsed_s={:.3f}  x_uT={:.3f}  y_uT={:.3f}  z_uT={:.3f}".format(
        elapsed, x, y, z))

    try:
        with open(log_path, "a") as f:
            f.write(row)
            row_count += 1
            if row_count % FLUSH_EVERY_N_ROWS == 0:
                f.flush()
    except OSError as e:
        print("Write failed (card removed/full?): {}".format(e))
        log_event("WRITE_ERROR", str(e))
        blink(4, 0.1, color=COLOR_ERROR)
        time.sleep(SAMPLE_INTERVAL)
        continue

    blink(1, 0.02)

    # Sleep for the remainder of the interval, accounting for how long
    # the read+write above took, so the sample rate stays close to
    # SAMPLE_INTERVAL even as work per loop varies.
    elapsed_this_loop = time.monotonic() - loop_start
    sleep_time = SAMPLE_INTERVAL - elapsed_this_loop
    if sleep_time > 0:
        time.sleep(sleep_time)
