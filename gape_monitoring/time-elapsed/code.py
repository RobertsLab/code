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
"""

import time
import board
import digitalio
import adafruit_tlv493d

# ----------------------------------------------------------------------
# CONFIGURATION -- tune these without touching the rest of the script
# ----------------------------------------------------------------------
SAMPLE_INTERVAL = 1.0     # seconds between magnetometer readings
SD_MOUNT_POINT = "/sd"
LOG_PREFIX = "gapelog_"
LOG_SUFFIX = ".csv"
FLUSH_EVERY_N_ROWS = 1    # flush to disk after every N rows (1 = safest)

# ----------------------------------------------------------------------
# Onboard status LED (optional heartbeat; safe no-op if board has none)
# ----------------------------------------------------------------------
try:
    led = digitalio.DigitalInOut(board.LED)
    led.direction = digitalio.Direction.OUTPUT
except AttributeError:
    led = None


def blink(times=1, duration=0.05):
    if led is None:
        return
    for _ in range(times):
        led.value = True
        time.sleep(duration)
        led.value = False
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
    print("FATAL: {}".format(e))
    while True:
        blink(3, 0.1)
        time.sleep(1)

# ----------------------------------------------------------------------
# Initialize the TLV493D magnetometer (I2C / STEMMA QT)
# ----------------------------------------------------------------------
try:
    i2c = board.I2C()
    tlv = adafruit_tlv493d.TLV493D(i2c)
    print("TLV493D magnetometer initialized")
except (OSError, ValueError, RuntimeError) as e:
    print("FATAL: could not initialize TLV493D: {}".format(e))
    print("Check wiring on the I2C / STEMMA QT bus.")
    while True:
        blink(5, 0.1)
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
        blink(4, 0.1)
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
