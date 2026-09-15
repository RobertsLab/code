`gape_monitoring`

### DESCRIPTION

This directory contains guides and code for implementing, and
analyzing, gape monitoring data collected using an Adafruit RP2040
Adalogger Feather and the Adafruit TLV493D three-axis magnetometer. All
components are STEMMA-QT compatible (i.e. uses a cable with a 4-pin
connector instead of requiring soldering).

Timestamping is handled by an Adafruit DS3231 Precision RTC (STEMMA QT)
when one is present: each row is stamped with a wall-clock
`YYYY-MM-DDTHH:MM:SS` timestamp read from the RTC. If the RTC is not
detected at boot (wiring issue, module not installed, etc.), the logger
does not halt -- it falls back to recording seconds elapsed since boot
(`time.monotonic()`) and logs an `RTC_INIT_ERROR` event so the gap is
visible after the fact. The analysis notebook reconstructs absolute
timestamps for elapsed-time rows from a start time recorded by hand at
deployment, or falls back to relative time if none is supplied.

### FILES

- `code.py`: CircuitPython script for an Adafruit Feather RP2040 Adalogger
  with a TLV493D magnetometer (STEMMA QT/I2C) and, optionally, a DS3231
  RTC (STEMMA QT, shares the I2C bus with the TLV493D). Deploy by saving
  it as `code.py` on the CIRCUITPY drive. Samples the magnetometer once
  per second and appends rows to a new `gapelog_NNN.csv` on the SD card
  each time the board boots, so repeated deployments never overwrite
  earlier data. Boot events and error conditions (SD card, RTC, sensor,
  or write failures) are also appended to a persistent `eventlog.csv`
  that accumulates across all deployments. An onboard NeoPixel blinks
  green on each successful sample and red on errors.

- `gaping_analysis.Rmd`: R Markdown notebook that reads one or more
  `gapelog_NNN.csv` files (one per logging session/boot), reconstructs
  wall-clock timestamps for elapsed-time rows from user-supplied session
  start times, computes field magnitude, and detects gape events as
  sustained drops in magnitude below a threshold (configurable as
  noise-based, relative, or absolute). Produces event tables (frequency,
  duration, time between events) and plots (magnitude with detected
  events, per-axis time series, duration/interval distributions, events
  per hour). All user-configurable settings live in a single parameters
  chunk near the top of the notebook.
