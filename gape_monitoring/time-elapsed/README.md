`gape_monitoring/time-elapsed`

## DESCRIPTION

Oyster valve-gaping monitoring using a three-axis magnetometer: one magnet
is mounted on each shell valve, and field magnitude falls as the shell
opens (magnet moves away from the sensor) and rises as it closes. This
directory holds the CircuitPython datalogger that runs on the sensor board
and the R Markdown notebook that turns its logs into gape events and plots.

This is the "time-elapsed" variant: the board has no real-time clock, so
the logger records seconds elapsed since boot rather than wall-clock
timestamps. The analysis notebook reconstructs absolute timestamps from a
start time recorded by hand at deployment, or falls back to relative time
if none is supplied. A future revision with an onboard RTC will write
timestamps directly and won't need this reconstruction step.

### FILES

- `code.py`: CircuitPython script for an Adafruit Feather RP2040 Adalogger
  with a TLV493D magnetometer (STEMMA QT/I2C). Deploy by saving it as
  `code.py` on the CIRCUITPY drive. Samples the magnetometer once per
  second and appends `elapsed_s,x_uT,y_uT,z_uT` rows to a new
  `gapelog_NNN.csv` on the SD card each time the board boots, so repeated
  deployments never overwrite earlier data. Boot events and error
  conditions (SD card, sensor, or write failures) are also appended to a
  persistent `eventlog.csv` that accumulates across all deployments. An
  onboard NeoPixel blinks green on each successful sample and red on
  errors.

- `gaping_analysis-time-elapsed.Rmd`: R Markdown notebook that reads one or
  more `gapelog_NNN.csv` files (one per logging session/boot), optionally
  reconstructs wall-clock timestamps from user-supplied session start
  times, computes field magnitude, and detects gape events as sustained
  drops in magnitude below a threshold (configurable as noise-based,
  relative, or absolute). Produces event tables (frequency, duration,
  time between events) and plots (magnitude with detected events, per-axis
  time series, duration/interval distributions, events per hour). All
  user-configurable settings live in a single parameters chunk near the
  top of the notebook.
