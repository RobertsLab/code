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

- `calibration_unit-A.Rmd`, `calibration_unit-B.Rmd`,
  `calibration_unit-C.Rmd`: one R Markdown notebook per logger unit, each
  turning a bench calibration recording into a field-strength-to-distance
  curve (following Vereycken et al. 2024, *Ecological Indicators*
  166:112437). Input paths are set in the `params:` block at the top of each
  file, or supplied at render time -- including interactively via
  `rmarkdown::render("calibration_unit-A.Rmd", params = "ask")`. Each
  notebook checks the `rp2040_uid` on every data row against the UID it was
  written for and refuses to knit against another unit's data. Magnet count
  is a grouping variable within each notebook, not a separate file, so 1, 2
  and 3 magnets are fit and compared side by side. Writes
  `calibration-coefficients_Unit-X.csv` and
  `calibration-hold-summary_Unit-X.csv` to `output/`.

- `_calibration-body.Rmd`: the shared analysis body pulled in as a knitr
  child document by all three unit notebooks -- data import, UID
  verification, manifest matching and QC, background subtraction, model
  fitting and comparison, plots, and export. Not knit directly (it has no
  YAML header). Edit this file to change the analysis for every unit at once.

- `calibration_manifest_template.csv`: template for the hand-written bench
  log that tells the calibration notebooks which stretch of the recording
  corresponds to which distance and magnet count. One row per hold:
  `n_magnets` (`0` = a no-magnet background hold), `distance_mm`,
  `start_time`, `end_time` (both `YYYY-MM-DDTHH:MM:SS`, matching the RTC),
  and free-text `notes`. Times are deliberately left blank -- the notebooks
  stop with an explanatory error until they are filled in.

- `gaping_analysis.Rmd`: R Markdown notebook that reads one or more
  `gapelog_NNN.csv` files (one per logging session/boot), reconstructs
  wall-clock timestamps for elapsed-time rows from user-supplied session
  start times, computes field magnitude, and detects gape events as
  sustained drops in magnitude below a threshold (configurable as
  noise-based, relative, or absolute). Produces event tables (frequency,
  duration, time between events) and plots (magnitude with detected
  events, per-axis time series, duration/interval distributions, events
  per hour). Optionally applies the calibration curves to convert field
  magnitude into **gape distance in mm** -- per-event maximum gape, mean
  gape, and gape amplitude relative to the closed-shell baseline, plus a
  distance time series. Set `apply_calibration <- FALSE` to skip that and
  work in field magnitude alone; every other result is unaffected either
  way. All user-configurable settings live in a single parameters chunk
  near the top of the notebook.

- `deployment_metadata_template.csv`: template recording which magnet count
  each `gapelog_NNN.csv` was deployed with (`source_file`, `n_magnets`,
  `oyster_id`, `notes`). Required by `gaping_analysis.Rmd` when
  `apply_calibration` is on -- the logger records the board UID but not the
  magnet count, so this file is the only place that association exists.

### CALIBRATION PROTOCOL

Each magnetometer needs its own calibration curve, and the curve is specific
to the magnet and geometry it was measured with -- re-calibrate after any
physical change to the magnet, its mounting, or the sensor-magnet layout.

1. Fix the magnetometer; move the magnet(s) to each distance in turn
   (0, 1, 2, ... 10 mm), holding each **still for 60 s** while the logger
   runs at 1 Hz. The notebooks discard 10 s at the start and 5 s at the end
   of each hold, leaving ~45 clean samples per point.
2. Leave a ~15-20 s pause between holds while repositioning, so the step
   boundaries are unambiguous in the recorded trace.
3. Repeat the full sweep with 1, 2 and 3 magnets.
4. Record a no-magnet background hold at the start and end of the session,
   with the magnets well away from the sensor.
5. Write down the wall-clock start/end of every hold as you go, then
   transcribe them into a copy of `calibration_manifest_template.csv`.

That is ~35 holds and roughly 45 minutes of recording per unit, in one
continuous session and so one `gapelog_NNN.csv` file.

**Magnet count is part of the calibration.** Different oysters may need
different magnet counts (leading-edge shell thickness, animal size), and more
magnets mean a stronger field at every distance -- so the same reading maps to
a different distance. The exported coefficients are keyed on
`(rp2040_uid, n_magnets, model)` for that reason. The logger records the board
UID but *not* the magnet count, so **the magnet count has to be written down
per deployment** and carried alongside the `gapelog_NNN.csv`; without it there
is no way to tell which calibration row applies. Applying a 1-magnet curve to
3-magnet data reads a true 10 mm gape as ~6 mm, and can return negative
distances at the closed end.

**Cross-check the magnet count when analysing a deployment.**
`gaping_analysis.Rmd` refuses to run if a log file has no metadata row, or if
no curve exists for that UID and magnet count. But a magnet count that is
*wrong yet calibrated* -- metadata says 1 magnet, 2 were deployed -- passes
every automated check, because the readings sit comfortably inside the
1-magnet curve's fitted range. It just compresses the distances: in testing,
a true 7.0 mm gape read as 4.9 mm and a 1.0 mm closed baseline read as
0.14 mm. Nothing in the data can recover the magnet count, so the notebook
prints what every calibrated magnet count would imply for that session's
closed and open field levels. Compare the `closed_mm` column against how far
the magnet actually sat from the sensor with the shell shut -- so measure and
write that down at deployment too.

Two caveats worth checking before committing to a magnet:

- **Saturation.** The TLV493D clips at +/- 130 mT. At 0-2 mm a stack of
  magnets can exceed that, and a clipped reading looks like a plateau rather
  than an error. The notebooks flag any hold with near-full-scale readings.
- **Range.** 0-10 mm covers gape distance, but the sensor-to-magnet
  separation also includes whatever the magnet and sensor are mounted on.
  Extend the sweep to 15-20 mm if the mounted geometry never actually gets
  down to a few mm.
