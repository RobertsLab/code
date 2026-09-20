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
  166:112437). Everything settable lives in one `## Variables` chunk under
  SETUP -- unit identity (`unit_label`, `unit_id`, `expected_uid`), input
  paths (`log_dir`, the unit's directory under `calibration/`, and
  `manifest_path`, the shared `calibration/calibration_metadata.csv`) and the
  analysis knobs, all as plain R assignments. To point a notebook at
  different files, edit that chunk. Each notebook checks the `rp2040_uid` on every data row against the UID it was
  written for and refuses to knit against another unit's data. Magnet count
  is a grouping variable within each notebook, so 1, 2 and 3 magnets are fit
  and compared side by side even though each was recorded to its own file.
  Writes `calibration-coefficients_Unit-X.csv`,
  `calibration-hold-summary_Unit-X.csv` and the rendered HTML to
  `outputs/calibration_unit-X/`. The coefficient file carries the ambient
  background vector (`bg_x_uT`, `bg_y_uT`, `bg_z_uT`) alongside the fitted
  `a`/`b`/`c`, because the curves are fit to background-corrected field and
  anything applying them has to correct its own readings the same way
  first.

  Each notebook is **self-contained**: data import, UID verification,
  metadata matching and QC, background subtraction, model fitting and
  comparison, plots and export all live in the file itself, so it can be
  knit on its own with no other `.Rmd` present. The three are identical
  apart from the title, the header comment and the `## Variables` chunk --
  which also means a change to the analysis has to be made in all three.
  `diff code/calibration_unit-A.Rmd code/calibration_unit-B.Rmd` should
  report only those two blocks.

- `calibration/calibration_metadata.csv`: the hand-written bench log that
  tells the calibration notebooks which stretch of which recording
  corresponds to which distance and magnet count. One row per hold, covering
  **all** units in one file -- each notebook filters it to its own `unit`:
  `unit` (`A`/`B`/`C`), `distance` (mm), `magnet.count` (`0` = a no-magnet
  background hold), `start.time`, `end.time` (both `YYYY-MM-DDTHH:MM:SS`,
  matching the RTC), and `data.file` (which `calibration-N.csv` the hold was
  recorded in). Samples are matched to holds on file *and* time, since two
  units recorded on the same day share wall-clock times.
  `calibration_manifest_template.csv` is the older per-unit form of this
  file, kept for reference.

- `gaping_analysis.Rmd`: R Markdown notebook that reads one or more
  `gapelog_NNN.csv` files (one per logging session/boot), reconstructs
  wall-clock timestamps for elapsed-time rows from user-supplied session
  start times, computes field magnitude, and detects gape events as
  sustained drops in magnitude below a threshold (configurable as
  noise-based, relative, or absolute). Field magnitude is
  **background-corrected** before anything else uses it: the per-axis ambient
  vector recorded at calibration is subtracted from the logged axes and the
  magnitude recomputed, so the curve is applied to the same quantity it was
  fit to. The corrected magnitude drives event detection as well as the
  distance conversion. With `apply_calibration <- FALSE` there is no
  background available, so magnitude falls back to raw and the notebook says
  so. Produces event tables (frequency,
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

A **hold** is the unit of measurement below: the magnet parked at one
distance, with one magnet count, held still while the logger samples at 1 Hz
for at least 60 s (~61 samples). One row of `calibration_metadata.csv` is one
hold, and each hold is averaged into a single point on the curve.

1. Fix the magnetometer; move the magnet(s) to each distance in turn
   (1, 2, ... 20 mm), holding each **still for at least 60 s** while the
   logger runs at 1 Hz -- ~61 samples per point.
2. **Log the holds only.** Stop logging while repositioning, so no
   adjustment period is ever recorded: each hold is a contiguous run of
   samples and the timestamp jumps to the next hold. Nothing then needs
   trimming off the ends of a hold, and every recorded sample is a genuine
   still-magnet reading. (Earlier recordings ran continuously and the
   notebooks discarded 10 s at the start and 5 s at the end of each hold;
   that trimming has been removed.)
3. Repeat the full sweep with 1, 2 and 3 magnets, one file per magnet count:
   `calibration-1.csv`, `calibration-2.csv`, `calibration-3.csv`.
4. Record a no-magnet background hold with the magnets well away from the
   sensor, as `calibration-0.csv`.
5. Write down the wall-clock start/end of every hold as you go, then
   transcribe them into `calibration/calibration_metadata.csv`.

That is ~61 holds and roughly two hours of recording per unit, split across
the four `calibration-N.csv` files in `calibration/unit-X/`.

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
- **Range.** 1-20 mm covers gape distance, but the sensor-to-magnet
  separation also includes whatever the magnet and sensor are mounted on.
  Extend the sweep further if the mounted geometry never actually gets
  down to a few mm.
