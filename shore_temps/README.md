# shore_temps

Temperature monitoring for the SHORE tanks. Each ESP32 "hub" reads a set of
DS18B20 temperature probes once a minute and uploads the readings over HTTPS
to a receiver on Gannet, which saves them to daily CSV files.

```
DS18B20 probes ──1-Wire──► ESP32 hub ──WiFi (UW MPSK), HTTPS POST──► Gannet receiver ──► data/raw/<hub>/<YYYY-MM-DD>.csv
```

There are two hubs, `shore-temps-a` and `shore-temps-b`. Both run the same
firmware; only `HUB_NAME` and `HUB_TOKEN` differ.

## Repository layout

| Path                      | What it is                                                         |
|---------------------------|--------------------------------------------------------------------|
| `src/main.cpp`            | Hub firmware (PlatformIO, Arduino framework)                       |
| `include/secrets.h.example` | Template for `secrets.h` (WiFi credentials and hub token)        |
| `include/gannet_ca.h`     | Root certificate the hub uses to verify Gannet's HTTPS certificate |
| `platformio.ini`          | Board, serial port, and library settings                           |
| `gannet/`                 | Receiver that runs on Gannet in Docker. Setup steps are in [gannet/readme.txt](gannet/readme.txt) |

## Hardware

- ESP32 board (`esp32dev` target; the hubs are Adafruit HUZZAH32s)
- DS18B20 probes, all on one 1-Wire bus on GPIO 4 (pin A5 on the HUZZAH32)

## How the hub works

- Reads every probe on the bus once a minute at 12-bit resolution. The bus is
  re-scanned each time, so a probe that is unplugged and plugged back in
  shows up again without a reboot.
- Sends each probe's ID and temperature in °C. A probe that fails to read is
  sent as `null`.
- Sets its clock from NTP (UTC) and adds a timestamp to each reading. Until
  the clock is set, readings go out without one, and Gannet records only
  when it received them.
- If Gannet can't be reached, keeps up to 180 readings (3 hours) in memory
  and sends them once the connection is back. The oldest are dropped when
  the backlog is full, and the whole backlog is lost if the hub reboots.
- Reboots itself if WiFi has been down for an hour.
- Prints its MAC address, WiFi diagnostics, and every reading to the serial
  monitor (115200 baud).

## Setting up a hub

Requires [PlatformIO](https://platformio.org/) (for example, the VS Code
extension).

1. Copy `include/secrets.h.example` to `include/secrets.h` and fill in
   `WIFI_SSID`, `WIFI_PASS`, and `HUB_TOKEN`. `secrets.h` is git-ignored;
   never commit it.
2. In `src/main.cpp`, set `HUB_NAME` to this hub's name (`shore-temps-a` or
   `shore-temps-b`). The name and token must match a pair in `HUB_TOKENS` on
   Gannet.
3. Connect the board over USB. `platformio.ini` expects `/dev/ttyUSB0`;
   change `upload_port` and `monitor_port` if yours is different.
4. Build, upload, and open the serial monitor:

   ```
   pio run -t upload
   pio device monitor
   ```

5. Check the serial output for `Uploaded 1 reading(s)`. If WiFi won't
   connect, the output shows whether the network was found and why the
   connection failed. UW IT may ask for the MAC address printed at startup.

## Data

Readings land on Gannet in `data/raw/<hub>/<YYYY-MM-DD>.csv`, one file per
hub per UTC day. These files are raw: probes are identified only by ID, and
temperatures are uncalibrated. Use the table below to map IDs to probe
numbers. See [gannet/readme.txt](gannet/readme.txt) for the column
definitions.

## Probe IDs

Each DS18B20 probe has a unique 64-bit ID. This ID is the `probe_id` in the
data files.

| Probe | Probe ID           |
|------:|--------------------|
|     1 | `285453BB79250B4D` |
|     2 | `285963BE79250B02` |
|     3 | `28F953870074746A` |
|     4 | `283EBEB579250B62` |
|     5 | `289C2BB779250B32` |
|     6 | `28F55BC479250B2A` |
|     7 | `28B27FB379250B13` |
|     8 | `28CC7DBC79250B81` |
|     9 | `282107B179250BAF` |
|    10 | `286B0FB679250BAD` |

## Troubleshooting

| Problem | Where to look |
|---------|---------------|
| Hub can't join WiFi | Serial output: whether the network was found and the disconnect reason |
| Uploads fail with a certificate error | Gannet's certificate may now chain to a different root. See the comment at the top of `include/gannet_ca.h` |
| Hub gets HTTP 401 or 400 | Token or hub name doesn't match `HUB_TOKENS` on Gannet. See the troubleshooting section in [gannet/readme.txt](gannet/readme.txt) |
| A probe is missing or reads `null` | Check its wiring. The bus is re-scanned every minute, so a reconnected probe comes back on its own |
