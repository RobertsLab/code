SHORE TANK TEMPERATURE RECEIVER (Gannet / Synology)
====================================================

The ESP32 temperature hubs send their probe readings to Gannet over HTTPS
once a minute. The receiver in this directory runs as a Docker container on
Gannet. It saves every reading to CSV files in data/.

    hub a ─┐   HTTPS POST         Gannet (Synology DSM)
    hub b ─┴──────────────────►   DSM reverse proxy (HTTPS, port 8443)
                                    └► receiver container (port 8080, Gannet only)
                                         └► data/raw/<hub>/<YYYY-MM-DD>.csv

PHP is not used. The receiver needs only Python's standard library, so there
is nothing to install with pip.


FILES
-----

  docker-compose.yml   Project definition that Container Manager reads.
  receiver/app.py      The receiver program.
  receiver/Dockerfile  Builds the container image (python:3.12-slim).
  .env.example         Template for the secrets file (.env).
  .gitignore           Keeps .env and data/ out of git.

Created on Gannet (never committed to git):

  .env                 Your DSM user ID and the hub tokens.
  data/                The temperature data.


REQUIREMENTS
------------

  - DSM 7.2 or later with the Container Manager package installed
    (Package Center). Check that Gannet's model supports Container Manager.
  - An admin account on Gannet, for the reverse proxy and firewall settings.
  - SSH access to Gannet. This is only needed to look up your user ID and to
    test the receiver. To turn SSH on: Control Panel > Terminal & SNMP.


SETUP
-----

1. Copy the files to Gannet

   Choose a folder on a shared folder, for example:

       /volume1/<share>/shore-temps/

   Use File Station or scp to copy everything in this gannet/ directory into
   that folder, including the hidden .env.example. Then create an empty
   subfolder called data. The result should look like this:

       shore-temps/
         docker-compose.yml
         .env.example
         receiver/app.py
         receiver/Dockerfile
         data/

2. Look up your user and group IDs

   The container runs as an ordinary DSM user, not as root. Use the user who
   owns the shore-temps folder. Over SSH, logged in as that user, run:

       id

   Write down the numbers after "uid=" and "gid=". On Synology the uid is
   usually 1026 or higher, and the gid is often 100 (the "users" group).

3. Create .env

   In the shore-temps folder:

       cp .env.example .env
       chmod 600 .env

   Edit .env and fill in:

       PUID=<uid from step 2>
       PGID=<gid from step 2>
       HUB_TOKENS=shore-temps-a:<token a>,shore-temps-b:<token b>

   Give each hub its own long random token. Generate one with:

       python3 -c "import secrets; print(secrets.token_urlsafe(32))"

   The hub names must match HUB_NAME in each hub's firmware, and each hub's
   firmware must contain its own token. Tokens shorter than 20 characters are
   refused. Never commit .env or put the tokens anywhere public.

4. Build and start the container

   In Container Manager:

     a. Project > Create
     b. Project name: shore-temps
     c. Path: the shore-temps folder from step 1
     d. Source: "Use existing docker-compose.yml"
     e. Skip the Web Station / web portal settings. They aren't needed,
        because the reverse proxy in step 6 handles HTTPS.
     f. Finish, and let it build and start the project.

   Open the container's log (Container Manager > Container >
   shore-temps-receiver > Log). You should see:

       Listening on port 8080; hubs: shore-temps-a, shore-temps-b; writing to /data/raw

   The container restarts automatically after a reboot or crash
   ("restart: unless-stopped").

5. Test on Gannet

   For security, the container only accepts connections from Gannet itself.
   Over SSH, run:

       curl http://127.0.0.1:8080/health

   Expected reply: {"status": "ok"}

6. Set up the HTTPS reverse proxy

   Go to Control Panel > Login Portal > Advanced > Reverse Proxy > Create:

       Source:       Protocol HTTPS, Hostname gannet.fish.washington.edu, Port 8443
       Destination:  Protocol HTTP,  Hostname localhost,                  Port 8080

   Then go to Control Panel > Security > Certificate > Settings and assign
   Gannet's certificate to the new reverse proxy entry.

   You can use a port other than 8443. If you do, use it in steps 7 and 8
   and in the hub firmware too.

7. Open the firewall

   If DSM's firewall is turned on (Control Panel > Security > Firewall),
   add a rule that allows TCP port 8443.

8. Test from the hubs' network

   On a phone or laptop joined to UW MPSK (the network the hubs use), open:

       https://gannet.fish.washington.edu:8443/health

   You should see {"status": "ok"}. If the page times out, the port is being
   blocked somewhere between UW MPSK and Gannet. Check the DSM firewall
   first, then ask UW IT whether the port is open.

   To test a full upload from any machine (replace <token a>):

       curl -H "Authorization: Bearer <token a>" \
            -d '{"hub":"shore-temps-a","uptime_s":1,"probes":[{"id":"28CC7DBC79250B81","temp_c":20.0}]}' \
            https://gannet.fish.washington.edu:8443/ingest

   Expected reply: {"stored_rows": 1}. Then delete the test row from
   today's file in data/raw/shore-temps-a/.


DATA
----

There is one file per hub per day. The date is the UTC date the reading
arrived:

    data/raw/shore-temps-a/2026-11-04.csv

    received_utc,reading_utc,hub,uptime_s,probe_id,temp_c
    2026-11-04T17:28:59Z,,shore-temps-a,61,28CC7DBC79250B81,18.69

  received_utc  When Gannet received the reading (UTC).
  reading_utc   When the hub took the reading (UTC). Empty if the hub sent no
                timestamp. Hubs only send one once their clock has been set
                from the network (NTP).
  hub           Hub name. Each token can only write data for its own hub.
  uptime_s      Seconds since the hub last started. Use it to spot duplicate
                readings and reboots.
  probe_id      The DS18B20's unique 64-bit ID. Tank names are NOT applied
                here.
  temp_c        Temperature in °C, uncalibrated. Empty if the probe failed
                to read.

These files are raw data. Apply the probe ID > tank map and the ice-bath
offsets in a separate processing step, and never edit these files in place.
Include data/ in Gannet's snapshot or backup schedule.


PROBE IDS
---------

Each probe's number and its ID, as it appears in the probe_id column:

  | Probe | Probe ID         |
  |------:|------------------|
  |     1 | 285453BB79250B4D |
  |     2 | 285963BE79250B02 |
  |     3 | 28F953870074746A |
  |     4 | 283EBEB579250B62 |
  |     5 | 289C2BB779250B32 |
  |     6 | 28F55BC479250B2A |
  |     7 | 28B27FB379250B13 |
  |     8 | 28CC7DBC79250B81 |
  |     9 | 282107B179250BAF |
  |    10 | 286B0FB679250BAD |


WHAT THE RECEIVER ACCEPTS
-------------------------

  GET  /health   No token needed. Returns {"status": "ok"}.
  POST /ingest   Header: Authorization: Bearer <token>
                 Body: one reading or a list of readings (a backlog), as
                 JSON, up to 256 KB:

                 {"hub": "shore-temps-a", "uptime_s": 61, "epoch": 1793642085,
                  "probes": [{"id": "28FF4A1B6E180301", "temp_c": 20.81}]}

                 "epoch" (the Unix time of the reading) is optional.
                 "temp_c" may be null.

  Replies: 200 stored, 400 bad data, 401 bad or missing token,
           413 body empty or too large.


MAINTENANCE
-----------

Change a token or add a hub
  Edit HUB_TOKENS in .env, then go to Container Manager > Project >
  shore-temps > Action > Build. Also update the token in that hub's firmware
  and reflash it.

Update the receiver code
  Copy the new receiver/app.py into place, then go to Container Manager >
  Project > shore-temps > Action > Build.

Stop the receiver
  Container Manager > Project > shore-temps > Action > Stop. The hubs keep
  their readings and resend them later, up to the limit of their memory.

Logs
  Container Manager > Container > shore-temps-receiver > Log. Only rejected
  requests and errors are logged; successful uploads are not. The log is
  capped at about 3 MB.


TROUBLESHOOTING
---------------

  Container stops at startup: "HUB_TOKENS must look like ..."
      .env is missing, doesn't define HUB_TOKENS, or the format is wrong
      (step 3).

  Container stops at startup: "Token for ... is too short"
      Generate a longer token (step 3).

  Container log shows "PermissionError" under /data
      PUID/PGID don't match a user who can write to data/. Check step 2,
      and check the folder's permissions in File Station.

  Hub gets 401
      The hub's token doesn't match HUB_TOKENS, or the container wasn't
      rebuilt after .env was edited.

  Hub gets 400 "hub name ... does not match token"
      The hub's HUB_NAME doesn't match the name paired with its token in
      .env.

  /health works on Gannet but not from UW MPSK
      Reverse proxy, certificate, or firewall problem (steps 6-8).
