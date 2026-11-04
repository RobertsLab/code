#include <Arduino.h>
// shore_temps — SHORE tank temperature hub (ESP32 + DS18B20 probes on one 1-Wire bus)
// Reads the probes once a minute and POSTs the readings to the receiver on Gannet (see gannet/readme.txt).
#include <WiFi.h>
#include <NetworkClientSecure.h>
#include <HTTPClient.h>
#include <OneWire.h>
#include <DallasTemperature.h>
#include <esp_mac.h>
#include <time.h>
#include <deque>
#include "secrets.h"                             // WIFI_SSID, WIFI_PASS, HUB_TOKEN; copy from secrets.h.example
#include "gannet_ca.h"                           // root certificate for Gannet's HTTPS

const char* HUB_NAME  = "shore-temps-b";        // "shore-temps-b" on the second hub; must match HUB_TOKENS on Gannet
const int   ONE_WIRE_PIN = 4;                    // GPIO 4 (A5 on HUZZAH32)
const char* RECEIVER_URL = "https://gannet.fish.washington.edu:8443/ingest";

const unsigned long READ_EVERY_MS   = 60UL * 1000;       // read probes once a minute
const unsigned long WIFI_GIVEUP_MS  = 60UL * 60 * 1000;  // reboot if WiFi is down this long (the backlog is lost on reboot)
const size_t MAX_BACKLOG  = 180;                 // readings held while Gannet is unreachable (3 h at one per minute)
const size_t MAX_PER_POST = 30;                  // readings per upload when catching up on a backlog
const time_t MIN_EPOCH    = 1704067200;          // 2024-01-01; an earlier clock means NTP hasn't set it yet

OneWire oneWire(ONE_WIRE_PIN);
DallasTemperature sensors(&oneWire);
NetworkClientSecure tls;

std::deque<String> backlog;                      // readings not yet accepted by Gannet, oldest first
unsigned long lastRead = 0, wifiLostAt = 0;

String idToString(const DeviceAddress a) {
  char buf[17];
  for (int i = 0; i < 8; i++) sprintf(buf + 2 * i, "%02X", a[i]);
  return String(buf);
}

String readProbes() {
  sensors.begin();                               // re-scan so a replugged probe reappears
  sensors.requestTemperatures();                 // ~750 ms at 12-bit resolution
  String out = "{\"hub\":\"" + String(HUB_NAME) + "\",\"uptime_s\":" + String(millis() / 1000);
  time_t now = time(nullptr);
  if (now >= MIN_EPOCH) out += ",\"epoch\":" + String((unsigned long)now);   // omitted until NTP has set the clock
  out += ",\"probes\":[";
  DeviceAddress a;
  int n = sensors.getDeviceCount(), written = 0;
  for (int i = 0; i < n; i++) {
    if (!sensors.getAddress(a, i)) continue;
    float c = sensors.getTempC(a);
    bool bad = (c == DEVICE_DISCONNECTED_C) || (c == 85.0);   // 85.0 = power-on value, not a reading
    if (written++) out += ",";
    out += "{\"id\":\"" + idToString(a) + "\",\"temp_c\":" + (bad ? String("null") : String(c, 2)) + "}";
  }
  out += "]}";
  Serial.println(out);
  return out;
}

void queueReading(const String& reading) {
  if (backlog.size() >= MAX_BACKLOG) {
    backlog.pop_front();
    Serial.println("Backlog full: dropped the oldest reading");
  }
  backlog.push_back(reading);
}

// Send queued readings to Gannet, oldest first, until the backlog is empty or an upload fails.
void uploadBacklog() {
  while (!backlog.empty() && WiFi.status() == WL_CONNECTED) {
    size_t n = min(backlog.size(), MAX_PER_POST);
    String body = "[";
    for (size_t i = 0; i < n; i++) {
      if (i) body += ",";
      body += backlog[i];
    }
    body += "]";

    HTTPClient http;
    http.begin(tls, RECEIVER_URL);
    http.setTimeout(15000);
    http.addHeader("Content-Type", "application/json");
    http.addHeader("Authorization", String("Bearer ") + HUB_TOKEN);
    int code = http.POST(body);
    String reply = code > 0 ? http.getString() : HTTPClient::errorToString(code);
    http.end();
    reply.trim();

    if (code == 200) {
      backlog.erase(backlog.begin(), backlog.begin() + n);
      Serial.printf("Uploaded %u reading(s); %u waiting\n", n, backlog.size());
    } else if (code == 400 || code == 413) {
      backlog.erase(backlog.begin(), backlog.begin() + n);   // Gannet will never accept these; don't let them block newer ones
      Serial.printf("Gannet rejected %u reading(s), HTTP %d: %s\n", n, code, reply.c_str());
    } else {
      Serial.printf("Upload failed (%d: %s); %u reading(s) waiting\n", code, reply.c_str(), backlog.size());
      return;                                    // keep everything and retry after the next reading
    }
  }
}

void setup() {
  Serial.begin(115200);
  delay(500);
  WiFi.mode(WIFI_STA);
  WiFi.setAutoReconnect(true);
  uint8_t mac[6];
  esp_read_mac(mac, ESP_MAC_WIFI_STA);           // read from eFuse; WiFi.macAddress() is all zeros until WiFi starts
  char macStr[18];
  sprintf(macStr, "%02X:%02X:%02X:%02X:%02X:%02X", mac[0], mac[1], mac[2], mac[3], mac[4], mac[5]);
  Serial.println("MAC address: " + String(macStr));     // UW IT may ask for this

  // WiFi diagnostics: is the network visible on 2.4 GHz, and why does a connection attempt fail?
  int found = WiFi.scanNetworks();
  bool seen = false;
  for (int i = 0; i < found; i++) {
    if (WiFi.SSID(i) == WIFI_SSID) {
      seen = true;
      Serial.printf("Found \"%s\": channel %d, signal %d dBm\n", WIFI_SSID, WiFi.channel(i), WiFi.RSSI(i));
    }
  }
  if (!seen) {
    Serial.printf("\"%s\" not found among %d 2.4 GHz networks:\n", WIFI_SSID, found);
    for (int i = 0; i < found; i++) Serial.printf("  \"%s\" (channel %d, %d dBm)\n", WiFi.SSID(i).c_str(), WiFi.channel(i), WiFi.RSSI(i));
  }
  WiFi.scanDelete();
  WiFi.onEvent([](WiFiEvent_t event, WiFiEventInfo_t info) {
    uint8_t reason = info.wifi_sta_disconnected.reason;
    Serial.printf("WiFi failed, reason %d: %s\n", reason, WiFi.disconnectReasonName((wifi_err_reason_t)reason));
  }, ARDUINO_EVENT_WIFI_STA_DISCONNECTED);

  WiFi.begin(WIFI_SSID, WIFI_PASS);
  configTime(0, 0, "pool.ntp.org", "time.nist.gov");   // UTC; syncs in the background once WiFi is up
  for (int i = 0; i < 60 && WiFi.status() != WL_CONNECTED; i++) delay(500);
  Serial.println("IP address: " + WiFi.localIP().toString());
  for (int i = 0; i < 20 && time(nullptr) < MIN_EPOCH; i++) delay(500);    // give NTP a moment before the first reading
  Serial.println(time(nullptr) >= MIN_EPOCH ? "Clock set from NTP" : "Clock not set yet; readings sent without a timestamp");

  tls.setCACert(GANNET_ROOT_CA);
  Serial.printf("Hub %s uploading to %s\n", HUB_NAME, RECEIVER_URL);
  sensors.setResolution(12);
  lastRead = millis();
  queueReading(readProbes());
  uploadBacklog();
}

void loop() {
  if (millis() - lastRead >= READ_EVERY_MS) {
    lastRead = millis();
    queueReading(readProbes());
    uploadBacklog();
  }

  if (WiFi.status() != WL_CONNECTED) {
    if (!wifiLostAt) wifiLostAt = millis();
    if (millis() - wifiLostAt > WIFI_GIVEUP_MS) ESP.restart();
  } else {
    wifiLostAt = 0;
  }
}
