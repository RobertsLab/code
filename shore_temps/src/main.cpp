#include <Arduino.h>
// shore_temps.ino — SHORE tank temperature hub (ESP32 + DS18B20 probes on one 1-Wire bus)
#include <WiFi.h>
#include <WebServer.h>
#include <ESPmDNS.h>
#include <OneWire.h>
#include <DallasTemperature.h>
#include <esp_mac.h>
#include "secrets.h"                             // WIFI_SSID, WIFI_PASS; copy from secrets.h.example

const char* HUB_NAME  = "shore-temps-a";        // "shore-temps-b" on the second hub
const int   ONE_WIRE_PIN = 4;                    // GPIO 4 (A5 on HUZZAH32)

const unsigned long READ_EVERY_MS   = 60UL * 1000;       // read probes once a minute
const unsigned long WIFI_GIVEUP_MS  = 10UL * 60 * 1000;  // reboot if WiFi is down this long

OneWire oneWire(ONE_WIRE_PIN);
DallasTemperature sensors(&oneWire);
WebServer server(80);

String latestJson = "{}";
unsigned long lastRead = 0, wifiLostAt = 0;

String idToString(const DeviceAddress a) {
  char buf[17];
  for (int i = 0; i < 8; i++) sprintf(buf + 2 * i, "%02X", a[i]);
  return String(buf);
}

void readProbes() {
  sensors.begin();                               // re-scan so a replugged probe reappears
  sensors.requestTemperatures();                 // ~750 ms at 12-bit resolution
  String out = "{\"hub\":\"" + String(HUB_NAME) + "\",\"uptime_s\":" + String(millis() / 1000) +
               ",\"probes\":[";
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
  latestJson = out;
  Serial.println(latestJson);
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
  WiFi.begin(WIFI_SSID, WIFI_PASS);
  for (int i = 0; i < 60 && WiFi.status() != WL_CONNECTED; i++) delay(500);
  Serial.println("IP address: " + WiFi.localIP().toString());
  MDNS.begin(HUB_NAME);                          // http://shore-temps-a.local/ on the same subnet

  sensors.setResolution(12);
  server.on("/", []() { server.send(200, "application/json", latestJson); });
  server.begin();
  readProbes();
  lastRead = millis();
}

void loop() {
  server.handleClient();
  if (millis() - lastRead >= READ_EVERY_MS) { readProbes(); lastRead = millis(); }

  if (WiFi.status() != WL_CONNECTED) {
    if (!wifiLostAt) wifiLostAt = millis();
    if (millis() - wifiLostAt > WIFI_GIVEUP_MS) ESP.restart();
  } else {
    wifiLostAt = 0;
  }
}
