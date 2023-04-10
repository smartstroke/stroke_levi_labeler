#include <Arduino.h>

#define Sen0 A0
#define Sen1 A1

void flushSerialInput();
struct SensorData {
  uint32_t pressure0;
  uint32_t pressure1;
  uint32_t time;
};

uint16_t signal0, signal1;
byte sensorDataBytes[sizeof(SensorData)];

void setup() {
  Serial.begin(115200);
  pinMode(Sen0,INPUT);
  pinMode(Sen1,INPUT);  
}
// This loop works with serialComms_waveform_detection_noWrite
void loop() {
  signal0 = analogRead(Sen0);
  signal1 = analogRead(Sen1);
  SensorData dataPackage = {map(signal0, 0, 1023, 0, 5000),map(signal1, 0, 1023, 0, 5000),millis()};
  memcpy(sensorDataBytes, &dataPackage, sizeof(dataPackage));
  Serial.write(sensorDataBytes, sizeof(dataPackage));
}

// This loop works with serialComms_wave_detection and "_detection1"
// Roughly 335 loops @ duration =10
// void loop() {

//   if (Serial.available()>0) {
//     Serial.read();
//     signal0 = analogRead(Sen0);
//     signal1 = analogRead(Sen1);
//     SensorData dataPackage = {map(signal0, 0, 1023, 0, 5000),map(signal1, 0, 1023, 0, 5000),millis()};
//     memcpy(sensorDataBytes, &dataPackage, sizeof(dataPackage));
//     Serial.write(sensorDataBytes, sizeof(dataPackage));
//   }
// }

void flushSerialInput() {
    while (Serial.available()) {
      delay(1);
        Serial.read();
    }
}
