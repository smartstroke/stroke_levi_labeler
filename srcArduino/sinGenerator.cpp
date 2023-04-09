#include <Arduino.h>

int waveIndex;
int waveDelta;
int waveOffset;
int waveMax;
int halfWaveMax;
int outputScale;

unsigned long startTime;
double signal_0;
double signal_1;
int mappedSignal_0;
int mappedSignal_1;

int bufferSize;

// the setup function runs once when you press reset
// or power the board
void setup() {
    waveIndex = 0;
    waveDelta = 1;
    waveOffset = (rand() / RAND_MAX) * 2500;
    waveMax = 1500;
    halfWaveMax = waveMax / 2;
    outputScale = 5000;
    startTime = millis();

    Serial.begin(19200);
    bufferSize = Serial.availableForWrite();
}

void updateWaveIndex() {
    waveIndex += waveDelta;
    if (waveIndex >= waveMax) {
        waveIndex = 0;
    }
}

void flushSerialInput() {
    while (Serial.available()) {
        Serial.read();
    }
}

// the loop function runs over and over again forever
void loop() {
    updateWaveIndex();

    // Throttle sample rate for battery life considerations
    delay(1);

    if (!Serial.available()) {
        return;
    }
    flushSerialInput();

    // Prep time data
    unsigned long time = millis() - startTime;

    // Prep first data stream
    double s1Argument = (((double)waveIndex) / ((double)halfWaveMax)) * PI;
    signal_0 = sin(s1Argument) + 1.0;
    mappedSignal_0 = signal_0 * outputScale;

    // Prep second data stream
    double s2Argument = ((double)(waveIndex + waveOffset) / ((double)halfWaveMax)) * PI;
    signal_1 = sin(s2Argument) + 1.0;
    mappedSignal_1 = signal_1 * outputScale;

    // Send data
    Serial.print(time);
    Serial.print(" ");
    Serial.print(mappedSignal_0);
    Serial.print(" ");
    Serial.print(mappedSignal_1);
    Serial.print(" ");
    Serial.print(time / 1000);
    Serial.print(" ");
    Serial.print(time % 1000);
    Serial.println();
}