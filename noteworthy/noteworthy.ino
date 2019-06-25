
/*
PlainFFT Library: example of use
Input vectors receive computed results from FFT
No warranty, no claims, just fun
Didier Longueville invenit et fecit October 2010
*/
#include <PlainFFT.h>

int signalPin = 0;
int enablePin = 3;
int greenLight = 8;
int redLight = 9;
int outputSignal = 1;

const uint32_t samples = 256;

double vReal[samples];
double vImag[samples];
//double outHFS[samples];


PlainFFT FFT = PlainFFT(); // Create FFT object
						   // These values can be changed in order to evaluate the functions
						   //double signalFrequency = 1500;
double samplingFrequency = 8928.5;
//uint8_t signalIntensity = 100;
// These are input and output vectors

uint8_t runOnce = 0x00;

int enableCycles = 0;
int signalCycles = 0;

int errorFlag = 0;
double x = 0.0;

#define SCL_INDEX 0x00
#define SCL_TIME 0x01
#define SCL_FREQUENCY 0x02


void setup() {
	Serial.begin(115200);
	Serial.println("Ready");
	while (!Serial) {};

	pinMode(enablePin, INPUT);
	pinMode(signalPin, INPUT);
	pinMode(redLight, OUTPUT);
	pinMode(greenLight, OUTPUT);
	pinMode(outputSignal, OUTPUT);

	//  digitalWrite(greenLight, HIGH);
	digitalWrite(greenLight, LOW);
	digitalWrite(redLight, LOW);

}

void loop() {

	Serial.println("Here?");
	x = 0;
	//while(1);
	// reset every time we go through the loop.
	enableCycles = 0;

	// make sure it is on for at least a certain number of cycles
	// the number will be calibrated through trial and error
	while (enableCycles < 1000) {
		enableCycles = (digitalRead(enablePin) == HIGH) ? enableCycles + 1 : 0;
		Serial.print("enableCycles = "); Serial.println(enableCycles, DEC);// Serial.print('\n');
	}

	digitalWrite(greenLight, HIGH);

	//  int tester;
	//  while(1){
	//    tester = analogRead(signalPin);
	//    Serial.print("signal pin = ");
	//    Serial.println(tester,DEC);
	//    };


	for (int i = 0; i < samples; i++) {
		vReal[i] = analogRead(signalPin);
	}
	printVector(vReal, samples, SCL_TIME);

	//  while(1);

	digitalWrite(greenLight, LOW);
	// PRINT LINES ADDED FOR DEBUGGING

	FFT.windowing(vReal, samples);  // Weigh data
	printVector(vReal, samples, SCL_TIME);

	Serial.println("BeforeFFT");
//	FFT.compute(vReal, vImag, samples, FFT_FORWARD); // Compute FFT
	compute(vReal, vImag, samples, FFT_FORWARD); // Compute FFT
	printVector(vReal, samples, SCL_INDEX);
	printVector(vImag, samples, SCL_INDEX);
	Serial.println("AfterFFT");

	Serial.println("Before");
//	FFT.complexToMagnitude(vReal, vImag, samples); // Compute magnitudes
	complexToMagnitude(vReal, vImag, samples); // Compute magnitudes
	printVector(vReal, (samples >> 1), SCL_FREQUENCY);
	Serial.println("After");

	while (1);

	//for (int i = 0; i < samples; i++) outHFS[i] = vReal[i];
	//scaleItDown(vReal, outHFS, samples, 2);
	//scaleItDown(vReal, outHFS, samples, 3);
	//scaleItDown(vReal, outHFS, samples, 4);
	//scaleItDown(vReal, outHFS, samples, 5);
	//scaleItDown(vReal, outHFS, samples, 6);

	double x = FFT.majorPeak(vReal, samples, samplingFrequency);
	//  Serial.println(x, 6);

	// if the peak is not high enough, send an error message.
	if (errorFlag == 1) {
		digitalWrite(redLight, HIGH);
		//NEED TO prepareDataPacket with an error
	}
	else prepareDataPacket();

	digitalWrite(greenLight, HIGH);
	digitalWrite(redLight, LOW);

}

void printVector(double *vD, uint32_t n, uint8_t scaleType) {
	double timeInterval = (1.0 / samplingFrequency);
	//Serial.print("n = ");Serial.println(n,DEC);

	for (uint16_t i = 0; i < n; i++) {
		// Print abscissa value
		switch (scaleType) {
		case SCL_INDEX:
			Serial.print(i, DEC);
			break;
		case SCL_TIME:
			Serial.print((i * timeInterval), 6);
			break;
		case SCL_FREQUENCY:
			Serial.print((i / (timeInterval * (samples - 1))), 6);
			break;
		}
		Serial.print(" ");
		// Print ordinate value
		Serial.print(vD[i], 6);
		Serial.println();
	}
	Serial.println();
}

void scaleItDown(const double *src, double *dest, int length, int scaleFactor) {
	double factor = 1.0 / scaleFactor;
	int index = -1; // needs to be incremented to zero

	for (int i = 0; i < length; i++) {
		index = (i % scaleFactor == 0) ? i + 1 : i;
		dest[index] += factor * src[i];
	}
}

void prepareDataPacket() {

}

void compute(double *vReal, double *vImag, uint32_t samples, uint8_t dir) {
	// Computes in-place complex-to-complex FFT 
	// Reverse bits
	uint32_t j = 0;
	for (uint32_t i = 0; i < (samples - 1); i++) {
		if (i < j) {
			swap(&vReal[i], &vReal[j]);
			swap(&vImag[i], &vImag[j]);
		}
		uint32_t k = (samples >> 1);
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	// Compute the FFT 
	double c1 = -1.0;
	double c2 = 0.0;
	uint16_t l2 = 1;
	for (uint8_t l = 0; l < exponent(samples); l++) {
		Serial.print("Inside the fft, loop "); Serial.print(l, DEC); Serial.print(" of "); Serial.println(exponent(samples), DEC);
		uint16_t l1 = l2;
		l2 <<= 1;
		double u1 = 1.0;
		double u2 = 0.0;
		for (j = 0; j < l1; j++) {
			for (uint32_t i = j; i < samples; i += l2) {
				uint32_t i1 = i + l1;
				double t1 = u1 * vReal[i1] - u2 * vImag[i1];
				double t2 = u1 * vImag[i1] + u2 * vReal[i1];
				vReal[i1] = vReal[i] - t1;
				vImag[i1] = vImag[i] - t2;
				vReal[i] += t1;
				vImag[i] += t2;
			}
			double z = (u1 * c1) - (u2 * c2);
			u2 = (u1 * c2) + (u2 * c1);
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == FFT_FORWARD) c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}
	Serial.println("Does it get here?");
	// Scaling for forward transform
	if (dir == FFT_FORWARD) {
		for (uint32_t i = 0; i < samples; i++) {
			vReal[i] /= samples;
			vImag[i] /= samples;
		}
	}
}


void swap(double *x, double *y) {
	double temp = *x;
	*x = *y;
	*y = temp;
}

uint8_t exponent(uint32_t value) {
	// computes the exponent of a powered 2 value
	uint8_t result = 0;
	while (((value >> result) & 1) != 1) result++;
	return(result);
}

void complexToMagnitude(double *vReal, double *vImag, uint32_t samples) {
	// vM is half the size of vReal and vImag
	for (uint16_t i = 0; i < samples; i++) vReal[i] = sqrt(sq(vReal[i]) + sq(vImag[i]));
}

