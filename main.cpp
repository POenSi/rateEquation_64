#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

int main() {

    const double PI = 3.141592653589793238463;
    double time = 250e-6;
    double switchTime = 200e-6;
    double pulseSaveStartTime = switchTime;
    double pulseSaveStopTime = switchTime + 50e-9;

    static int stepSkipping = 100000;
    static int pulseStepSkipping = 1;

    static unsigned long stepCount = 100000000;
    double stepSize = time / stepCount;
    auto switchStep = (unsigned long) (switchTime / stepSize);
    auto pulseSaveStartStep = (unsigned long) (pulseSaveStartTime / stepSize);
    auto pulseSaveStopStep = (unsigned long) (pulseSaveStopTime / stepSize);
    cout << "stepSize = " << stepSize * 1e9 << " ns. \n";

    cout << "switchStep = " << switchStep << "\n";

    double concentration = 1.38e20 * (pow(100, 3));
    double fluorescenceLifeTime = 230e-6;
    double absorptionCrossSection = 7.7e-20 / (pow(100, 2));//808
    double emissionCrossSection = 2.8e-19 / (pow(100, 2));//1064
    double relaxationTime = 200e-12;
    double aeLength = 75e-3;
    double aeDiameter = 5e-3;
    double aeRefractionIndex = 1.819;
    double speedOfLight = 299792458;
    double planckConstant = 6.626070e-34;
    double pumpWaveLength = 808e-9;
    double lasingWaveLength = 1064e-9;
    double singleMatrixPower = 2e3;
    unsigned char numberOfMatrix = 3;
    double pumpAbsorptionCoefficient = 0.9;
    double absorbedPumpPower = singleMatrixPower * numberOfMatrix * pumpAbsorptionCoefficient;
    double luminescenceEfficiency = 1e-2;
    double amplifiedLuminescenceCoefficient = 1e-50;

    double geometricResonatorLength = 0.2;
    double equivalentResonatorLength = geometricResonatorLength + (aeRefractionIndex - 1) * aeLength;
    double pumpFrequency = speedOfLight / pumpWaveLength;
    double lasingFrequency = speedOfLight / lasingWaveLength;
    double aeVolume = PI * pow(aeDiameter, 2) / 4 * aeLength;

    double pumpRate = absorbedPumpPower / (planckConstant * pumpFrequency) / aeVolume / concentration;


    double exitMirrorReflection = 0.1;
    double mirrorReflection = 0.99;
    double internalLoss = 0.01;
    double closedSwitchInternalLoss = 1;
    double logarithmicRoundTripLoss =
            (-log(exitMirrorReflection) + (-log(mirrorReflection))) / 2 + (-log(1 - internalLoss));
    double closedSwitchLogarithmicRoundTripLoss =
            (-log(exitMirrorReflection) + (-log(mirrorReflection))) / 2 + (-log(1 - closedSwitchInternalLoss));
    double photonInResonatorLifeTime = equivalentResonatorLength / (speedOfLight * logarithmicRoundTripLoss);
    double closedSwitchPhotonInResonatorLifeTime =
            equivalentResonatorLength / (speedOfLight * closedSwitchLogarithmicRoundTripLoss);


    auto *__restrict__ n1 = (double *) malloc(stepCount * sizeof(double));
    auto *__restrict__ n2 = (double *) malloc(stepCount * sizeof(double));
    auto *__restrict__ n = (double *) malloc(stepCount * sizeof(double));
    auto *__restrict__ photonsInResonator = (double *) malloc(stepCount * sizeof(double));


    n1[0] = concentration;
    n2[0] = 0;

    photonsInResonator[0] = 1;
    double B = emissionCrossSection * aeLength * speedOfLight / (aeVolume * equivalentResonatorLength);

    double photonInResonatorLifeTimeCurrent;
    photonInResonatorLifeTimeCurrent = closedSwitchPhotonInResonatorLifeTime;
    for (unsigned long step = 1; step < stepCount; step++) {


        if (step == switchStep) {
            cout << "closedPhotonLifeTime " << photonInResonatorLifeTimeCurrent << " \n";
            photonInResonatorLifeTimeCurrent = photonInResonatorLifeTime;
            cout << "openedPhotonLifeTime " << photonInResonatorLifeTimeCurrent << " \n";
        }


        n2[step] = n2[step - 1] + (pumpRate * n1[step - 1] -
                                   B * (photonsInResonator[step - 1] +
                                        n2[step - 1] * fluorescenceLifeTime * amplifiedLuminescenceCoefficient) *
                                   n2[step - 1] -
                                   n2[step - 1] * fluorescenceLifeTime) * stepSize;
        n1[step] = concentration - n2[step];
        photonsInResonator[step] = photonsInResonator[step - 1] +
                                   (aeVolume * B * (photonsInResonator[step - 1] + n2[step - 1] * fluorescenceLifeTime *
                                                                                   amplifiedLuminescenceCoefficient) *
                                    n2[step - 1] -
                                    photonsInResonator[step - 1] / photonInResonatorLifeTimeCurrent) * stepSize;

    }

    ofstream nFile;
    ofstream n1File;
    ofstream n2File;
    ofstream photonsFile;
    ofstream pulsePhotonsFile;
    nFile.open("n.txt");
    n1File.open("n1.txt");
    n2File.open("n2.txt");
    photonsFile.open("photons.txt");
    pulsePhotonsFile.open("pulsePhotonsFile.txt");


    const clock_t start = clock();


    double maximum = 0;
    unsigned long maximumPosition = 0;
    for (unsigned long step = 0; step <= stepCount; step++) {
        if (photonsInResonator[step] > maximum) {
            maximum = photonsInResonator[step];
            maximumPosition = step;
        }
    }

    unsigned long leftHalfHeightPosition = 0;
    unsigned long rightHalfHeightPosition = 0;

    for (unsigned long step = maximumPosition; step > 0; step--) {
        if (photonsInResonator[step] < maximum / 2) {
            leftHalfHeightPosition = step;
            break;
        }
    }

    for (unsigned long step = maximumPosition; step < stepCount; step++) {
        if (photonsInResonator[step] < maximum / 2) {
            rightHalfHeightPosition = step;
            break;
        }
    }


    unsigned long pulseDurationInSteps = rightHalfHeightPosition - leftHalfHeightPosition;

    double pulseDuration = pulseDurationInSteps * stepSize;
    cout << "pulseDuration " << pulseDuration << "\n";
    cout << "maximum " << maximum << " position " << maximumPosition << " \n";


    double pulseDelay = (maximumPosition - switchStep) * stepSize;
    cout << "pulseDelay " << pulseDelay << " \n";


    clock_t now = clock();
    clock_t delta = now - start;
    double milliSecondsElapsed = static_cast<double>(delta) / CLOCKS_PER_SEC * 1e3;
    cout << "calculation time " << milliSecondsElapsed << " ms";


    for (unsigned long step = 0; step < (stepCount - stepSkipping); step = step + stepSkipping) {
        nFile << step << "\t" << n[step] << "\n";
        n1File << step << "\t" << n1[step] << "\n";
        n2File << step << "\t" << n2[step] << "\n";
        photonsFile << step << "\t" << photonsInResonator[step] << "\n";
    }


    for (unsigned long step = pulseSaveStartStep;
         step < (pulseSaveStopStep - pulseStepSkipping); step = step + pulseStepSkipping) {
        pulsePhotonsFile << step << "\t" << photonsInResonator[step] << "\n";
    }


    nFile.close();
    n1File.close();
    n2File.close();
    photonsFile.close();
    pulsePhotonsFile.close();

    free(n);
    free(n1);
    free(n2);
    free(photonsInResonator);

    return 0;
}