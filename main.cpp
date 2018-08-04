#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <dos.h>
#include <unistd.h>

using namespace std;


const double PI = 3.141592653589793238463;
const double SPEED_OF_LIGHT = 299792458;
const double PLANCK_CONSTANT = 6.626070e-34;

double concentration = 1.38e20 * (pow(100, 3));
double fluorescenceLifeTime = 230e-6;
double absorptionCrossSection = 7.7e-20 / (pow(100, 2));//808
double emissionCrossSection = 2.8e-19 / (pow(100, 2));//1064
double relaxationTime = 200e-12;
double aeLength = 75e-3;
double pumpDiameter = 5e-3;
double lasingModeDiameter = 4.8e-3;
double aeRefractionIndex = 1.819;

double pumpWaveLength = 808e-9;
double lasingWaveLength = 1064e-9;
double singleMatrixPower = 2e3;
unsigned char numberOfMatrix = 3;
double pumpAbsorptionCoefficient = 0.9;

double luminescenceEfficiency = 1e-2;
double amplifiedLasingLuminescenceCoefficient = 1e-10;
double amplifiedNonLasingCoefficient = 1e-5;


double calc(double internalLoss, double closedSwitchInternalLoss, double amplifiedLasingLuminescenceCoefficient);


int main() {
    for (int i = 0; i < 5; i++) {
        usleep(10);
        cout << "=======================  " << i << "  =================\n";
        double internalLoss = i * 0.05;
        double closedSwitchInternalLoss = 0.9;
        double amplifiedLasingLuminescenceCoefficient = 1e-10;
        double p = calc(internalLoss, closedSwitchInternalLoss, amplifiedLasingLuminescenceCoefficient);
        cout << p << "\n";
        usleep(10);
    }
    usleep(10000);
    return 0;
}


double calc(double internalLoss, double closedSwitchInternalLoss, double amplifiedLasingLuminescenceCoefficient) {

    cout << "internalLoss " << internalLoss << "\n";
    cout << "closedSwitchInternalLoss " << closedSwitchInternalLoss << "\n";
    usleep(10);
    const clock_t start = clock();
    double absorbedPumpPower = singleMatrixPower * numberOfMatrix * pumpAbsorptionCoefficient;
    double pumpPowerInLasingMode = absorbedPumpPower * (pow(lasingModeDiameter / pumpDiameter, 2));
    double time = 250e-6;
    double switchTime = 220e-6;
    double pulseSaveStartTime = switchTime;
    double pulseSaveStopTime = switchTime + 200e-9;

    static int stepSkipping = 1000;
    static int pulseStepSkipping = 1;

    static unsigned long stepCount = 1000000;
    double stepSize = time / stepCount;
    auto switchStep = (unsigned long) (switchTime / stepSize);
    auto pulseSaveStartStep = (unsigned long) (pulseSaveStartTime / stepSize);
    auto pulseSaveStopStep = (unsigned long) (pulseSaveStopTime / stepSize);
    cout << "stepSize = " << stepSize * 1e9 << " ns. \n";
    cout << "switchStep = " << switchStep << "\n";


    double geometricResonatorLength = 0.8;
    double equivalentResonatorLength = geometricResonatorLength + (aeRefractionIndex - 1) * aeLength;
    double pumpFrequency = SPEED_OF_LIGHT / pumpWaveLength;
    double lasingFrequency = SPEED_OF_LIGHT / lasingWaveLength;
    double lasingModeVolume = PI * pow(lasingModeDiameter, 2) / 4 * aeLength;

    double pumpRate = pumpPowerInLasingMode / (PLANCK_CONSTANT * pumpFrequency) / lasingModeVolume / concentration;


    double exitMirrorReflection = 0.04;
    double backMirrorReflection = 0.995;

    double logarithmicRoundTripLoss =
            (-log(exitMirrorReflection) + (-log(backMirrorReflection)) + (-log(1 - internalLoss))) / 2;
    cout << "logarithmicRoundTripLoss " << logarithmicRoundTripLoss << "\n";
    double photonInResonatorLifeTime = equivalentResonatorLength / (logarithmicRoundTripLoss * SPEED_OF_LIGHT);
    cout << "photonInResonatorLifeTime " << photonInResonatorLifeTime * 1e9 << " ns.\n";
    double closedSwitchLogarithmicRoundTripLoss =
            (-log(exitMirrorReflection) + (-log(backMirrorReflection)) + (-log(1 - closedSwitchInternalLoss))) / 2;
    cout << "closedSwitchLogarithmicRoundTripLoss " << closedSwitchLogarithmicRoundTripLoss << "\n";
    double closedSwitchPhotonInResonatorLifeTime =
            equivalentResonatorLength / (closedSwitchLogarithmicRoundTripLoss * SPEED_OF_LIGHT);
    cout << "closedSwitchPhotonInResonatorLifeTime " << closedSwitchPhotonInResonatorLifeTime * 1e9 << " ns.\n";

    auto *n1 = (double *) malloc(stepCount * sizeof(double));
    auto *n2 = (double *) malloc(stepCount * sizeof(double));
    auto *photonsInResonator = (double *) malloc(stepCount * sizeof(double));


    n1[0] = concentration;
    n2[0] = 0;

    photonsInResonator[0] = 1;
    double B = emissionCrossSection * aeLength * SPEED_OF_LIGHT / (lasingModeVolume * equivalentResonatorLength);

    double photonInResonatorLifeTimeCurrent;
    photonInResonatorLifeTimeCurrent = closedSwitchPhotonInResonatorLifeTime;
    for (unsigned long step = 1; step <= stepCount; step++) {
        if (step == switchStep) {
            photonInResonatorLifeTimeCurrent = photonInResonatorLifeTime;
        }
        double tmp = B * photonsInResonator[step - 1] * n2[step - 1];
        n2[step] = n2[step - 1] + (pumpRate * n1[step - 1] -
                                   tmp -
                                   n2[step - 1] / fluorescenceLifeTime *
                                   amplifiedLasingLuminescenceCoefficient -
                                   n2[step - 1] / fluorescenceLifeTime) * stepSize;
        n1[step] = concentration - n2[step];
        photonsInResonator[step] = photonsInResonator[step - 1] +
                                   (lasingModeVolume * tmp +
                                    n2[step - 1] / fluorescenceLifeTime *
                                    amplifiedLasingLuminescenceCoefficient -
                                    photonsInResonator[step - 1] / photonInResonatorLifeTimeCurrent) * stepSize;

    }


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
    cout << "pulseDuration " << pulseDuration * 1e9 << " ns.\n";


    double pulseDelay = (maximumPosition - switchStep) * stepSize;
    cout << "pulseDelay " << pulseDelay * 1e9 << " ns.\n";

    cout << "n1 " << n1[switchStep] << "\n";
    cout << "n2 " << n2[switchStep] << "\n";
    cout << "delta " << (n2[switchStep]) / concentration * 100 << "\n";


    clock_t now = clock();
    clock_t delta = now - start;
    double milliSecondsElapsed = static_cast<double>(delta) / CLOCKS_PER_SEC * 1e3;
    cout << "calculation time " << milliSecondsElapsed << " ms\n";


    double photonsToWattCoefficient =
            PLANCK_CONSTANT * lasingFrequency * (1 - exitMirrorReflection) /
            photonInResonatorLifeTime;

    double pulseIntegral = 0;
    for (unsigned long step = pulseSaveStartStep;
         step < (pulseSaveStopStep - pulseStepSkipping);
         step = step + pulseStepSkipping) {
        double pulsePower = photonsInResonator[step] * photonsToWattCoefficient;
        pulseIntegral += pulsePower;
        //  pulsePhotonsFile << step << "\t" << pulsePower << "\n";
    }

    double pulseEnergy = pulseIntegral;
    cout << "energy " << pulseEnergy * stepSize << " J\n";
/*



ofstream n1File;
    ofstream n2File;
    ofstream photonsFile;
    ofstream pulsePhotonsFile;

    n1File.open("n1.txt");
    n2File.open("n2.txt");
    photonsFile.open("photons.txt");
    pulsePhotonsFile.open("pulsePhotonsFile.txt");

    for (unsigned long step = 0;
         step < (stepCount - stepSkipping);
         step = step + stepSkipping) {
        n1File << step << "\t" << n1[step] << "\n";
        n2File << step << "\t" << n2[step] << "\n";
        photonsFile << step << "\t" << photonsInResonator[step] << "\n";
    }

*/


/*

    for (unsigned long step = 0;
         step < (stepCount - stepSkipping);
         step = step + stepSkipping) {
        n1File << step << "\t" << n1[step] << "\n";
        n2File << step << "\t" << n2[step] << "\n";
        photonsFile << step << "\t" << photonsInResonator[step] << "\n";
    }


    for (unsigned long step = pulseSaveStartSKtep;
         step < (pulseSaveStopStep - pulseStepSkipping);
         step = step + pulseStepSkipping
            ) {
        pulsePhotonsFile << step << "\t" << photonsInResonator[step] << "\n";
    }



    n1File.close();
    n2File.close();


    photonsFile.close();
    pulsePhotonsFile.close();


    free(n1);
    free(n2);
    free(photonsInResonator);

    */
    usleep(1000);
    return pulseDuration;
}
