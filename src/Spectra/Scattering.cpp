#include <SpectralEvaluation/Spectra/Scattering.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/Fit/Vector.h>

/*******************************************************************************
* DOASIS - DOAS Intelligent System
* ---------------------------------
* Copyright (c) 2004 Stefan Kraus
*                    Institute of Environmental Physics
*                    University of Heidelberg
*
* Any unauthorized use of this software is prohibited!
*
* The application, its components and the source codes are NOT for public use!
* Only non-commerical use is allowed without a written permission of the authors
* or the Institute of Environmental Physics.
* Do not make this source code public available!
*
* A detailed license agreement is given in the
*
* InstallerDOASIS\Setup Files\license.rtf
*
* file. If you want to use this application or its source codes, contact the
* author
*	mailto:stefan.kraus@iup.uni-heidelberg.de
* or the Institute of Environmental Physics
*  mailto:sekretariat@iup.uni-heidelberg.de
*/

using namespace MathFit;

static const double RAD_DEG = (M_PI / 180.0);

// Naturkonstanten (NC)
static const double NC_C = 2.99792458e8;			// Lichtgeschwindigkeit [m/s]
static const double NC_H = 6.6260755e-34;			// Plancksches Wirkungsquantum h [J s]
static const double NC_Kb = 1.380658e-23;			// Boltzmannsche Konstante kB [J/K]

// Konstantendefinitionen
static const int GASE = 2;
static const double DEPOL_RAY = 0.035;     // Depolarisationsfaktor rayleigh
static const double DEPOL_RAM = (6.0 / 7.0);   // Depolarisationsfaktor raman

// spezielle Faktoren
static const double HC_K = ((NC_H * NC_C) / NC_Kb);	// hc/k
static const double CALIBRATE_SIGMA = ((256.0 * M_PI * M_PI * M_PI * M_PI * M_PI) / 27.0 * 1.0e-20 * HC_K);

namespace Doasis
{
CVector Scattering::CalcVRamanSpectrum(CVector& vWavelength, CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA, double v)
{
    if (vWavelength.GetSize() <= 0 || vOrigSpec.GetSize() != vWavelength.GetSize())
    {
        return CVector();
    }

    // Grösse anpassen
    CVector vRamanSpec(vWavelength.GetSize());

    // Statische Datenfelder initialisieren (first is N2, second O2)
    double fNuTilde[] = { 233100, 155500 }; // vibrationsenergieeigenwertabstände [m^-1]
    double fGamma2[] = { 0.518, 1.35 };	// Anisotropie der Polarisierbarkeit [1e-36 m^6]
    double fB0[] = { 198.96, 143.77 };		// Rotationskonstante im niedrigsten Schwingungszustand [m^-1]
    int iGj[2][2] = { {6, 3},{0, 1} };		// Entartungsfaktoren der Grundzustände (gerade/ungerade)
    int iSpin[] = { 1, 0 };				// Kernspin

    // Umkehrfunktion der Dispersionsrelation berechnen
    CVector vIndex(vWavelength.GetSize());
    vIndex.Wedge(0, 1);
    CCubicSplineFunction csfInverseWave;
    csfInverseWave.SetData(vWavelength, vIndex);

    // Hilfsgrößen anlegen
    int iNChannels = vWavelength.GetSize();
    int	iCount = 2 * GASE * iNChannels;				// Gesamtzahl der zu überlagernden Wellenlängen

    CVector vScatWave(iCount);
    CVector vScatInt(iCount);

    // Verhältnis der Phasenfunktionen phi := raman/ramleigh berechnen
    double fTmpC = std::cos(fSZA * RAD_DEG);
    double fTmpC2 = fTmpC * fTmpC;
    double fPhiRam = ((1.0 + DEPOL_RAM + (1.0 - DEPOL_RAM) * fTmpC2) / (1.0 + 0.5 * DEPOL_RAM));
    double fPhiRay = ((1.0 + DEPOL_RAY + (1.0 - DEPOL_RAY) * fTmpC2) / (1.0 + 0.5 * DEPOL_RAY));
    double fPhi = fPhiRam / fPhiRay;

    // Mischungsverhältnis der Gase O2 und N2
    double fMixRatio[GASE];  // Mischungsverhältnis der Gase
    fMixRatio[0] = fMixing;
    fMixRatio[1] = 1 - fMixRatio[0];

    // Zustandssumme berechnen (es fehlt der Faktor k/hc) [m/K]
    double fZ[GASE];
    double fZv[GASE];
    double fSpinFactor[GASE][2];

    for (int i = 0; i < GASE; i++)
    {
        // Zustandssumme
        fZ[i] = (((2.0 * iSpin[i] + 1) * (2.0 * iSpin[i] + 1) * fTemp) / (2.0 * fB0[i])); // JL: ok, k/hc weggelassen
        fZv[i] = std::exp(-HC_K / fTemp * fNuTilde[i]) / (1 - std::exp(-HC_K / fTemp * fNuTilde[i]));	// added by JL

        // Spinfaktoren
        double fTmp = (fMixRatio[i]) / fZ[i] / fZv[i];
        fSpinFactor[i][0] = fTmp * iGj[i][0];					// Entartungsfaktor
        fSpinFactor[i][1] = fTmp * (iGj[i][0] + iGj[i][1]);		// Summe der Entartungsfaktoren
    }

    // Zwischenergebnisse
    double fDeltaNu[3 * GASE];  // Corrected 2023-09-21, used to be 2*GASE 
    double fSigma[3 * GASE];    // Corrected 2023-09-21, used to be 2*GASE 

    // raman-Streuung berechnen (äußere Iteration geht über die Drehimpulse)
    vRamanSpec.Zero();

    for (int j = 0; j < iJMax; j++)
    {
        // Placzek-Teller Koeffizienten incl. Entartungsfaktor (2*j+1)
        /*
        double fStokes = ((9.0 / 8.0) + (3.0 / 4.0) * j - (3.0 / 8.0) / (2 * j + 3));
        double fQ = 0;
        double fAntiStokes = ((-3.0 / 8.0) + (3.0 / 4.0) * j - (3.0 / 8.0) / (2 * j - 1));
        */
        // keine ahnung wie die originaleinträge zustandekommen, deswegen noch einmal, adaptiert von H.Haug 1996
        double fStokes = (3 * (j + 1) * (j + 2) / (2 * (2 * j + 1) * (2 * j + 3))) / (2 * j + 1);
        double fQ = (j * (j + 1) / ((2 * j - 1) * (2 * j + 3))) / (2 * j + 1);
        double fAntiStokes = (3 * (j - 1) / (2 * (2 * j + 1) * (2 * j - 1))) / (2 * j + 1);

        double fTmpE = ((-HC_K) * (j * (j + 1.0)) / fTemp);
        double fTmpEV = ((-HC_K) * (v + 1 / 2) / fTemp);			// vibration
        int iNuIndex = 0;
        int iSigmaIndex = 0;

        // Schleife über die Gase
        for (int g = 0; g < GASE; g++)
        {
            // Verschiebung in der Frequenz [nm^-1]
            fDeltaNu[iNuIndex++] = -fNuTilde[g] + (-fB0[g] * (4.0 * j + 6.0) * 1.0e-9);			// Stokes
            fDeltaNu[iNuIndex++] = -fNuTilde[g];												// Q
            fDeltaNu[iNuIndex++] = -fNuTilde[g] + (fB0[g] * (4.0 * j - 2.0) * 1.0e-9);			// Antistokes

            // Streuquerschnitt
            double fTmpS = (fSpinFactor[g][0] * std::exp(fTmpE * fB0[g] + fTmpEV * fNuTilde[g]));
            fSigma[iSigmaIndex++] = fGamma2[GASE] * fTmpS * fStokes;							// Stokes
            fSigma[iSigmaIndex++] = fGamma2[GASE] * fTmpS * fQ;
            fSigma[iSigmaIndex++] = fGamma2[GASE] * fTmpS * fAntiStokes;						// Antistokes

            // Spinfaktor zwischen gerade und ungerade wechseln
            fSpinFactor[g][0] = fSpinFactor[g][1] - fSpinFactor[g][0];
        }

        // Effekt der raman-Streuung der Quantenzahl j berechnen
        int iIndex = 0;
        for (int l = 0; l < iNChannels; l++)
        {
            iNuIndex = 0;
            iSigmaIndex = 0;
            double fLambda = vWavelength[l];		// Wellenlänge
            double fIntensity = (CALIBRATE_SIGMA * vOrigSpec[l] * fPhi);	// kalibrierte Intensität

            // über die Stokeslinien und die verschiedenen Gase iterieren
            for (int k = (3 * GASE); k > 0; k--)
            {
                // Verschiebung berechnen
                double fLambdaScat = (fLambda / (1.0 + fLambda * fDeltaNu[iNuIndex++]));
                double fLambdaScat2 = fLambdaScat * fLambdaScat;

                // gestreute Wellenlänge merken
                vScatWave[iIndex] = fLambdaScat;
                // gestreute Intensität berechnen
                vScatInt[iIndex] = fIntensity * fSigma[iSigmaIndex++] / (fLambdaScat2 * fLambdaScat2);
                iIndex++;
            }
        }

        // Zielintensitäten berechnen	
        for (int l = 0; l < iCount; l++)
        {
            // Wellenlänge auslesen und in einen Index und den Nachkommaanteil konvertieren
            double fLambda = csfInverseWave.GetValue(vScatWave[l]);
            double fIntensity = vScatInt[l];
            int iIndexLoop = (int)fLambda;

            if ((iIndexLoop >= 0) && (iIndexLoop < iNChannels - 1))
            {
                fLambda -= iIndexLoop;
                // Zielintensität berechnen
                vRamanSpec[iIndexLoop] += fIntensity * (1.0 - fLambda);
                vRamanSpec[iIndexLoop + 1] += fIntensity * fLambda;
            }
        }
    }

    // ok
    return vRamanSpec;
}

CVector Scattering::CalcRamanSpectrum(CVector& vWavelength, CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA)
{
    if (vWavelength.GetSize() <= 0 || vOrigSpec.GetSize() != vWavelength.GetSize())
    {
        return CVector();
    }

    // Grösse anpassen
    CVector vRamanSpec(vWavelength.GetSize());

    // Statische Datenfelder initialisieren
    double fGamma2[] = { 0.518, 1.35 };	// Anisotropie der Polarisierbarkeit [1e-36 m^6]
    double fB0[] = { 198.96, 143.77 };		// Rotationskonstante im niedrigsten Schwingungszustand [m^-1]
    int iGj[2][2] = { {6, 3},{0, 1} };		// Entartungsfaktoren der Grundzustände (gerade/ungerade)
    int iSpin[2] = { 1, 0 };				// Kernspin

    // Umkehrfunktion der Dispersionsrelation berechnen
    CVector vIndex(vWavelength.GetSize());
    vIndex.Wedge(0, 1);
    CCubicSplineFunction csfInverseWave;
    csfInverseWave.SetData(vWavelength, vIndex);

    // Hilfsgrößen anlegen
    int iNChannels = vWavelength.GetSize();
    int	iCount = 2 * GASE * iNChannels;				// Gesamtzahl der zu überlagernden Wellenlängen

    CVector vScatWave(iCount);
    CVector vScatInt(iCount);

    // Verhältnis der Phasenfunktionen phi := raman/ramleigh berechnen
    double fTmpC = std::cos(fSZA * RAD_DEG);
    double fTmpC2 = fTmpC * fTmpC;
    double fPhiRam = ((1.0 + DEPOL_RAM + (1.0 - DEPOL_RAM) * fTmpC2) / (1.0 + 0.5 * DEPOL_RAM));
    double fPhiRay = ((1.0 + DEPOL_RAY + (1.0 - DEPOL_RAY) * fTmpC2) / (1.0 + 0.5 * DEPOL_RAY));
    double fPhi = fPhiRam / fPhiRay;

    // Mischungsverhältnis der Gase O2 und N2
    double fMixRatio[GASE]; // Mischungsverhältnis der Gase
    fMixRatio[0] = fMixing;
    fMixRatio[1] = 1 - fMixRatio[0];

    // Zustandssumme berechnen (es fehlt der Faktor k/hc) [m/K]
    double fZ[GASE];
    double fSpinFactor[GASE][2];

    int i;
    for (i = 0; i < GASE; i++)
    {
        // Zustandssumme
        fZ[i] = (((2.0 * iSpin[i] + 1) * (2.0 * iSpin[i] + 1) * fTemp) / (2.0 * fB0[i]));

        // Spinfaktoren
        double fTmp = (fGamma2[i] * fMixRatio[i]) / fZ[i];
        fSpinFactor[i][0] = fTmp * iGj[i][0];					// Entartungsfaktor
        fSpinFactor[i][1] = fTmp * (iGj[i][0] + iGj[i][1]);		// Summe der Entartungsfaktoren
    }

    // Zwischenergebnisse
    double fDeltaNu[2 * GASE];
    double fSigma[2 * GASE];

    // raman-Streuung berechnen (äußere Iteration geht über die Drehimpulse)
    vRamanSpec.Zero();

    int j;
    for (j = 0; j < iJMax; j++)
    {
        // Placzek-Teller Koeffizienten incl. Entartungsfaktor (2*j+1)
        double fStokes = ((9.0 / 8.0) + (3.0 / 4.0) * j - (3.0 / 8.0) / (2 * j + 3));
        double fAntiStokes = ((-3.0 / 8.0) + (3.0 / 4.0) * j - (3.0 / 8.0) / (2 * j - 1));
        double fTmpE = ((-HC_K) * (j * (j + 1.0)) / fTemp);
        int iNuIndex = 0;
        int iSigmaIndex = 0;

        // Schleife über die Gase
        int g;
        for (g = 0; g < GASE; g++)
        {
            // Verschiebung in der Frequenz [nm^-1]
            fDeltaNu[iNuIndex++] = (-fB0[g] * (4.0 * j + 6.0) * 1.0e-9);		// Stokes
            fDeltaNu[iNuIndex++] = (fB0[g] * (4.0 * j - 2.0) * 1.0e-9);			// Antistokes
            // Streuquerschnitt
            double fTmpS = (fSpinFactor[g][0] * std::exp(fTmpE * fB0[g]));
            fSigma[iSigmaIndex++] = fTmpS * fStokes;							// Stokes
            fSigma[iSigmaIndex++] = fTmpS * fAntiStokes;						// Antistokes
            // Spinfaktor zwischen gerade und ungerade wechseln
            fSpinFactor[g][0] = fSpinFactor[g][1] - fSpinFactor[g][0];
        }

        // Effekt der raman-Streuung der Quantenzahl j berechnen
        int iIndex = 0;
        int l;
        for (l = 0; l < iNChannels; l++)
        {
            iNuIndex = 0;
            iSigmaIndex = 0;
            double fLambda = vWavelength[l];		// Wellenlänge
            double fIntensity = (CALIBRATE_SIGMA * vOrigSpec[l] * fPhi);	// kalibrierte Intensität

            // über die Stokeslinien und die verschiedenen Gase iterieren
            int k;
            for (k = (2 * GASE); k > 0; k--)
            {
                // Verschiebung berechnen
                double fLambdaScat = (fLambda / (1.0 + fLambda * fDeltaNu[iNuIndex++]));
                double fLambdaScat2 = fLambdaScat * fLambdaScat;

                // gestreute Wellenlänge merken
                vScatWave[iIndex] = fLambdaScat;
                // gestreute Intensität berechnen
                vScatInt[iIndex] = fIntensity * fSigma[iSigmaIndex++] / (fLambdaScat2 * fLambdaScat2);
                iIndex++;
            }
        }

        // Zielintensitäten berechnen	
        for (l = 0; l < iCount; l++)
        {
            // Wellenlänge auslesen und in einen Index und den Nachkommaanteil konvertieren
            double fLambda = csfInverseWave.GetValue(vScatWave[l]);
            double fIntensity = vScatInt[l];
            int iIndexLoop = (int)fLambda;

            if ((iIndexLoop >= 0) && (iIndexLoop < iNChannels - 1))
            {
                fLambda -= iIndexLoop;
                // Zielintensität berechnen
                vRamanSpec[iIndexLoop] += fIntensity * (1.0 - fLambda);
                vRamanSpec[iIndexLoop + 1] += fIntensity * fLambda;
            }
        }
    }

    // ok
    return vRamanSpec;
}

CVector Scattering::CalcRingSpectrum(CVector& vWavelength, CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA)
{
    if (vWavelength.GetSize() <= 0 || vWavelength.GetSize() != vOrigSpec.GetSize())
    {
        return CVector();
    }

    // zunächst das Originalspektrum in eine Energie umrechnen (durch lambda teilen)
    CVector vEnergy(vWavelength.GetSize());
    vEnergy.Copy(vOrigSpec);
    vEnergy.DivSimpleSafe(vWavelength);

    // Ramanspektrum berechnen
    CVector vRingSpec = CalcRamanSpectrum(vWavelength, vEnergy, fTemp, iJMax, fMixing, fSZA);

    // Spektren dividieren
    vRingSpec.DivSimpleSafe(vEnergy);

    return vRingSpec;
}

novac::CSpectrum Scattering::CalcRingSpectrum(novac::CSpectrum& specOrig, double fTemp, int iJMax, double fMixing, double fSZA)
{
    // if there is no wavelength information, we can't get the spectrum
    if (!specOrig.IsWavelengthValid())
    {
        return novac::CSpectrum();
    }

    // CVector vRing = CalcRingSpectrum(specOrig.Wavelength.ToVector(), specOrig.Intensity.ToVector(), fTemp, iJMax, fMixing, fSZA);
    CVector wavelength(specOrig.m_wavelength.data(), specOrig.m_length, 1, false);
    CVector intensity(specOrig.m_data, specOrig.m_length, 1, false);
    CVector vRing = CalcRingSpectrum(wavelength, intensity, fTemp, iJMax, fMixing, fSZA);

    novac::CSpectrum specRing;
    specRing.m_length = vRing.GetSize();
    memcpy(specRing.m_data, vRing.GetSafePtr(), vRing.GetSize() * sizeof(double));
    specRing.m_wavelength = std::vector<double>(begin(specOrig.m_wavelength), end(specOrig.m_wavelength));

    return specRing;
}
}

