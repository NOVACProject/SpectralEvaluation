#pragma once
#include <SpectralEvaluation/Spectra/Spectrum.h>

#ifndef M_PI
#define M_PI 3.141592
#endif

namespace MathFit
{
    class CVector;
}

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

namespace Doasis
{
    class Scattering
    {
    public:

        Scattering() = default;

        /// <summary>
        /// Calculates the Raman spectrum of the given spectrum
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="vWavelength">The wavelength grid</param>
        /// <param name="vOrigSpec">The spectrum's intensity values</param>
        /// <param name="fTemp">The temperature</param>
        /// <param name="iJMax">The JMAX parameter</param>
        /// <param name="fMixing">The mixing ratio</param>
        /// <param name="fSZA">The SZA</param>
        /// <returns>A vector containg the Raman spectrum's intensity values</returns>
        static MathFit::CVector CalcRamanSpectrum(MathFit::CVector& vWavelength, MathFit::CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA);

        /// <summary>
        /// Calculates the Raman spectrum for a given spectrum.
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Raman spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Raman spectrum</param>
        /// <param name="iJMax">The JMax of the Raman spectrum</param>
        /// <param name="fMixing">The mixing value of the Raman spectrum</param>
        /// <returns>The Raman spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRamanSpectrum(const CSpectrum& specOrig, double fTemp, int iJMax, double fMixing)
        {
            return CalcRamanSpectrum(specOrig, fTemp, iJMax, fMixing, 90);
        }

        /// <summary>
        /// Calculates the Raman spectrum for a given spectrum.
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Raman spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Raman spectrum</param>
        /// <param name="iJMax">The JMax of the Raman spectrum</param>
        /// <returns>The Raman spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRamanSpectrum(const CSpectrum& specOrig, double fTemp, int iJMax)
        {
            return CalcRamanSpectrum(specOrig, fTemp, iJMax, 0.8, 90);
        }

        /// <summary>
        /// Calculates the Raman spectrum for a given spectrum.
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Raman spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Raman spectrum</param>
        /// <returns>The Raman spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * JMax: 30
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRamanSpectrum(const CSpectrum& specOrig, double fTemp)
        {
            return CalcRamanSpectrum(specOrig, fTemp, 30, 0.8, 90);
        }

        /// <summary>
        /// Calculates the Raman spectrum for a given spectrum.
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Raman spectrum is needed</param>
        /// <returns>The Raman spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * Temperature: 250K
        /// * JMax: 30
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRamanSpectrum(const CSpectrum& specOrig)
        {
            return CalcRamanSpectrum(specOrig, 250, 30, 0.8, 90);
        }

        /// <summary>
        /// Calculates the Ring spectrum
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalized spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="vWavelength">The wavelength grid</param>
        /// <param name="vOrigSpec">The intensity values of the spectrum</param>
        /// <param name="fTemp">The temperature</param>
        /// <param name="iJMax">The JMAX parameter</param>
        /// <param name="fMixing">The mixing ratio</param>
        /// <param name="fSZA">The SZA</param>
        /// <returns>A vector that contains the Ring spectrum</returns>
        static MathFit::CVector CalcRingSpectrum(MathFit::CVector& vWavelength, MathFit::CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA);

        /// <summary>
        /// Calculates the Ring spectrum for a given spectrum.
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalied spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Ring spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Ring spectrum</param>
        /// <param name="iJMax">The JMax of the Ring spectrum</param>
        /// <param name="fMixing">The mixing value of the Ring spectrum</param>
        /// <param name="fSZA">The SZA of the Ring spectrum</param>
        /// <returns>The Ring spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        static CSpectrum CalcRingSpectrum(CSpectrum& specOrig, double fTemp, int iJMax, double fMixing, double fSZA);

        /// <summary>
        /// Calculates the Ring spectrum for a given spectrum.
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalied spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Ring spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Ring spectrum</param>
        /// <param name="iJMax">The JMax of the Ring spectrum</param>
        /// <param name="fMixing">The mixing value of the Ring spectrum</param>
        /// <returns>The Ring spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRingSpectrum(CSpectrum& specOrig, double fTemp, int iJMax, double fMixing)
        {
            return CalcRingSpectrum(specOrig, fTemp, iJMax, fMixing, 90);
        }

        /// <summary>
        /// Calculates the Ring spectrum for a given spectrum.
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalied spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Ring spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Ring spectrum</param>
        /// <param name="iJMax">The JMax of the Ring spectrum</param>
        /// <returns>The Ring spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRingSpectrum(CSpectrum& specOrig, double fTemp, int iJMax)
        {
            return CalcRingSpectrum(specOrig, fTemp, iJMax, 0.8, 90);
        }

        /// <summary>
        /// Calculates the Ring spectrum for a given spectrum.
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalied spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Ring spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Ring spectrum</param>
        /// <returns>The Ring spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * JMax: 30
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRingSpectrum(CSpectrum& specOrig, double fTemp)
        {
            return CalcRingSpectrum(specOrig, fTemp, 30, 0.8, 90);
        }

        /// <summary>
        /// Calculates the Ring spectrum for a given spectrum.
        /// </summary>
        /// <remarks>
        /// The Ring spectrum is calculated by dividing the Raman spectrum of the wavelength normalied spectrum by the
        /// wavelength normalized spectrum.
        /// </remarks>
        /// <seealso cref="CalcRamanSpectrum"/>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Ring spectrum is needed</param>
        /// <returns>The Ring spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        /// <remarks>
        /// This method assumes the following preset values:
        /// * Temperature: 250K
        /// * JMax: 30
        /// * Mixing: 0.8
        /// * SZA: 90
        /// </remarks>
        static CSpectrum CalcRingSpectrum(CSpectrum& specOrig)
        {
            return CalcRingSpectrum(specOrig, 250, 30, 0.8, 90);
        }

    private:


        /// <summary>
        /// Calculates the vibrational raman spectrum of the given spectrum. Test function by J. Lampel
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <seealso cref="Haug, H., Vibrationsramanstreuung 1996"/>
        /// <param name="vWavelength">The wavelength grid</param>
        /// <param name="vOrigSpec">The spectrum's intensity values</param>
        /// <param name="fTemp">The temperature</param>
        /// <param name="iJMax">The JMAX parameter</param>
        /// <param name="fMixing">The mixing ratio</param>
        /// <param name="fSZA">The SZA</param>
        /// <returns>A vector containg the Raman spectrum's intensity values</returns>
        static MathFit::CVector CalcVRamanSpectrum(MathFit::CVector& vWavelength, MathFit::CVector& vOrigSpec, double fTemp, int iJMax, double fMixing, double fSZA, double v);

        /// <summary>
        /// Calculates the Raman spectrum for a given spectrum.
        /// </summary>
        /// <seealso cref="Bussemer M.,Der Ring Effekt: Ursachen und Einfluß auf die Spektroskopischen Messungen Stratosphaerischer Spurenstoffe, Thesis, IUP, University of Heidelberg, 1993"/>
        /// <param name="specOrig">The spectrum for which the Raman spectrum is needed</param>
        /// <param name="fTemp">The temperature of the Raman spectrum</param>
        /// <param name="iJMax">The JMax of the Raman spectrum</param>
        /// <param name="fMixing">The mixing value of the Raman spectrum</param>
        /// <param name="fSZA">The SZA of the Raman spectrum</param>
        /// <returns>The Raman spectrum or <b>null</b> if the original spectrum has no wavelength information set.</returns>
        static CSpectrum CalcRamanSpectrum(const CSpectrum& specOrig, double fTemp, int iJMax, double fMixing, double fSZA);

    };
}
