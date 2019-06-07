#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Utils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

#include <iostream>

namespace Evaluation
{

    CFitWindow::CFitWindow()
    {
        Clear();
    }

    CFitWindow::~CFitWindow()
    {
    }

    CFitWindow::CFitWindow(const CFitWindow &wnd)
    {
        *this = wnd;
    }


    CFitWindow &CFitWindow::operator =(const CFitWindow &w2)
    {
        this->channel = w2.channel;
        this->fitHigh = w2.fitHigh;
        this->fitLow = w2.fitLow;
        this->fitType = w2.fitType;
        this->shiftSky = w2.shiftSky;
        this->interlaceStep = w2.interlaceStep;
        this->name = std::string(w2.name);
        this->nRef = w2.nRef;
        this->polyOrder = w2.polyOrder;
        this->UV = w2.UV;
        this->specLength = w2.specLength;
        this->startChannel = w2.startChannel;

        for (int i = 0; i < w2.nRef; ++i)
        {
            this->ref[i] = w2.ref[i];
        }
        this->fraunhoferRef = w2.fraunhoferRef;
        this->findOptimalShift = w2.findOptimalShift;
        return *this;
    }

    void CFitWindow::Clear()
    {
        fitHigh = 460;
        fitLow = 320;
        channel = 0;
        specLength = 2048;
        startChannel = 0;
        fitType = FIT_HP_DIV;
        shiftSky = true;
        interlaceStep = 1;
        name = "SO2";
        nRef = 0;
        polyOrder = 5;
        UV = true;
        for (int i = 0; i < MAX_N_REFERENCES; ++i)
        {
            ref[i].m_path = "";
            ref[i].m_specieName = "";
        }

        fraunhoferRef.m_path = "";
        fraunhoferRef.m_specieName = "SolarSpec";
        findOptimalShift = false;
    }

    bool ReadReferences(CFitWindow& window)
    {
        // For each reference in the fit-window, read it in and make sure that it exists...
        for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
        {
            // Read in the cross section
            if (window.ref[referenceIndex].ReadCrossSectionDataFromFile())
            {
                std::cout << "Failed to read cross section file: " << window.ref[referenceIndex].m_path.c_str() << std::endl;
                return false;
            }
        }

        if (window.fraunhoferRef.m_path.size() > 4)
        {
            if (window.fraunhoferRef.ReadCrossSectionDataFromFile())
            {
                std::cout << "Failed to read Fraunhofer reference file: " << window.fraunhoferRef.m_path.c_str() << std::endl;
                return false;
            }
        }

        return true;
    }

    void HighPassFilterReferences(CFitWindow& window)
    {
        if (window.fitType != Evaluation::FIT_HP_DIV && window.fitType != Evaluation::FIT_HP_SUB)
        {
            return;
        }

        for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
        {
            // Local handle for more convenient syntax.
            CReferenceFile& thisReference = window.ref[referenceIndex];

            if (!thisReference.m_isFiltered)
            {
                if (EqualsIgnoringCase(thisReference.m_specieName, "ring"))
                {
                    HighPassFilter_Ring(*thisReference.m_data);
                }
                else
                {
                    HighPassFilter(*thisReference.m_data);
                }
            }
        }
    }

    void ScaleReferencesToMolecCm2(CFitWindow& window)
    {
        for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
        {
            // Local handle for more convenient syntax.
            CReferenceFile& thisReference = window.ref[referenceIndex];

            if (thisReference.m_isFiltered)
            {
                // Convert from ppmm to moleculues / cm2
                Multiply(*thisReference.m_data, (1.0 / 2.5e15));
            }
        }
    }

}