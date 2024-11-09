#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

#include <iostream>

namespace novac
{

CFitWindow::CFitWindow(const CFitWindow& other)
    : fitLow(other.fitLow),
    fitHigh(other.fitHigh),
    channel(other.channel),
    reference(other.reference),
    fraunhoferRef(other.fraunhoferRef),
    polyOrder(other.polyOrder),
    includeIntensitySpacePolyominal(other.includeIntensitySpacePolyominal),
    ringCalculation(other.ringCalculation),
    specLength(other.specLength),
    startChannel(other.startChannel),
    name(other.name),
    fitType(other.fitType),
    shiftSky(other.shiftSky),
    skyShift(other.skyShift),
    skySqueeze(other.skySqueeze),
    interlaceStep(other.interlaceStep),
    offsetRemovalRange(other.offsetRemovalRange),
    findOptimalShift(other.findOptimalShift),
    child(begin(other.child), end(other.child))
{}

CFitWindow::CFitWindow(CFitWindow&& other)
    : fitLow(other.fitLow),
    fitHigh(other.fitHigh),
    channel(other.channel),
    fraunhoferRef(other.fraunhoferRef),
    polyOrder(other.polyOrder),
    includeIntensitySpacePolyominal(other.includeIntensitySpacePolyominal),
    ringCalculation(other.ringCalculation),
    specLength(other.specLength),
    startChannel(other.startChannel),
    name(other.name),
    fitType(other.fitType),
    shiftSky(other.shiftSky),
    skyShift(other.skyShift),
    skySqueeze(other.skySqueeze),
    interlaceStep(other.interlaceStep),
    offsetRemovalRange(other.offsetRemovalRange),
    findOptimalShift(other.findOptimalShift),
    child(begin(other.child), end(other.child))
{
    reference = std::move(other.reference);
}

CFitWindow& CFitWindow::operator=(const CFitWindow& other)
{
    this->channel = other.channel;
    this->fitHigh = other.fitHigh;
    this->fitLow = other.fitLow;
    this->fitType = other.fitType;
    this->shiftSky = other.shiftSky;
    this->skyShift = other.skyShift;
    this->skySqueeze = other.skySqueeze;
    this->interlaceStep = other.interlaceStep;
    this->name = other.name;
    this->polyOrder = other.polyOrder;
    this->includeIntensitySpacePolyominal = other.includeIntensitySpacePolyominal;
    this->ringCalculation = other.ringCalculation;
    this->offsetRemovalRange = other.offsetRemovalRange;
    this->specLength = other.specLength;
    this->startChannel = other.startChannel;

    this->reference = std::vector<CReferenceFile>(begin(other.reference), end(other.reference));

    this->fraunhoferRef = other.fraunhoferRef;
    this->findOptimalShift = other.findOptimalShift;
    this->child = std::vector<CFitWindow>(begin(other.child), end(other.child));
    return *this;
}

CFitWindow& CFitWindow::operator=(CFitWindow&& other)
{
    this->channel = other.channel;
    this->fitHigh = other.fitHigh;
    this->fitLow = other.fitLow;
    this->fitType = other.fitType;
    this->shiftSky = other.shiftSky;
    this->skyShift = other.skyShift;
    this->skySqueeze = other.skySqueeze;
    this->interlaceStep = other.interlaceStep;
    this->name = std::move(other.name);
    this->polyOrder = other.polyOrder;
    this->includeIntensitySpacePolyominal = other.includeIntensitySpacePolyominal;
    this->ringCalculation = other.ringCalculation;
    this->offsetRemovalRange = other.offsetRemovalRange;
    this->specLength = other.specLength;
    this->startChannel = other.startChannel;

    this->reference = std::move(other.reference);
    this->fraunhoferRef = other.fraunhoferRef;
    this->findOptimalShift = other.findOptimalShift;
    this->child = std::move(other.child);

    return *this;
}

void CFitWindow::Clear()
{
    fitHigh = 460;
    fitLow = 320;
    channel = 0;
    specLength = 2048;
    startChannel = 0;
    fitType = FIT_TYPE::FIT_HP_DIV;
    shiftSky = true;
    interlaceStep = 1;
    name = "SO2";
    polyOrder = 5;
    ringCalculation = RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
    includeIntensitySpacePolyominal = false;
    offsetRemovalRange = IndexRange(50, 200);
    reference.clear();
    reference.reserve(10);

    child.clear();

    fraunhoferRef.m_path = "";
    fraunhoferRef.m_specieName = "SolarSpec";
    findOptimalShift = false;
}

void ReadReferences(CFitWindow& window)
{
    // For each reference in the fit-window, read it in and make sure that it exists...
    for (CReferenceFile& ref : window.reference)
    {
        ref.ReadCrossSectionDataFromFile();
    }

    if (window.fraunhoferRef.m_path.size() > 4)
    {
        window.fraunhoferRef.ReadCrossSectionDataFromFile();
    }

    // If children are defined, then read them as well
    for (CFitWindow& c : window.child)
    {
        ReadReferences(c);
    }
}

void ScaleReferencesToMolecCm2(CFitWindow& window)
{
    // If children are defined, then handle them as well
    for (CFitWindow& c : window.child)
    {
        ScaleReferencesToMolecCm2(c);
    }

    for (CReferenceFile& thisReference : window.reference)
    {
        if (thisReference.m_isFiltered)
        {
            // Convert from ppmm to moleculues / cm2
            Multiply(*thisReference.m_data, (1.0 / 2.5e15));
        }
    }
}

void AddAsReference(CFitWindow& window, const std::vector<double>& referenceData, const std::string& name, int linkShiftToIdx)
{
    CReferenceFile newReference;

    newReference.m_data = std::make_unique<CCrossSectionData>(referenceData);
    newReference.m_specieName = name;
    newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
    newReference.m_columnValue = 1.0;
    newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
    newReference.m_squeezeValue = 1.0;

    if (linkShiftToIdx >= 0)
    {
        newReference.m_shiftOption = SHIFT_TYPE::SHIFT_LINK;
        newReference.m_shiftValue = linkShiftToIdx;
    }
    else
    {
        newReference.m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
        newReference.m_shiftValue = 0.0;
    }

    window.reference.push_back(newReference);
}

size_t AddAsSky(CFitWindow& window, const std::vector<double>& referenceData, SHIFT_TYPE shiftOption)
{
    CReferenceFile newReference;

    newReference.m_data = std::make_unique<novac::CCrossSectionData>(referenceData);
    newReference.m_specieName = "sky";
    newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
    newReference.m_columnValue = -1.0;
    newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
    newReference.m_squeezeValue = 1.0;
    newReference.m_shiftOption = shiftOption;
    newReference.m_shiftValue = 0.0;

    window.reference.push_back(newReference);

    return window.reference.size();
}
}
