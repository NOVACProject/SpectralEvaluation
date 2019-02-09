#include "FitWindow.h"

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
    this->channel       = w2.channel;
    this->fitHigh       = w2.fitHigh;
    this->fitLow        = w2.fitLow;
    this->fitType       = w2.fitType;
    this->shiftSky      = w2.shiftSky;
    this->interlaceStep = w2.interlaceStep;
    this->name          = std::string(w2.name);
    this->nRef          = w2.nRef;
    this->polyOrder     = w2.polyOrder;
    this->UV            = w2.UV;
    this->specLength    = w2.specLength;
    this->startChannel  = w2.startChannel;

    for (int i = 0; i < w2.nRef; ++i)
    {
        this->ref[i] = w2.ref[i];
    }
    this->fraunhoferRef     = w2.fraunhoferRef;
    this->findOptimalShift  = w2.findOptimalShift;
    return *this;
}

void CFitWindow::Clear()
{
    fitHigh         = 460;
    fitLow          = 320;
    channel         = 0;
    specLength      = 2048;
    startChannel    = 0;
    fitType         = FIT_HP_DIV;
    shiftSky        = true;
    interlaceStep   = 1;
    name            = "SO2";
    nRef            = 0;
    polyOrder       = 5;
    UV              = true;
    for (int i = 0; i < MAX_N_REFERENCES; ++i)
    {
        ref[i].m_path = "";
        ref[i].m_specieName = "";
    }

    fraunhoferRef.m_path        = "";
    fraunhoferRef.m_specieName  = "SolarSpec";
    findOptimalShift = false;
}

}