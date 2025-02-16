#pragma once

#include <string>

namespace novac
{
/** CReferenceFitResult is a class to store the evaluated parameters of a
    reference file  from evaluating a spectrum */
class CReferenceFitResult
{
public:
    CReferenceFitResult()
        : m_column(0.0), m_columnError(0.0), m_shift(0.0), m_shiftError(0.0), m_squeeze(0.0), m_squeezeError(0.0), m_specieName("")
    {
    }

    CReferenceFitResult(double column, double columnError, double shift, double shiftError, double squeeze, double squeezeError)
        : m_column(column),
        m_columnError(columnError),
        m_shift(shift),
        m_shiftError(shiftError),
        m_squeeze(squeeze),
        m_squeezeError(squeezeError),
        m_specieName("")
    {
    }

    CReferenceFitResult(const CReferenceFitResult& other)
        : m_column(other.m_column),
        m_columnError(other.m_columnError),
        m_shift(other.m_shift),
        m_shiftError(other.m_shiftError),
        m_squeeze(other.m_squeeze),
        m_squeezeError(other.m_squeezeError),
        m_specieName(other.m_specieName)
    {
    }

    /** The resulting column */
    double m_column = 0.0;

    /** The error in the resulting column */
    double m_columnError = 0.0;

    /** The shift that was applied to the reference in the evaluation */
    double m_shift = 0.0;

    /** The uncertainty in the applied shift */
    double m_shiftError = 0.0;

    /** The squeeze that was applied to the reference in the evaluation */
    double m_squeeze = 0.0;

    /** The uncertainty in the applied squeeze */
    double m_squeezeError = 0.0;

    /** The name of the specie that the reference identifies */
    std::string m_specieName = "";
};

}