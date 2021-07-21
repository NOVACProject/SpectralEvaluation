#include <SpectralEvaluation/Math/PolynomialFit.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/Fit/PolynomialFunction.h>
#include <numeric>

namespace novac
{

// ------------------------- Polynomial function fitting -------------------------

class CVectorDataFunction : public MathFit::IParamFunction
{
private:
    MathFit::CVector xVector;
    MathFit::CVector yVector;
public:
    CVectorDataFunction(const MathFit::CVector& x, const MathFit::CVector& y)
        : xVector(x), yVector(y)
    {
        if (xVector.GetSize() != yVector.GetSize())
        {
            throw std::invalid_argument("Cannot setup a vector data function with not equal length of the x- and y- vectors");
        }
    }

    virtual MathFit::TFitData GetLinearBasisFunction(MathFit::TFitData /*fXValue*/, int /*iParamID*/, bool /*bFixedID*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetLinearBasisFunction is not implemented");
    }

    virtual MathFit::TFitData GetValue(MathFit::TFitData /*fXValue*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetValue is not implemented");
    }

    virtual MathFit::CVector& GetValues(MathFit::CVector& vXValues, MathFit::CVector& vYTargetVector) override
    {
        if (vXValues.GetSize() != this->xVector.GetSize())
        {
            throw std::invalid_argument("The CVectorDataFunction is constructed to only use vectors of equal length.");
        }

        const int size = this->xVector.GetSize();
        for (int i = 0; i < size; i++)
        {
            vYTargetVector.SetAt(i, this->yVector.GetAt(i));
        }

        return vYTargetVector;
    }

    virtual MathFit::TFitData GetSlope(MathFit::TFitData /*fXValue*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetSlope is not implemented");
    }
};

// Copied from https://en.wikipedia.org/wiki/LU_decomposition
/* INPUT: A - array of pointers to rows of a square matrix having dimension N
    *        Tol - small tolerance number to detect failure when the matrix is near degenerate
    * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
    *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
    *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
    *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
    */
int LUPDecompose(double** A, int N, double Tol, int* P) {

    int i, j, k, imax;
    double maxA, * ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done 
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
    * OUTPUT: x - solution vector of A*x=b
    */
void LUPSolve(double** A, int* P, double* b, int N, double* x) {

    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

PolynomialFit::PolynomialFit(int order)
    : polynomialOrder(order)
{
    this->functionToFit = new MathFit::CPolynomialFunction(order);
}

PolynomialFit::~PolynomialFit()
{
    delete this->functionToFit;
}

bool PolynomialFit::FitCubicPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients)
{
    polynomialCoefficients.resize(4);

    double A[16]{}; // the Vandermode Matrix
    for (size_t rowIdx = 0; rowIdx < 4; ++rowIdx)
    {
        A[rowIdx * 4] = 1.0;
        A[rowIdx * 4 + 1] = xData[rowIdx];
        A[rowIdx * 4 + 2] = std::pow(xData[rowIdx], 2.0);
        A[rowIdx * 4 + 3] = std::pow(xData[rowIdx], 3.0);
    }
    double B[4]{ yData[0],  yData[1],  yData[2],  yData[3] };
    // Multiply by transpose (TODO: this can be done manually in the setup above)
    double AtA1[4]{};
    double AtA2[4]{};
    double AtA3[4]{};
    double AtA4[4]{};
    double* AtA[4]{ AtA1, AtA2, AtA3, AtA4 }; // Transpose(A) * A
    double Atb[4]{}; // Transpose(A) * B
    for (size_t rowIdx = 0; rowIdx < 4; ++rowIdx)
    {
        Atb[rowIdx] =
            A[rowIdx + 0] * yData[0] +
            A[rowIdx + 4] * yData[1] +
            A[rowIdx + 8] * yData[2] +
            A[rowIdx + 12] * yData[3];

        for (size_t colIdx = 0; colIdx < 4; ++colIdx)
        {
            AtA[rowIdx][colIdx] =
                A[rowIdx + 0] * A[colIdx + 0] +
                A[rowIdx + 4] * A[colIdx + 4] +
                A[rowIdx + 8] * A[colIdx + 8] +
                A[rowIdx + 12] * A[colIdx + 12];
        }
    }

    // Solve the normal equations AtA * x = AtB
    int permutations[8]{};
    LUPDecompose(AtA, 4, 1e-9, permutations);
    LUPSolve(AtA, permutations, Atb, 4, polynomialCoefficients.data());

    return true;
}

bool PolynomialFit::FitPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients)
{
    if (xData.size() != yData.size())
    {
        return false;
    }
    if (xData.size() == 4 && polynomialOrder == 3)
    {
        return FitCubicPolynomial(xData, yData, polynomialCoefficients);
    }

    const bool autoReleaseData = false;
    MathFit::CVector xVector{ &xData[0], (int)xData.size(), 1, autoReleaseData };
    MathFit::CVector yVector{ &yData[0], (int)yData.size(), 1, autoReleaseData };

    try
    {
        // This is a home-made replacement for the spline. In order to speed up the fitting somewhat.
        CVectorDataFunction cubicSplineRepresentation{ xVector, yVector };

        MathFit::CStandardMetricFunction diff(cubicSplineRepresentation, *(this->functionToFit));

        MathFit::CLeastSquareFit fit(diff);
        fit.SetFitRange(xVector);
        fit.SetMaxFitSteps(500);
        fit.PrepareMinimize();

        (void)fit.Minimize();
        fit.FinishMinimize();

        // Now extract the coefficients
        {
            const auto& modelVector = this->functionToFit->GetCoefficients();
            polynomialCoefficients.resize(modelVector.GetSize());
            for (int orderIdx = 0; orderIdx < modelVector.GetSize(); ++orderIdx)
            {
                polynomialCoefficients[orderIdx] = modelVector.GetAt(orderIdx);
            }
        }

        return true;
    }
    catch (MathFit::CFitException& e)
    {
        std::cout << "Fit failed: " << e.mMessage << std::endl;

        return false;
    }
}

}
