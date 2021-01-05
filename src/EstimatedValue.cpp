#include <SpectralEvaluation/EstimatedValue.h>

double EstimatedValue::RelativeError() const
{
    return this->error / this->value;
}
