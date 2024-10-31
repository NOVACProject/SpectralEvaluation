#include <SpectralEvaluation/EstimatedValue.h>

namespace novac
{

double EstimatedValue::RelativeError() const
{
    return this->error / this->value;
}

}  // namespace novac