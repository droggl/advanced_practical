#include "random_utils.h"

bool randomBoolean() {
    static std::default_random_engine generator(std::random_device{}());
    static std::bernoulli_distribution distribution(0.5);
    return distribution(generator);
}

// Explicit instantiation of the template function for int and double
template int randomFrom(const int min, const int max);
template double randomFrom(const double min, const double max);
