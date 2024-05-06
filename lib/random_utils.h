#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>

bool randomBoolean();

template <typename T>
T randomFrom(const T min, const T max) {
    static std::random_device rdev;
    static std::default_random_engine re(rdev());
    typedef typename std::conditional<
        std::is_floating_point<T>::value,
        std::uniform_real_distribution<T>,
        std::uniform_int_distribution<T>>::type dist_type;
    dist_type uni(min, max);
    return static_cast<T>(uni(re));
}

#endif // RANDOM_UTILS_H
