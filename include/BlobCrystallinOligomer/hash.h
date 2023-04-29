// hash.h

#ifndef HASH_H
#define HASH_H

#include <iostream>

/* Copied from a stack exchange question (which is copied from the BOOST
   library) for allowing pairs to be hashed.
*/
namespace std {

// Originally was outside of std
template <typename T>
void hash_combine(std::size_t& seed, const T& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <typename S, typename T>
struct hash<pair<S, T>> {
    size_t operator()(const pair<S, T>& v) const {
        size_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }
};
} // namespace std

#endif // HASH_H
