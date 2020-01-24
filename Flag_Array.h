#pragma once

#include <cassert>
#include <vector>

template <int N>
struct Bitset
{
    static_assert(N == 1 || N == 2 || N == 4 || N == 8);

    Bitset() : size(0), _size(0) {}
    Bitset(int n, int val = 0)
        : size(n), _size(n * N / 8 + 1), arr(std::vector<uint8_t>(_size, val)) {}

    int size;
    int _size;
    std::vector<uint8_t> arr;

    inline void Set(int v, int x)
    {
        int index = v * N / 8;
        assert(v >= 0 && v < size);
        assert(index >= 0 && index < _size);
        constexpr uint8_t mask = ~(UINT8_MAX << N);
        uint8_t shift = (v % (8 / N)) * N;
        arr[index] &= ~(mask << shift);
        arr[index] |= x << shift;
    }

    inline uint8_t Query(int v) const
    {
        int index = v * N / 8;
        assert(v >= 0 && v < size);
        assert(index >= 0 && index < _size);
        constexpr uint8_t mask = ~(UINT8_MAX << N);
        uint8_t shift = (v % (8 / N)) * N;
        return (arr[index] >> shift) & mask;
    }
};

template <>
struct Bitset<1>
{
    static constexpr int N = 1;

    Bitset() : size(0), _size(0) {}
    Bitset(int n, int val = 0) : size(n), _size(n / 8 + 1), arr(std::vector<uint8_t>(_size, val)) {}

    int size;
    int _size;
    std::vector<uint8_t> arr;

    inline void Set(int v)
    {
        int index = v / 8;
        assert(v >= 0 && v < size);
        assert(index >= 0 && index < _size);
        arr[index] |= 1 << (v % 8);
    }

    inline void Unset(int v)
    {
        int index = v / 8;
        assert(v >= 0 && v < size);
        assert(index >= 0 && index < _size);
        arr[index] &= ~(1 << (v % 8));
    }

    inline uint8_t Query(int v) const
    {
        int index = v / 8;
        assert(v >= 0 && v < size);
        assert(index >= 0 && index < _size);
        return arr[index] & (1 << (v % 8));
    }
};

template <int N>
bool BitsetIntersect(int n, Bitset<N>& s1, Bitset<N>& s2)
{
    assert(s1.arr.size() == s2.arr.size() && n * N / 8 + 1 == s1.arr.size());

    for (int i = 0, _end = s1.arr.size(); i < _end; i++)
    {
        if (s1.arr[i] & s2.arr[i])
            return true;
    }

    return false;
}

struct Wide_Flag_Arr
{
    Wide_Flag_Arr(int n, bool val = false) : arr(std::vector<bool>(n, val)) {}

    std::vector<bool> arr;

    inline void Set(int v)
    {
        arr[v] = true;
    }

    inline void Unset(int v)
    {
        arr[v] = false;
    }

    inline bool Query(int v) const
    {
        return arr[v];
    }
};

/*
    If N is small, cache space isn't a problem, reduce instruction count.
    However, if N is large, save cache space by using compact arrays.
 */
#define SMALL_N 0

#if (SMALL_N)
    typedef Wide_Flag_Arr FlagArr;
#else   
    typedef Bitset<1> FlagArr;
#endif