//
// Created by hari on 8/14/15.
//

#ifndef EXPFLOAT_UTILS_H
#define EXPFLOAT_UTILS_H

#include <stdint.h>

inline uint64_t rdtsc() {
    uint32_t lo, hi;
    __asm__ __volatile__ (
    "xorl %%eax, %%eax\n"
            "cpuid\n"
            "rdtsc\n"
    : "=a" (lo), "=d" (hi)
    :
    : "%ebx", "%ecx");
    return (uint64_t)hi << 32 | lo;
}


#endif //EXPFLOAT_UTILS_H
