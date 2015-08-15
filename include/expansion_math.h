//
// Created by hari on 8/14/15.
//

#ifndef EXPFLOAT_EXPANSION_MATH_H
#define EXPFLOAT_EXPANSION_MATH_H

inline void
fast_two_sum(float a, float b, float& x, float& y) {
    float v;

    x = a+b;
    if ((a>b) == (a>-b)) {
        v  = x - a;
        y = b - v;
    } else {
        v  = x - b;
        y = a - v;
    }
}

inline void
two_sum(float a, float b, float* x, float* y) {
    float av, bv, ar, br;

    *x = a+b;

    bv = *x - a;
    av = *x - bv;

    br = b - bv;
    ar = a - av;

    *y = ar+ br;
}

inline void
grow_expansion (float& e1, float& e2, float b) {
    float q, h;

    fast_two_sum(b, e2, q, h);
    fast_two_sum(q, e1, e1, e2);
}

#endif //EXPFLOAT_EXPANSION_MATH_H
