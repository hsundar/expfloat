//
// Created by hari on 8/14/15.
//

#ifndef EXPFLOAT_EXPANSION_MATH_H
#define EXPFLOAT_EXPANSION_MATH_H

#define S 16

// @brief computes the sum of two single precision floating point values and computes
// their double precision sum in the expansion form

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
two_sum(float a, float b, float& x, float& y) {
    float av, bv, ar, br;

    x = a+b;

    bv = x - a;
    av = x - bv;

    br = b - bv;
    ar = a - av;

    y = ar+ br;
}

inline void
grow_expansion (float& e1, float& e2, float b) {
    float q, h;

    two_sum(b, e2, q, h);
    two_sum(q, e1, e1, e2);
}


inline void
fast_expansion_sum (float& e1, float& e2, float f1, float f2) {
    float q,h;

    two_sum(f2, e2, q, h);
    if (e1 < f1) { // want branch to fail
        two_sum(q,e1,q,h);
        two_sum(q,f1,e1,e2);
    } else {
        two_sum(q,f1,q,h);
        two_sum(q,e1,e1,e2);

    }

}

inline void
split (float a, float& a_hi, float& a_lo) {
    float c, ab;

    c = ((1 << S) + 1) * a;
    ab = c - a;
    a_hi = c - ab;
    a_lo = a - a_hi;
}

inline void
two_product (float a, float b, float&x, float& y) {
    float a_hi, a_lo, b_hi, b_lo, err1, err2;

    x = a*b;

    split (a, a_hi, a_lo);
    split (b, b_hi, b_lo);

    err1 = x - (a_hi * b_hi);
    err2 = err1 - (a_lo * b_hi);
    err1 = err2 - (a_hi * b_lo);

    y = (a_lo * b_lo) - err1;
}

// scale (e1,e2) by a
inline void
scale_expansion(float* e1, float* e2, float a) {
  float q,h,T,t;

  two_product(*e2, a, q, h);
  two_product(*e1, a, T,t);

  two_sum(q,t,q,h);
  two_sum(T,q,*e1,*e2);
}

inline void
daxpy (float *x1, float *x2, float a, float b) {
  scale_expansion(x1, x2, a);
  grow_expansion(*x1, *x2, b);
}

// todo Add code for dot product, matvec, dgemm, rk4

inline void exp_dot (float* a, float* b, unsigned int n, float& r1, float& r2) {
    float x=0.0, y=0.0;
    r1=0.0; r2=0.0;
    for (int i = 0; i < n; ++i) {
        // two_product(a[i], b[i], x, y);
        // fast_expansion_sum(r1, r2, x, y);
        grow_expansion(r1, r2, a[i]*b[i]);
    }
}

template <typename T>
inline T dot (T* a, T* b, int n) {
    T res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += a[i]*b[i];
    }
    return res;
}



#endif //EXPFLOAT_EXPANSION_MATH_H
