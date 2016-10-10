#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>

#include <utils.h>
#include <expansion_math.h>
#include <test_apps.h>
#include <iostream>
#include <algorithm>
#include <limits>

#include <quadmath.h>
#include <math.h>

#include <inttypes.h>


#define N 1000000
#define CPU_SPEED 2.2E+09

typedef unsigned __int128 uint128_t;
typedef unsigned long uint64;


/*      UINT64_MAX 18446744073709551615ULL */
#define P10_UINT64 10000000000000000000ULL   /* 19 zeroes */
#define E10_UINT64 19

#define STRINGIZER(x)   # x
#define TO_STRING(x)    STRINGIZER(x)

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}


static int print_u128_u(uint128_t u128)
{
    int rc;
    if (u128 > UINT64_MAX)
    {
        uint128_t leading  = u128 / P10_UINT64;
        uint64_t  trailing = u128 % P10_UINT64;
        rc = print_u128_u(leading);
        rc += printf("%." TO_STRING(E10_UINT64) PRIu64, trailing);
    }
    else
    {
        uint64_t u64 = u128;
        rc = printf("%" PRIu64, u64);
    }
    return rc;
}

int main(int argc, char* argv[]) {
  int i;
  float a, b, x, y;
  double da, db, dz;
  float e1 = 0.0, e2 = 0.0, sum = 0.0;
  double dsum = 0.0;
  float *arr1 = (float *) malloc(N * sizeof(float));
  float *arr2 = (float *) malloc(N * sizeof(float));
  double *darr1 = (double *) malloc(N * sizeof(double));
  double *darr2 = (double *) malloc(N * sizeof(double));

  printf("\n=========== a+b ===============\n");

  // unsigned long long t1, t2, t3, t4;
  double t1, t2, t3, t4;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0, 1);

  a = dist(mt);
  b = dist(mt);

  fast_two_sum(a, b, x, y);

  printf("a: %.8f, b: %.8f, x: %.8f, y: %.8e \n", a, b, x, y);

  da = a;
  db = b;

  dz = da + db;

  printf("dz: %.14g, delta: %.9g\n", dz, dz - x);

  printf("\n=========== sum(a) ==============\n");

  for (i = 0; i < N; ++i) {
    da = dist(mt);
    a = da;
    arr1[i] = a;
    darr1[i] = da;

    da = dist(mt);
    a = da;
    arr2[i] = a;
    darr2[i] = da;
  }

  t1 = rdtsc();
  for (i = 0; i < N; ++i)
    sum = sum + arr1[i];
  t2 = rdtsc();
  for (i = 0; i < N; ++i)
    dsum = dsum + darr1[i];
  t3 = rdtsc();
  for (i = 0; i < N; ++i)
    grow_expansion(e1, e2, arr1[i]);
  t4 = rdtsc();


  printf("sum: %.8f, dsum: %.15g \n", sum, dsum);
  printf("expansion: %.8f, %.8e\n", e1, e2);
  printf("delta: %.9g\n", dsum - e1);

  printf("times: %g, %g, %g\n", (t2 - t1) / CPU_SPEED, (t3 - t2) / CPU_SPEED, (t4 - t3) / CPU_SPEED);

  printf("\n=========== a*b ===============\n");

  a = dist(mt);
  b = dist(mt);

  two_product(a, b, x, y);

  printf("a: %.8f, b: %.8f, x: %.8f, y: %.8e \n", a, b, x, y);

  da = a;
  db = b;

  dz = da * db;

  printf("dz: %.14g, delta: %.9g\n", dz, dz - x);

  printf("\n========== (a,b) ==============\n");

  t1 = rdtsc();
  sum = dot(arr1, arr2, N);
  t2 = rdtsc();
  dsum = dot(darr1, darr2, N);
  t3 = rdtsc();
  exp_dot(arr1, arr2, N, e1, e2);
  t4 = rdtsc();


  printf("sum: %.8f, dsum: %.15g \n", sum, dsum);
  printf("expansion: %.8f, %.8e\n", e1, e2);
  printf("delta: %.9g\n", dsum - e1);

  printf("times: %g, %g, %g\n", (t2 - t1) / CPU_SPEED, (t3 - t2) / CPU_SPEED, (t4 - t3) / CPU_SPEED);


  free(arr1);
  free(arr2);
  free(darr1);
  free(darr2);

  /*
  printf("\n=========== Runge-Kutta ============\n");

  unsigned int n = 10000, r = 9;
  double *qres, *qrhs, *q;
  float *fqres, *fqrhs, *fq;

  qres = new double[n * r];
  qrhs = new double[n * r];
  q = new double[n * r];

  fqres = new float[n * r];
  fqrhs = new float[n * r];
  fq = new float[n * r];

  // Initialize arrays ...
  for (int i = 0; i < n * r; ++i) {
    fqres[i] = qres[i] = dist(mt);
    fqrhs[i] = qrhs[i] = dist(mt);
    fq[i] = q[i] = dist(mt);
  }

  t1 = rdtsc();
  test_rk45<double, 1>(n, qres, qrhs, q);
  t2 = rdtsc();
  test_rk45<float, 1>(n, fqres, fqrhs, fq);
  t3 = rdtsc();
  test_rk45_exp<1>(n, fqres, fqrhs, fq);
  t4 = rdtsc();

  std::cout << "time for rk4_double: " << (t2 - t1) / CPU_SPEED << "s" << std::endl;
  std::cout << "time for rk4_float:  " << (t3 - t2) / CPU_SPEED << "s" << std::endl;
  std::cout << "time for rk4_exp:    " << (t4 - t3) / CPU_SPEED << "s" << std::endl;



  delete[] q;
  delete[] qrhs;
  delete[] qres;
  delete[] fq;
  delete[] fqrhs;
  delete[] fqres;

  printf("\n============================\n");
*/


  //! @hari - 8 Oct 2016 - for Fall NSF Large 2016.
  printf("\n=========== H-float ============\n");

/*
  float16_t    half_a;	// 16-bit half-precision binary format
  float32_t    sing_a;  // 32-bit single-precision binary format
  float64_t    doub_a;	// 64-bit double-precision binary format
  extFloat80_t extf_a;  // 80-bit double-extended-precision binary format (old Intel or Motorola format)
  float128_t 	 quad_a;  // 128-bit quadruple-precision binary format
*/

  unsigned long n = atol(argv[1]); //  1000000;
  __float128  *sin2pi_q = new __float128[n];

  char buf[128];
  int width = 46;   

  // sin2pi_q[0] = sinq(0.0q);
  for (int i=0; i<n; ++i) {
    sin2pi_q[i] = sinq(M_PIq * 2 * i / n / atoi(argv[2]) );
    // quadmath_snprintf (buf, sizeof(buf), "%+-.80Qe", hello[i]);
    // printf("%s\n", buf);
  }

  // 128 bits hex - 0x ff ff ff ff ff ff ff ff 
  // loop through ... 
  unsigned long      cnts[12] = {1};
  unsigned __int128  mask[12]; /* = {0x000000000000000000000000000003FF, //  10
                                     0x000000000000000000000000000FFC00, //  20
                                     0x0000000000000000000000003FF00000, //  30
                                     0x0000000000000000000000FFC0000000, //  40
                                     0x00000000000000000003FF0000000000, //  50
                                     0x00000000000000000FFC000000000000, //  60
                                     0x000000000000003FF000000000000000, //  70
                                     0x000000000000FFC00000000000000000, //  80
                                     0x0000000003FF00000000000000000000, //  90
                                     0x0000000FFC0000000000000000000000, // 100
                                     0x00003FF0000000000000000000000000, // 110
                                     0x0000C000000000000000000000000000, // 112
                                    }; */

  unsigned __int128  bits;
  unsigned __int128  base[12];

  bits = *((__int128 *)(sin2pi_q));
  printf("Masks ...\n");
  for (int j=0; j<11; ++j) {
    mask[j] = 0x03FFu;
    mask[j] = mask[j] << j*10; 
    base[j] = bits & mask[j];

    // printf ("%d: mask = ", j);
    // print_u128_u(mask[j]);
    // printf("\n");
    printf("%d:   \t  %016" PRIx64 "%016" PRIx64 "\n",j,(uint64)(mask[j]>>64),(uint64)mask[j]);

  }
  mask[11] = 0x03u;
  mask[11] = mask[11] << 110; 
  base[11] = bits & mask[11];
  printf("%d: \t  %016" PRIx64 "%016" PRIx64 "\n",11,(uint64)(mask[11]>>64),(uint64)mask[11]);

  quadmath_snprintf (buf, sizeof(buf), "%+-.32Qe", sin2pi_q[1] - sin2pi_q[0]);
  printf("sin(2 pi x) -- %s\n", buf);
  quadmath_snprintf (buf, sizeof(buf), "%+-.32Qe", sin2pi_q[2] - sin2pi_q[1]);
  printf("sin(2 pi x) -- %s\n", buf);
  quadmath_snprintf (buf, sizeof(buf), "%+-.32Qe", sin2pi_q[3] - sin2pi_q[2]);
  printf("sin(2 pi x) -- %s\n", buf);
  quadmath_snprintf (buf, sizeof(buf), "%+-.32Qe", sin2pi_q[4] - sin2pi_q[3]);
  printf("sin(2 pi x) -- %s\n", buf);
  quadmath_snprintf (buf, sizeof(buf), "%+-.32Qe", sin2pi_q[5] - sin2pi_q[4]);
  printf("sin(2 pi x) -- %s\n", buf);
 
  __float128 r1,r2,r3,r4;
  double d, d_prev;
  
  unsigned long  d_cnt, f_cnt;


  // approx by 2x double 
  d = static_cast<double> ( sin2pi_q[0] );
  r1 = d;
  r2 = sin2pi_q[0] - r1;  
  d_cnt = 2;
  d_prev = d;
  for (int i=1; i<n; ++i) {
    d = static_cast<double> ( sin2pi_q[i] - sin2pi_q[i-1]);
    double err = std::abs(d-d_prev);
    if ( err < 1e-16 ) { // almost_equal(d, d_prev, 10) ) {
      r2 = sin2pi_q[i]  - sin2pi_q[i-1] - r1;
      d_cnt++;  
    } else {
      // double err = std::abs(d-d_prev);
      // printf("error is %g\n", err);
      r1 = d;
      r2 = sin2pi_q[i]  - sin2pi_q[i-1] - r1;
      d_prev = d;
      d_cnt+=2;  
    }	
  }

  printf("Storage by doubles instead of quads.\n");
  printf("quad storage:\t %ld bytes \n", n*16); 
  printf("double storage:\t %ld bytes \n", d_cnt*8); 

  // approx by 4x float
  

  delete [] sin2pi_q;

  printf("\n============================\n");

	
  printf("\n========double by float=========\n");

  double  *sin2pi = new double[n];

  for (int i=0; i<n; ++i) {
    sin2pi[i] = sin(M_PI * 2 * i / n  );
  }
  
  double d1, d2;
  float f,f_prev;
  
  f = static_cast<float> ( sin2pi[0] );
  d1 = f;
  d2 = sin2pi[0] - d1;  
  f_cnt = 2;
  f_prev = f;
  for (int i=1; i<n; ++i) {
    f = static_cast<float> ( sin2pi[i] - sin2pi[i-1]);
    double err = std::abs(f-f_prev);
    if ( err < 1e-8 ) { // almost_equal(d, d_prev, 10) ) {
      d2 = sin2pi[i] - d1;
      f_cnt++;  
    } else {
      // double err = std::abs(d-d_prev);
      // printf("error is %g\n", err);
      d1 = f;
      d2 = sin2pi[i] - d1;
      f_prev = f;
      f_cnt+=2;  
    }	
  }
  printf("Storage by floats instead of doubles.\n");
  printf("double storage:\t %ld bytes \n", n*8); 
  printf("float storage:\t %ld bytes \n", f_cnt*4); 

  delete [] sin2pi;
  printf("\n============================\n");
  return 0;
}

