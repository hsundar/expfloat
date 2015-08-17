#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>

#include <utils.h>
#include <expansion_math.h>
#include <test_apps.h>
#include <iostream>

#define N 100000000
#define CPU_SPEED 2.2E+09

int main() {
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


  printf("\n============================\n");

  free(arr1);
  free(arr2);
  free(darr1);
  free(darr2);


  printf("\n=========== Runge-Kutta ============\n");

  unsigned int n = 100000, r = 9;
  double *qres, *qrhs, *q;
  float *fqres, *fqrhs, *fq;

  qres = new double[n * r];
  qrhs = new double[n * r];
  q = new double[n * r];

  // Initialize arrays ...
  for (int i = 0; i < n * r; ++i) {
    qres[i] = dist(mt);
    qrhs[i] = dist(mt);
    q[i] = dist(mt);
  }

  t1 = rdtsc();
  test_rk4<double, 9>(n, qres, qrhs, q);
  t2 = rdtsc();
  test_rk4<float, 9>(n, fqres, fqrhs, fq);
  t3 = rdtsc();
  test_rk4_exp<9>(n, fqres, fqrhs, fq);
  t4 = rdtsc();

  std::cout << "time for rk4_double: " << (t2 - t1) / CPU_SPEED << "s" << std::endl;
  std::cout << "time for rk4_float:  " << (t3 - t2) / CPU_SPEED << "s" << std::endl;
  std::cout << "time for rk4_exp:    " << (t4 - t3) / CPU_SPEED << "s" << std::endl;

  delete[] q;
  delete[] qrhs;
  delete[] qres;

  printf("\n============================\n");

  return 0;
}

