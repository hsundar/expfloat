//
// Created by hari on 8/17/15.
//

#ifndef EXPFLOAT_TEST_APPS_H
#define EXPFLOAT_TEST_APPS_H

template<typename T, int r>
void test_rk4(unsigned int n, T *qres, T *qrhs, T *q) {
  T rk4a = 1.0 / 6.0;
  T rk4b = 1.0 / 3.0;
  T dt = 0.01;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < r; ++j) {
      qres[i * r + j] = rk4a * qres[i * r + j] + dt * qrhs[i * r + j];
      q[j * n + i] += rk4b * qres[i * r + j];
    }
  }
}

template<int r>
void test_rk4_exp(unsigned int n, float *qres, float *qrhs, float *q) {

}

#endif //EXPFLOAT_TEST_APPS_H
