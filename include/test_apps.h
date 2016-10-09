//
// Created by hari on 8/17/15.
//

#ifndef EXPFLOAT_TEST_APPS_H
#define EXPFLOAT_TEST_APPS_H

/* low storage 4th order 5 Stage RK scheme
   a and b are NOT the standard RK coefficients */
const double _lsrk45a[5] = {
    0.0,
    -567301805773.0 / 1357537059087.0,
    -2404267990393.0 / 2016746695238.0,
    -3550918686646.0 / 2091501179385.0,
    -1275806237668.0 / 842570457699.0
};

const double        _lsrk45b[5] = {
  1432997174477.0 / 9575080441755.0,
  5161836677717.0 / 13612068292357.0,
  1720146321549.0 / 2090206949498.0,
  3134564353537.0 / 4481467310338.0,
  2277821191437.0 / 14882151754819.0
};

const double        _lsrk45c[5] = {
  0.0,
  1432997174477.0 / 9575080441755.0,
  2526269341429.0 / 6820363962896.0,
  2006345519317.0 / 3224310063776.0,
  2802321613138.0 / 2924317926251.0
};

/* classical 4th order Runge-Kutta scheme */
const double        _crk4a[3] = { .5, .5, 1. };
const double        _crk4b[4] = { 1. / 6., 1. / 3., 1. / 3., 1. / 6. };
const double        _crk4c[4] = { 0., .5, .5, 1. };

/*---------------------------------------------*/

template<typename T, int r>
void test_rk45(unsigned int n, T *qres, T *qrhs, T *q) {
  T rk4a, rk4b, rk4c;
  T dt = 0.01;

  for (int t = 0; t < 10000; ++t) {

    // Low-storage 5-stage 4th-order RK
    for (int k = 0; k < 5; ++k) {
      rk4a = _lsrk45a[k];
      rk4b = _lsrk45b[k];
      rk4c = _lsrk45c[k];

      // for understanding.
      // ---- time_local = time + rk4c * dt;
      // ---- compute RHS ...

      for (int i = 0; i < n; i++) { // spatial unknowns
        for (int j = 0; j < r; ++j) { // vector fields
          qres[i * r + j] = rk4a * qres[i * r + j] + dt * qrhs[i * r + j];
          q[j * n + i] += rk4b * qres[i * r + j];
        }
      }
    }
  }
}

template<int r>
void test_rk45_exp(unsigned int n, float *qres, float *qrhs, float *q) {
  float rk4a, rk4b;
  float dt = 0.01;

  // @todo, move outside the function ...
  float *qtmp = new float[n*r];

  for (int t = 0; t < 10000; ++t) {

    // Low-storage 5-stage 4th-order RK
    for (int k = 0; k < 5; ++k) {
      rk4a = _lsrk45a[k];
      rk4b = _lsrk45b[k];

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < r; ++j) {
          // create extra temp storage of size n*r

          // following sum in exp mode
          daxpy(qres + i*r + j, qtmp + i*r + j, rk4a, dt*qrhs[i*r+j]);

          // use exp for this scaling ?
          q[j * n + i] += rk4b * qres[i * r + j];
        }
      }
    }
  }

  delete [] qtmp;
}

#endif //EXPFLOAT_TEST_APPS_H
