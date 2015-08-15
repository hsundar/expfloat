#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <utils.h>
#include <expansion_math.h>

#define N 100000000
#define CPU_SPEED 2.2E+09

int main() {
    int i;
    float a, b, x, y;
    double da, db, dz;
    float e1 = 0.0, e2 = 0.0, sum = 0.0;
    double dsum = 0.0;
    float *arr = (float *) malloc(N * sizeof(float));
    double *darr = (double *) malloc(N * sizeof(double));

    // unsigned long long t1, t2, t3, t4;
    double t1, t2, t3, t4;

    srand(time(NULL));

    a = rand();
    a /= RAND_MAX;
    b = rand();
    b /= RAND_MAX;

    two_sum(a, b, &x, &y);

    printf("a: %.8f, b: %.8f, x: %.8f, y: %.8e \n", a, b, x, y);

    da = a;
    db = b;

    dz = da + db;

    printf("dz: %.14g, delta: %.9g\n", dz, dz - x);

    printf("\n============================\n");

    for (i = 0; i < N; ++i) {
        da = rand();
        a = da;
        a /= RAND_MAX;
        da /= RAND_MAX;
        arr[i] = a;
        darr[i] = da;
    }

    t1 = rdtsc();
    for (i = 0; i < N; ++i)
        sum = sum + arr[i];
    t2 = rdtsc();
    for (i = 0; i < N; ++i)
        dsum = dsum + darr[i];
    t3 = rdtsc();
    for (i = 0; i < N; ++i)
        grow_expansion(e1, e2, arr[i]);
    t4 = rdtsc();


    printf("sum: %.8f, dsum: %.15g \n", sum, dsum);
    printf("expansion: %.8f, %.8e\n", e1, e2);
    printf("delta: %.9g\n", dsum - e1);

    printf("times: %g, %g, %g\n", (t2 - t1) / CPU_SPEED, (t3 - t2) / CPU_SPEED, (t4 - t3) / CPU_SPEED);

    return 0;
}
