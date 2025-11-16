#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>

static const size_t N = 100000000; // 1e8

double elapsed_ms(const std::chrono::high_resolution_clock::time_point& start,
                  const std::chrono::high_resolution_clock::time_point& end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

int main() {
    std::cout << "N = " << N << "\n";

    // ----------- vector -----------
    std::vector<double> va(N), vb(N), vd(N);
    for (size_t i = 0; i < N; i++) {
        vb[i] = 1.1;
        vd[i] = 2.2;
    }

    auto start_v = std::chrono::high_resolution_clock::now();
    double c = 3.3;
    for (size_t i = 0; i < N; i++) {
        va[i] = vb[i] * c + vd[i];
    }
    auto end_v = std::chrono::high_resolution_clock::now();
    std::cout << "vector:  " << elapsed_ms(start_v, end_v) << " ms\n";


    // ----------- raw pointer -----------
    double* __restrict pa = (double*)aligned_alloc(64, sizeof(double) * N);
    double* __restrict pb = (double*)aligned_alloc(64, sizeof(double) * N);
    double* __restrict pd = (double*)aligned_alloc(64, sizeof(double) * N);

    for (size_t i = 0; i < N; i++) {
        pb[i] = 1.1;
        pd[i] = 2.2;
    }

    auto start_p = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
        pa[i] = pb[i] * c + pd[i];
    }
    auto end_p = std::chrono::high_resolution_clock::now();
    std::cout << "pointer: " << elapsed_ms(start_p, end_p) << " ms\n";

    free(pa);
    free(pb);
    free(pd);

    return 0;
}
