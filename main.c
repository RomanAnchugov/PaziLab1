#include <stdio.h>
#include "jacobi_curve/jacobi_curve.h"
#include " gost_curve_params.h"

void print_mpi(gcry_mpi_t value) {
    unsigned char *buffer;
    gcry_error_t err = gcry_mpi_aprint(GCRYMPI_FMT_HEX, &buffer, NULL, value);
    if (err != 0) {
        printf(" error: %d\n", err);
    }
    printf("MPI %s\n", buffer);
    gcry_free(buffer);
}

int main() {
    jacobian_curve test_curve;
    gcry_mpi_t p, q, t, a, x_base, y_base;

    gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, p_str, 0, NULL);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, q_str, 0, NULL);
    gcry_mpi_scan(&t, GCRYMPI_FMT_HEX, t_str, 0, NULL);
    gcry_mpi_scan(&a, GCRYMPI_FMT_HEX, a_str, 0, NULL);
    gcry_mpi_scan(&x_base, GCRYMPI_FMT_HEX, x_base_str, 0, NULL);
    gcry_mpi_scan(&y_base, GCRYMPI_FMT_HEX, y_base_str, 0, NULL);

    calculate_jacobi_curve(&test_curve, p, q, t, a, x_base, y_base);

    print_mpi(test_curve.e);
    print_mpi(test_curve.d);
    print_mpi(test_curve.x_base);
    print_mpi(test_curve.y_base);

    return 0;
}
