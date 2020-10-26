#include <stdio.h>
#include "jacobi_curve/jacobi_curve.h"
#include "gost_curve_params.h"
#include "point/point.h"

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

    point test1, test2;
    create_neutral_point(&test1);
    create_point(&test2, p, q, t);
    print_point(test1);
    print_point(test2);
    point test_res;
    add_point(&test_res, &test1, &test2, &test_curve);
    print_point(test_res);

    gcry_mpi_t two;
    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, NULL);
    montgomery(&test_res, &test2, two, &test_curve);
    print_point(test_res);

    return 0;
}
