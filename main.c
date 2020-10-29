#include <stdio.h>
#include "jacobi_curve/jacobi_curve.h"
#include "gost_curve_params.h"
#include "point/point.h"


void print_line() {
    printf("-----------\n\n");
}

void print_is_equals(const point *p1, const point *p2, const jacobian_curve *curve) {
    printf(is_equals(p1, p2, curve) ? "YES\n" : "NO\n");
}

void print_is_on_curve(const point *p, const jacobian_curve *curve) {
    printf("Does the point lie on the curve? ");
    printf(is_on_curve(p, curve) ? "YES\n" : "NO\n");
}

int main() {

    //curve params
    gcry_mpi_t p, q, t, a, x_base, y_base, mpi_k;
    jacobian_curve curve;

    //P
    point point_p;

    gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, p_str, 0, NULL);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, q_str, 0, NULL);
    gcry_mpi_scan(&t, GCRYMPI_FMT_HEX, t_str, 0, NULL);
    gcry_mpi_scan(&a, GCRYMPI_FMT_HEX, a_str, 0, NULL);
    gcry_mpi_scan(&x_base, GCRYMPI_FMT_HEX, x_base_str, 0, NULL);
    gcry_mpi_scan(&y_base, GCRYMPI_FMT_HEX, y_base_str, 0, NULL);


    calculate_jacobi_curve(&curve, p, q, t, a, x_base, y_base);
    create_point(&point_p, curve.x_base, curve.y_base, curve.z_base);

    //TEST1: Q = k*P
    //kP
    point point_k_mult_p;

    //Вычисление кратной точки(539_16 = 1337_10 - random)
    gcry_mpi_scan(&mpi_k, GCRYMPI_FMT_HEX, "539", 0, NULL);
    montgomery(&point_k_mult_p, &point_p, mpi_k, &curve);

    printf("Test 1\n\n");
    printf("Point k*P coords:\n");
    print_point_in_affine(&point_k_mult_p, &curve);
    print_is_on_curve(&point_k_mult_p, &curve);
    printf("\n");

    //TEST2: qP=E
    point point_q_mult_p, point_e;

    create_neutral_point(&point_e);
    montgomery(&point_q_mult_p, &point_p, curve.q, &curve);
    print_line();

    printf("Test 2\n\n");
    printf("Point q*P coors:\n");
    print_point_in_affine(&point_q_mult_p, &curve);
    print_is_on_curve(&point_q_mult_p, &curve);
    printf("\n");

    printf("Are points q*P and E equal? ");
    print_is_equals(&point_q_mult_p, &point_e, &curve);
    printf("\n");
    print_line();

    //TEST3: (q+1)P=P and (q-1)P=-P

    gcry_mpi_t one, q_plus_one, q_minus_one, k1, k2, k1_plus_k2;
    //(q+1)P, qP, -P, (q-1)P
    point point_q_plus_one_mult_p, point_minus_p, point_q_minus_one_mult_p;

    gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, NULL);
    q_plus_one = gcry_mpi_new(0);
    q_minus_one = gcry_mpi_new(0);

    gcry_mpi_addm(q_plus_one, curve.q, one, curve.p);
    gcry_mpi_subm(q_minus_one, curve.q, one, curve.p);

    montgomery(&point_q_plus_one_mult_p, &point_p, q_plus_one, &curve);
    negative_point(&point_minus_p, &point_p);
    montgomery(&point_q_minus_one_mult_p, &point_p, q_minus_one, &curve);

    printf("Test 3\n\n");
    printf("Point (q + 1)*P coords:\n");
    print_point_in_affine(&point_q_plus_one_mult_p, &curve);
    print_is_on_curve(&point_q_plus_one_mult_p, &curve);
    printf("\n");

    printf("Are points (q + 1)*P and P equal? ");
    print_is_equals(&point_q_plus_one_mult_p, &point_p, &curve);
    printf("\n");

    printf("Point (q - 1)*P coords:\n");
    print_point_in_affine(&point_q_minus_one_mult_p, &curve);
    print_is_on_curve(&point_q_minus_one_mult_p, &curve);
    printf("\n");

    printf("Point -P coords:\n");
    print_point_in_affine(&point_minus_p, &curve);
    print_is_on_curve(&point_minus_p, &curve);
    printf("\n");

    printf("Are points (q - 1)*P and -P equal? ");
    print_is_equals(&point_minus_p, &point_minus_p, &curve);
    printf("\n");
    print_line();

    //TEST4: k1 * P + k2 * P = (k1 + k2) * P
    //k1*P, k2*P, (k1+k2)P, k1*P+k2*P
    point point_k1_mult_p, point_k2_mult_p, point_k1_plus_k2_mult_p, point_k1_mult_p_plus_k2_mult_p;

    k1_plus_k2 = gcry_mpi_new(0);
    gcry_mpi_scan(&k1, GCRYMPI_FMT_HEX, "3A57CB", 0, NULL);
    gcry_mpi_scan(&k2, GCRYMPI_FMT_HEX, "332C37", 0, NULL);
    gcry_mpi_addm(k1_plus_k2, k1, k2, curve.p);

    montgomery(&point_k1_mult_p, &point_p, k1, &curve);
    montgomery(&point_k2_mult_p, &point_p, k2, &curve);
    add_point(&point_k1_mult_p_plus_k2_mult_p, &point_k1_mult_p, &point_k2_mult_p, &curve);
    montgomery(&point_k1_plus_k2_mult_p, &point_p, k1_plus_k2, &curve);

    printf("Test 4\n\n");
    printf("Point k1*P + k2*P coors:\n");
    print_point_in_affine(&point_k1_mult_p_plus_k2_mult_p, &curve);
    print_is_on_curve(&point_k1_mult_p_plus_k2_mult_p, &curve);
    printf("\n");

    printf("Point (k1 + k2)*P:\n");
    print_point_in_affine(&point_k1_plus_k2_mult_p, &curve);
    print_is_on_curve(&point_k1_plus_k2_mult_p, &curve);
    printf("\n");

    printf("Are points k1*P + k2*P and (k1 + k2)*P equal? ");
    print_is_equals(&point_k1_mult_p_plus_k2_mult_p, &point_k1_plus_k2_mult_p, &curve);
    printf("\n");
    print_line();


    //RELEASE
    release_point(&point_e);
    release_point(&point_p);
    release_point(&point_q_plus_one_mult_p);
    release_point(&point_q_mult_p);
    release_point(&point_minus_p);
    release_point(&point_q_minus_one_mult_p);
    release_point(&point_k_mult_p);
    release_point(&point_k1_mult_p);
    release_point(&point_k2_mult_p);
    release_point(&point_k1_plus_k2_mult_p);
    release_point(&point_k1_mult_p_plus_k2_mult_p);
    gcry_mpi_release(a);
    gcry_mpi_release(x_base);
    gcry_mpi_release(y_base);
    gcry_mpi_release(one);
    gcry_mpi_release(q_plus_one);
    gcry_mpi_release(q_minus_one);
    gcry_mpi_release(mpi_k);
    gcry_mpi_release(k1);
    gcry_mpi_release(k2);
    gcry_mpi_release(k1_plus_k2);

    return 0;
}
