//
// Created by r.anchugov on 26.10.2020.
//

#include "point.h"

void create_neutral_point(point *p) {
    gcry_mpi_t x;
    gcry_mpi_t y;
    gcry_mpi_t z;

    gcry_mpi_scan(&x, GCRYMPI_FMT_HEX, "0", 0, NULL);
    gcry_mpi_scan(&y, GCRYMPI_FMT_HEX, "1", 0, NULL);
    gcry_mpi_scan(&z, GCRYMPI_FMT_HEX, "1", 0, NULL);

    create_point(p, x, y, z);
}

void create_point(point *p, gcry_mpi_t x, gcry_mpi_t y, gcry_mpi_t z) {
    p->x = gcry_mpi_copy(x);
    p->y = gcry_mpi_copy(y);
    p->z = gcry_mpi_copy(z);
}

void release_point(point *p) {
    gcry_mpi_release(p->x);
    gcry_mpi_release(p->y);
    gcry_mpi_release(p->z);
}

void print_point(const point p) {
    unsigned char *buf_x;
    unsigned char *buf_y;
    unsigned char *buf_z;
    gcry_error_t err;

    err = gcry_mpi_aprint(GCRYMPI_FMT_HEX, &buf_x, NULL, p.x);
    if (err != 0) {
        printf("X error: %d\n", err);
        return;
    }

    err = gcry_mpi_aprint(GCRYMPI_FMT_HEX, &buf_y, NULL, p.y);
    if (err != 0) {
        printf("Y error: %d\n", err);
        return;
    }

    err = gcry_mpi_aprint(GCRYMPI_FMT_HEX, &buf_z, NULL, p.z);
    if (err != 0) {
        printf("Z error: %d\n", err);
        return;
    }

    printf("X: %s\nY: %s\nZ: %s\n", buf_x, buf_y, buf_z);
    gcry_free(buf_x);
    gcry_free(buf_y);
    gcry_free(buf_z);
}

void copy(point *p_res, const point *p) {
    create_point(p_res, p->x, p->y, p->z);
}

//Реализация алгоритма сложения двух точке
void add_point(point *p_res, const point *p1, const point *p2, const jacobian_curve *curve) {
    gcry_mpi_t buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9, buf10, buf11, buf12, buf13, buf14, buf15, buf16;
    gcry_mpi_t buf17, buf18, buf19, buf20, buf21, buf22, buf23, buf24, buf25, buf26, buf27, buf28, two;

    buf1 = gcry_mpi_new(0);
    buf2 = gcry_mpi_new(0);
    buf3 = gcry_mpi_new(0);
    buf4 = gcry_mpi_new(0);
    buf5 = gcry_mpi_new(0);
    buf6 = gcry_mpi_new(0);
    buf7 = gcry_mpi_new(0);
    buf8 = gcry_mpi_new(0);
    buf9 = gcry_mpi_new(0);
    buf10 = gcry_mpi_new(0);
    buf11 = gcry_mpi_new(0);
    buf12 = gcry_mpi_new(0);
    buf13 = gcry_mpi_new(0);
    buf14 = gcry_mpi_new(0);
    buf15 = gcry_mpi_new(0);
    buf16 = gcry_mpi_new(0);
    buf17 = gcry_mpi_new(0);
    buf18 = gcry_mpi_new(0);
    buf19 = gcry_mpi_new(0);
    buf20 = gcry_mpi_new(0);
    buf21 = gcry_mpi_new(0);
    buf22 = gcry_mpi_new(0);
    buf23 = gcry_mpi_new(0);
    buf24 = gcry_mpi_new(0);
    buf25 = gcry_mpi_new(0);
    buf26 = gcry_mpi_new(0);
    buf27 = gcry_mpi_new(0);
    buf28 = gcry_mpi_new(0);

    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, NULL);

    //x3 = x1*z1*y2 + y1*x2*z2
    // buf1 = x1*z1
    gcry_mpi_mulm(buf1, p1->x, p1->z, curve->p);
    // buf2 = x1*z1*y2
    gcry_mpi_mulm(buf2, buf1, p2->y, curve->p);
    // buf3 = y1*x2
    gcry_mpi_mulm(buf3, p1->y, p2->x, curve->p);
    // buf4 = y1*x2*z2
    gcry_mpi_mulm(buf4, buf3, p2->z, curve->p);
    // buf5 = x1*z1*y2 + y1*x2*z2
    gcry_mpi_addm(buf5, buf2, buf4, curve->p);

    //y3 = (z1^2*z2^2 + e*x1^2*x2^2)*(y1*y2 - 2*d*x1*x2*z1*z2) + 2*e*x1*x2*z1*z2*(x1^2*z2^2 + z1^2*x2^2)
    // buf6 = z1^2
    gcry_mpi_mulm(buf6, p1->z, p1->z, curve->p);
    // buf7 = z2^2
    gcry_mpi_mulm(buf7, p2->z, p2->z, curve->p);
    // buf8 = z1^2*z2^2
    gcry_mpi_mulm(buf8, buf6, buf7, curve->p);
    // buf9 = x1^2
    gcry_mpi_mulm(buf9, p1->x, p1->x, curve->p);
    // buf10 = x2^2
    gcry_mpi_mulm(buf10, p2->x, p2->x, curve->p);
    // buf11 = x1^2*x2^2
    gcry_mpi_mulm(buf11, buf9, buf10, curve->p);
    // buf12 = e*x1^2*x2^2
    gcry_mpi_mulm(buf12, buf11, curve->e, curve->p);
    // buf13 = z1^2*z2^2 + e*x1^2*x2^2
    gcry_mpi_addm(buf13, buf8, buf12, curve->p);
    // buf14 = y1*y2
    gcry_mpi_mulm(buf14, p1->y, p2->y, curve->p);
    // buf15 = x1*x2
    gcry_mpi_mulm(buf15, p1->x, p2->x, curve->p);
    // buf16 = z1*z2
    gcry_mpi_mulm(buf16, p1->z, p2->z, curve->p);
    // buf17 = x1*x2*z1*z2
    gcry_mpi_mulm(buf17, buf15, buf16, curve->p);
    // buf18 = 2*x1*x2*z1*z2
    gcry_mpi_mulm(buf18, buf17, two, curve->p);
    // buf19 = 2*d*x1*x2*z1*z2
    gcry_mpi_mulm(buf19, buf18, curve->d, curve->p);
    // buf20 = y1*y2 - 2*d*x1*x2*z1*z2
    gcry_mpi_subm(buf20, buf14, buf19, curve->p);
    // buf21 = (z1^2*z2^2 + e*x1^2*x2^2)*(y1*y2 - 2*d*x1*x2*z1*z2)
    gcry_mpi_mulm(buf21, buf13, buf20, curve->p);
    // buf22 = 2*e*x1*x2*z1*z2
    gcry_mpi_mulm(buf22, buf18, curve->e, curve->p);
    // buf23 = x1^2*z2^2
    gcry_mpi_mulm(buf23, buf9, buf7, curve->p);
    // buf24 = z1^2*x2^2
    gcry_mpi_mulm(buf24, buf6, buf10, curve->p);
    // buf25 = x1^2*z2^2 + z1^2*x2^2
    gcry_mpi_addm(buf25, buf23, buf24, curve->p);
    // buf26 = 2*e*x1*x2*z1*z2*(x1^2*z2^2 + z1^2*x2^2)
    gcry_mpi_mulm(buf26, buf22, buf25, curve->p);
    gcry_mpi_addm(buf27, buf21, buf26, curve->p);

    //z3 = z1^2*z2^2 - e*x1^2*x2^2
    // buf28 = z1^2*z2^2 - e*x1^2*x2^2
    gcry_mpi_subm(buf28, buf8, buf12, curve->p);

    p_res->x = buf5;
    p_res->y = buf27;
    p_res->z = buf28;

    gcry_mpi_release(two);
    gcry_mpi_release(buf1);
    gcry_mpi_release(buf2);
    gcry_mpi_release(buf3);
    gcry_mpi_release(buf4);
    gcry_mpi_release(buf7);
    gcry_mpi_release(buf9);
    gcry_mpi_release(buf11);
    gcry_mpi_release(buf12);
    gcry_mpi_release(buf13);
    gcry_mpi_release(buf14);
    gcry_mpi_release(buf15);
    gcry_mpi_release(buf16);
    gcry_mpi_release(buf17);
    gcry_mpi_release(buf18);
    gcry_mpi_release(buf19);
    gcry_mpi_release(buf20);
    gcry_mpi_release(buf21);
    gcry_mpi_release(buf22);
    gcry_mpi_release(buf23);
    gcry_mpi_release(buf24);
    gcry_mpi_release(buf25);
    gcry_mpi_release(buf26);
}

//Реализация "лесенки Монтгомери"
void montgomery(point *p_res, const point *p, gcry_mpi_t scalar, const jacobian_curve *curve) {
    int scalar_bits, i;
    point point_p, point_e;

    copy(&point_p, p);
    create_neutral_point(&point_e);

    scalar_bits = (int) gcry_mpi_get_nbits(scalar);

    for (i = scalar_bits - 1; i >= 0; i--) {
        if (gcry_mpi_test_bit(scalar, i)) {
            add_point(&point_e, &point_e, &point_p, curve);
            add_point(&point_p, &point_p, &point_p, curve);
        } else {
            add_point(&point_p, &point_p, &point_e, curve);
            add_point(&point_e, &point_e, &point_e, curve);
        }
    }
    copy(p_res, &point_e);

    release_point(&point_e);
    release_point(&point_p);
}

// Перевод из проективных в аффинные координаты
void to_affine_coordinates(point *p_res, const point *p, const jacobian_curve *curve) {
    gcry_mpi_t buf1, buf2, buf3, buf4, zero;

    zero = gcry_mpi_new(0);
    buf1 = gcry_mpi_new(0);
    buf2 = gcry_mpi_new(0);
    buf3 = gcry_mpi_new(0);
    buf4 = gcry_mpi_new(0);

    // x = x/z
    // buf1 = 1/z
    gcry_mpi_invm(buf1, p->z, curve->p);
    // buf2 = x/z
    gcry_mpi_mulm(buf2, p->x, buf1, curve->p);

    // y = y/(z^2)
    // buf3 = 1/(z^2)
    gcry_mpi_mulm(buf3, buf1, buf1, curve->p);
    // buf4 = y/(z^2)
    gcry_mpi_mulm(buf4, p->y, buf3, curve->p);

    p_res->x = buf2;
    p_res->y = buf4;
    p_res->z = zero;
}

