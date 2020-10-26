//
// Created by r.anchugov on 26.10.2020.
//

#ifndef PAZI_LAB_1_POINT_H
#define PAZI_LAB_1_POINT_H

#include <gcrypt.h>
#include "../jacobi_curve/jacobi_curve.h"

typedef struct point point;

struct point {
    gcry_mpi_t x, y, z;
};

void create_neutral_point(point *p);

void create_point(point *p, gcry_mpi_t x, gcry_mpi_t y, gcry_mpi_t z);

void release_point(point *p);

void print_point(point p);

void copy(point *p_res, const point *p);

void add_point(point *p_res, const point *p1, const point *p2, const jacobian_curve *curve);

void montgomery(point *p_res, const point *p, gcry_mpi_t scalar, const jacobian_curve *curve);

void to_affine_coordinates(point *p_res, const point *p, const jacobian_curve *curve);

#endif //PAZI_LAB_1_POINT_H
