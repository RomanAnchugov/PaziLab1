//
// Created by r.anchugov on 25.10.2020.
//

#ifndef PAZI_LAB_1_JACOBI_CURVE_H
#define PAZI_LAB_1_JACOBI_CURVE_H

#include <gcrypt.h>

typedef struct jacobian_curve jacobian_curve;

struct jacobian_curve {
    gcry_mpi_t e, d, p, q, t, x_base, y_base, z_base;
};
#endif //PAZI_LAB_1_JACOBI_CURVE_H
