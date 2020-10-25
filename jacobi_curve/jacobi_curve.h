//
// Created by r.anchugov on 25.10.2020.
//

#ifndef PAZI_LAB_1_JACOBI_CURVE_H
#define PAZI_LAB_1_JACOBI_CURVE_H

#include <gcrypt.h>

typedef struct jacobian_curve jacobian_curve;

/**
 * e,d - Коэфы прямой
 * p - Характеристика простого поля, над которым определяется эллиптическая кривая (модуль)
 * q - Порядок подгруппы простого порядка группы точек эллиптической кривой
 * t - точка второго порядка в форме Вейерштрасса(t,0)(преобразованная)
 * *_base - Координаты точки B на кривой в форме Вейерштрасса(преобразованные)
 */
struct jacobian_curve {
    gcry_mpi_t e, d, p, q, t, x_base, y_base, z_base;
};

void calculate_jacobi_curve(jacobian_curve *curve, gcry_mpi_t p, gcry_mpi_t q, gcry_mpi_t t,
                         gcry_mpi_t a, gcry_mpi_t x_base, gcry_mpi_t y_base);
#endif //PAZI_LAB_1_JACOBI_CURVE_H
