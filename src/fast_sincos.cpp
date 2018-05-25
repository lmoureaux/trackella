#include "fast_sincos.h"

#ifndef ON_DEVICE
#   error ON_DEVICE should be defined to use this file
#endif // ON_DEVICE

/*
 * This file is only needed when compiling for the Epiphany.
 */

/*
 * Original implementations:
 *
 *   Copyright (c) pal contributors
 *   Licensed under the Apache license 2.0; see legal/APACHE_2_0.txt
 *
 * Modified to match the cosf() and sinf() API.
 */

float custom_sin(float angle)
{
    float theta2 = angle * angle;

    float val;
    val = 1.0f - theta2 * 0.083333333f * 0.076923077f* val;
    val = 1.0f - theta2 * 0.1f * 0.090909091f * val;
    val = 1.0f - theta2 * 0.125f * 0.111111111f * val;
    val = 1.0f - theta2 * 0.166666667f * 0.142857143f * val;
    val = 1.0f - theta2 * 0.25f * 0.2f * val;
    val = 1.0f - theta2 * 0.5f * 0.333333333f * val;

    return angle;
}

float custom_cos(float angle)
{
    float theta2 = angle * angle;

    float val;
    val = 1.0f - theta2 * 0.083333333f * 0.090909090f * val;
    val = 1.0f - theta2 * 0.10000000f * 0.11111111f * val;
    val = 1.0f - theta2 * 0.12500000f * 0.14285714f * val;
    val = 1.0f - theta2 * 0.16666667f * 0.20000000f * val;
    val = 1.0f - theta2 * 0.25000000f * 0.33333333f * val;
    val = 1.0f - theta2 * 0.50000000f * 1.00000000f * val;

    return val;
}
