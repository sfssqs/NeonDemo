//
// Created by Xia,Xing on 18/4/27.
//

#ifndef HELLO_NEON_YUV_UTILS_H
#define HELLO_NEON_YUV_UTILS_H

#include <stdint.h>

void bgr888_to_yuv444_float(unsigned char *yuv, unsigned char *bgr, int pixel_num);

void bgr888_to_yuv444_c(unsigned char *yuv, unsigned char *bgr, int pixel_num);

void bgr888_to_yuv444_neon(unsigned char *__restrict__ yuv, unsigned char *__restrict__ bgr,
                           int pixel_num);

void reference_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int n);

void neon_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int n);

void neon_asm_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int numPixels);

#endif //HELLO_NEON_YUV_UTILS_H
