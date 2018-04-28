#include "yuv-utils.h"

#if defined(HAVE_NEON) && defined(HAVE_NEON_X86)
/*
 * The latest version and instruction for NEON_2_SSE.h is at:
 *    https://github.com/intel/ARM_NEON_2_x86_SSE
 */
#include "NEON_2_SSE.h"
#elif defined(HAVE_NEON)

#include <arm_neon.h>

#endif

void bgr888_to_yuv444_neon(unsigned char *__restrict__ yuv,
                           unsigned char *__restrict__ bgr,
                           int pixel_num) {

    const uint8x8_t u8_zero = vdup_n_u8(0);
    const uint16x8_t u16_rounding = vdupq_n_u16(128);
    const int16x8_t s16_rounding = vdupq_n_s16(128);
    const int8x16_t s8_rounding = vdupq_n_s8(128);

    int count = pixel_num / 16;

    int i;
    for (i = 0; i < count; ++i) {
        // Load bgr
        uint8x16x3_t pixel_bgr = vld3q_u8(bgr);

        uint8x8_t high_r = vget_high_u8(pixel_bgr.val[0]);
        uint8x8_t low_r = vget_low_u8(pixel_bgr.val[0]);
        uint8x8_t high_g = vget_high_u8(pixel_bgr.val[1]);
        uint8x8_t low_g = vget_low_u8(pixel_bgr.val[1]);
        uint8x8_t high_b = vget_high_u8(pixel_bgr.val[2]);
        uint8x8_t low_b = vget_low_u8(pixel_bgr.val[2]);
        int16x8_t signed_high_r = vreinterpretq_s16_u16(vaddl_u8(high_r, u8_zero));
        int16x8_t signed_low_r = vreinterpretq_s16_u16(vaddl_u8(low_r, u8_zero));
        int16x8_t signed_high_g = vreinterpretq_s16_u16(vaddl_u8(high_g, u8_zero));
        int16x8_t signed_low_g = vreinterpretq_s16_u16(vaddl_u8(low_g, u8_zero));
        int16x8_t signed_high_b = vreinterpretq_s16_u16(vaddl_u8(high_b, u8_zero));
        int16x8_t signed_low_b = vreinterpretq_s16_u16(vaddl_u8(low_b, u8_zero));

        // NOTE:
        // declaration may not appear after executable statement in block
        uint16x8_t high_y;
        uint16x8_t low_y;
        uint8x8_t scalar = vdup_n_u8(76);
        int16x8_t high_u;
        int16x8_t low_u;
        int16x8_t signed_scalar = vdupq_n_s16(-43);
        int16x8_t high_v;
        int16x8_t low_v;
        uint8x16x3_t pixel_yuv;
        int8x16_t u;
        int8x16_t v;

        // 1. Multiply transform matrix (Y′: unsigned, U/V: signed)
        high_y = vmull_u8(high_r, scalar);
        low_y = vmull_u8(low_r, scalar);

        high_u = vmulq_s16(signed_high_r, signed_scalar);
        low_u = vmulq_s16(signed_low_r, signed_scalar);

        signed_scalar = vdupq_n_s16(127);
        high_v = vmulq_s16(signed_high_r, signed_scalar);
        low_v = vmulq_s16(signed_low_r, signed_scalar);

        scalar = vdup_n_u8(150);
        high_y = vmlal_u8(high_y, high_g, scalar);
        low_y = vmlal_u8(low_y, low_g, scalar);

        signed_scalar = vdupq_n_s16(-84);
        high_u = vmlaq_s16(high_u, signed_high_g, signed_scalar);
        low_u = vmlaq_s16(low_u, signed_low_g, signed_scalar);

        signed_scalar = vdupq_n_s16(-106);
        high_v = vmlaq_s16(high_v, signed_high_g, signed_scalar);
        low_v = vmlaq_s16(low_v, signed_low_g, signed_scalar);

        scalar = vdup_n_u8(29);
        high_y = vmlal_u8(high_y, high_b, scalar);
        low_y = vmlal_u8(low_y, low_b, scalar);

        signed_scalar = vdupq_n_s16(127);
        high_u = vmlaq_s16(high_u, signed_high_b, signed_scalar);
        low_u = vmlaq_s16(low_u, signed_low_b, signed_scalar);

        signed_scalar = vdupq_n_s16(-21);
        high_v = vmlaq_s16(high_v, signed_high_b, signed_scalar);
        low_v = vmlaq_s16(low_v, signed_low_b, signed_scalar);

        // 2. Scale down (">>8") to 8-bit values with rounding ("+128") (Y′: unsigned, U/V: signed)
        // 3. Add an offset to the values to eliminate any negative values (all results are 8-bit unsigned)

        high_y = vaddq_u16(high_y, u16_rounding);
        low_y = vaddq_u16(low_y, u16_rounding);

        high_u = vaddq_s16(high_u, s16_rounding);
        low_u = vaddq_s16(low_u, s16_rounding);

        high_v = vaddq_s16(high_v, s16_rounding);
        low_v = vaddq_s16(low_v, s16_rounding);

        pixel_yuv.val[0] = vcombine_u8(vqshrn_n_u16(low_y, 8), vqshrn_n_u16(high_y, 8));

        u = vcombine_s8(vqshrn_n_s16(low_u, 8), vqshrn_n_s16(high_u, 8));

        v = vcombine_s8(vqshrn_n_s16(low_v, 8), vqshrn_n_s16(high_v, 8));

        u = vaddq_s8(u, s8_rounding);
        pixel_yuv.val[1] = vreinterpretq_u8_s8(u);

        v = vaddq_s8(v, s8_rounding);
        pixel_yuv.val[2] = vreinterpretq_u8_s8(v);

        // Store
        vst3q_u8(yuv, pixel_yuv);

        bgr += 3 * 16;
        yuv += 3 * 16;
    }

    // Handle leftovers
    for (i = count * 16; i < pixel_num; ++i) {
        uint8_t r = bgr[i * 3];
        uint8_t g = bgr[i * 3 + 1];
        uint8_t b = bgr[i * 3 + 2];

        uint16_t y_tmp = 76 * r + 150 * g + 29 * b;
        int16_t u_tmp = -43 * r - 84 * g + 127 * b;
        int16_t v_tmp = 127 * r - 106 * g - 21 * b;

        y_tmp = (y_tmp + 128) >> 8;
        u_tmp = (u_tmp + 128) >> 8;
        v_tmp = (v_tmp + 128) >> 8;

        yuv[i * 3] = (uint8_t) y_tmp;
        yuv[i * 3 + 1] = (uint8_t) (u_tmp + 128);
        yuv[i * 3 + 2] = (uint8_t) (v_tmp + 128);
    }
}

// Refer to: https://en.wikipedia.org/wiki/YUV (Full swing for BT.601)
void bgr888_to_yuv444_c(unsigned char *yuv, unsigned char *bgr, int pixel_num) {
    int i;
    for (i = 0; i < pixel_num; ++i) {
        uint8_t r = bgr[i * 3];
        uint8_t g = bgr[i * 3 + 1];
        uint8_t b = bgr[i * 3 + 2];

        // 1. Multiply transform matrix (Y′: unsigned, U/V: signed)
        uint16_t y_tmp = 76 * r + 150 * g + 29 * b;
        int16_t u_tmp = -43 * r - 84 * g + 127 * b;
        int16_t v_tmp = 127 * r - 106 * g - 21 * b;

        // 2. Scale down (">>8") to 8-bit values with rounding ("+128") (Y′: unsigned, U/V: signed)
        y_tmp = (y_tmp + 128) >> 8;
        u_tmp = (u_tmp + 128) >> 8;
        v_tmp = (v_tmp + 128) >> 8;

        // 3. Add an offset to the values to eliminate any negative values (all results are 8-bit unsigned)
        yuv[i * 3] = (uint8_t) y_tmp;
        yuv[i * 3 + 1] = (uint8_t) (u_tmp + 128);
        yuv[i * 3 + 2] = (uint8_t) (v_tmp + 128);
    }
}

void bgr888_to_yuv444_float(unsigned char *yuv, unsigned char *bgr, int pixel_num) {
    int i;
    for (i = 0; i < pixel_num; ++i) {
        uint8_t r = bgr[i * 3];
        uint8_t g = bgr[i * 3 + 1];
        uint8_t b = bgr[i * 3 + 2];

        uint8_t y = 0.299 * r + 0.587 * g + 0.114 * b;
        uint8_t u = -0.169 * r - 0.331 * g + 0.5 * b + 128;
        uint8_t v = 0.5 * r - 0.419 * g - 0.081 * b + 128;

        yuv[i * 3] = y;
        yuv[i * 3 + 1] = u;
        yuv[i * 3 + 2] = v;
    }
}


void reference_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int n) {
    int i;
    for (i = 0; i < n; i++) {
        int r = *src++; // load red
        int g = *src++; // load green
        int b = *src++; // load blue

        // build weighted average:
        int y = (r * 77) + (g * 151) + (b * 28);

        // undo the scale by 256 and write to memory:
        *dest++ = (y >> 8);
    }
}

//使用NEON Intrinsics优化
void neon_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int n) {
    int i;
    //读取8字节的预设值到64位寄存器
    uint8x8_t rfac = vdup_n_u8(77);// 转换权值 R
    uint8x8_t gfac = vdup_n_u8(151);// 转换权值 G
    uint8x8_t bfac = vdup_n_u8(28);// 转换权值 B
    n /= 8;

    for (i = 0; i < n; i++) {
        uint16x8_t temp;
        uint8x8x3_t rgb = vld3_u8 (src);//一次读取3个unit8x8到3个64位寄存器
        uint8x8_t result;

        temp = vmull_u8(rgb.val[0], rfac); // temp=rgb.val[0]*rfac
        temp = vmlal_u8(temp, rgb.val[1], gfac);// temp=temp+rgb.val[1]*gfac
        temp = vmlal_u8(temp, rgb.val[2], bfac);//temp=temp+rgb.val[2]*bfac

        result = vshrn_n_u16 (temp, 8); // 128位寄存器每16位右移第二个参数位
        vst1_u8 (dest, result); // 转存运算结果到dest
        src += 8 * 3;
        dest += 8;
    }
}

//NEON汇编代码优化:
void neon_asm_convert(uint8_t *__restrict dest, uint8_t *__restrict src, int numPixels) {
//    asm volatile("lsr          %2, %2, #3      \n"
//            "# build the three constants: \n"
//            "mov         r4, #28          \n" // Blue channel multiplier
//            "mov         r5, #151         \n" // Green channel multiplier
//            "mov         r6, #77          \n" // Red channel multiplier
//            "vdup.8      d4, r4           \n"
//            "vdup.8      d5, r5           \n"
//            "vdup.8      d6, r6           \n"
//            ".loop:                       \n"
//            "# load 8 pixels:             \n"
//            "vld4.8      {d0-d3}, [%1]!   \n"
//            "# do the weight average:     \n"
//            "vmull.u8    q7, d0, d4       \n"
//            "vmlal.u8    q7, d1, d5       \n"
//            "vmlal.u8    q7, d2, d6       \n"
//            "# shift and store:           \n"
//            "vshrn.u16   d7, q7, #8       \n" // Divide q3 by 256 and store in the d7
//            "vst1.8      {d7}, [%0]!      \n"
//            "subs        %2, %2, #1       \n" // Decrement iteration count
//            "bne         .loop            \n" // Repeat unil iteration count is not zero
//    :
//    : "r"(dest), "r"(src), "r"(numPixels)
//    : "r4", "r5", "r6"
//    );
}