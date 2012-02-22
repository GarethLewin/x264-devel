/*****************************************************************************
 * mc-c.c: x86 motion compensation
 *****************************************************************************
 * Copyright (C) 2003-2012 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at licensing@x264.com.
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "common/common.h"
#include "mc.h"

#define DECL_SUF( func, args )\
    void func##_mmx2 args;\
    void func##_sse2 args;\
    void func##_ssse3 args;

DECL_SUF( x264_pixel_avg_16x16, ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_16x8,  ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_8x16,  ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_8x8,   ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_8x4,   ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_4x16,  ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_4x8,   ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_4x4,   ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))
DECL_SUF( x264_pixel_avg_4x2,   ( pixel *, intptr_t, pixel *, intptr_t, pixel *, intptr_t, int ))

#define MC_WEIGHT(w,type) \
    void x264_mc_weight_w##w##_##type( pixel *, intptr_t, pixel *, intptr_t, const x264_weight_t *, int );

#define MC_WEIGHT_OFFSET(w,type) \
    void x264_mc_offsetadd_w##w##_##type( pixel *, intptr_t, pixel *, intptr_t, const x264_weight_t *, int ); \
    void x264_mc_offsetsub_w##w##_##type( pixel *, intptr_t, pixel *, intptr_t, const x264_weight_t *, int ); \
    MC_WEIGHT(w,type)

MC_WEIGHT_OFFSET( 4, mmx2 )
MC_WEIGHT_OFFSET( 8, mmx2 )
MC_WEIGHT_OFFSET( 12, mmx2 )
MC_WEIGHT_OFFSET( 16, mmx2 )
MC_WEIGHT_OFFSET( 20, mmx2 )
MC_WEIGHT_OFFSET( 12, sse2 )
MC_WEIGHT_OFFSET( 16, sse2 )
MC_WEIGHT_OFFSET( 20, sse2 )
#if HIGH_BIT_DEPTH
MC_WEIGHT_OFFSET( 8, sse2 )
#endif
MC_WEIGHT( 8, sse2  )
MC_WEIGHT( 4, ssse3 )
MC_WEIGHT( 8, ssse3 )
MC_WEIGHT( 12, ssse3 )
MC_WEIGHT( 16, ssse3 )
MC_WEIGHT( 20, ssse3 )
#undef MC_OFFSET
#undef MC_WEIGHT

void x264_mc_copy_w4_mmx  ( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_mc_copy_w8_mmx  ( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_mc_copy_w8_sse2 ( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_mc_copy_w16_mmx ( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_mc_copy_w16_sse2( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_mc_copy_w16_aligned_sse2( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_prefetch_fenc_420_mmx2( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_prefetch_fenc_422_mmx2( pixel *, intptr_t, pixel *, intptr_t, int );
void x264_prefetch_ref_mmx2( pixel *, intptr_t, int );
void x264_plane_copy_core_mmx2( pixel *, intptr_t, pixel *, intptr_t, int w, int h);
void x264_plane_copy_c( pixel *, intptr_t, pixel *, intptr_t, int w, int h );
void x264_plane_copy_interleave_core_mmx2( pixel *dst,  intptr_t i_dst,
                                           pixel *srcu, intptr_t i_srcu,
                                           pixel *srcv, intptr_t i_srcv, int w, int h );
void x264_plane_copy_interleave_core_sse2( pixel *dst,  intptr_t i_dst,
                                           pixel *srcu, intptr_t i_srcu,
                                           pixel *srcv, intptr_t i_srcv, int w, int h );
void x264_plane_copy_interleave_core_avx( pixel *dst,  intptr_t i_dst,
                                          pixel *srcu, intptr_t i_srcu,
                                          pixel *srcv, intptr_t i_srcv, int w, int h );
void x264_plane_copy_interleave_c( pixel *dst,  intptr_t i_dst,
                                   pixel *srcu, intptr_t i_srcu,
                                   pixel *srcv, intptr_t i_srcv, int w, int h );
void x264_plane_copy_deinterleave_mmx( pixel *dstu, intptr_t i_dstu,
                                       pixel *dstv, intptr_t i_dstv,
                                       pixel *src,  intptr_t i_src, int w, int h );
void x264_plane_copy_deinterleave_sse2( pixel *dstu, intptr_t i_dstu,
                                        pixel *dstv, intptr_t i_dstv,
                                        pixel *src,  intptr_t i_src, int w, int h );
void x264_plane_copy_deinterleave_ssse3( uint8_t *dstu, intptr_t i_dstu,
                                         uint8_t *dstv, intptr_t i_dstv,
                                         uint8_t *src,  intptr_t i_src, int w, int h );
void x264_plane_copy_deinterleave_avx( uint16_t *dstu, intptr_t i_dstu,
                                       uint16_t *dstv, intptr_t i_dstv,
                                       uint16_t *src,  intptr_t i_src, int w, int h );
void x264_store_interleave_chroma_mmx2( pixel *dst, intptr_t i_dst, pixel *srcu, pixel *srcv, int height );
void x264_store_interleave_chroma_sse2( pixel *dst, intptr_t i_dst, pixel *srcu, pixel *srcv, int height );
void x264_store_interleave_chroma_avx ( pixel *dst, intptr_t i_dst, pixel *srcu, pixel *srcv, int height );
void x264_load_deinterleave_chroma_fenc_mmx ( pixel *dst, pixel *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fenc_sse2( pixel *dst, pixel *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fenc_ssse3( uint8_t *dst, uint8_t *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fenc_avx( uint16_t *dst, uint16_t *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fdec_mmx ( pixel *dst, pixel *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fdec_sse2( pixel *dst, pixel *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fdec_ssse3( uint8_t *dst, uint8_t *src, intptr_t i_src, int height );
void x264_load_deinterleave_chroma_fdec_avx( uint16_t *dst, uint16_t *src, intptr_t i_src, int height );
void *x264_memcpy_aligned_mmx ( void *dst, const void *src, size_t n );
void *x264_memcpy_aligned_sse2( void *dst, const void *src, size_t n );
void x264_memzero_aligned_mmx ( void *dst, size_t n );
void x264_memzero_aligned_sse2( void *dst, size_t n );
void x264_integral_init4h_sse4( uint16_t *sum, uint8_t *pix, intptr_t stride );
void x264_integral_init8h_sse4( uint16_t *sum, uint8_t *pix, intptr_t stride );
void x264_integral_init8h_avx ( uint16_t *sum, uint8_t *pix, intptr_t stride );
void x264_integral_init4v_mmx  ( uint16_t *sum8, uint16_t *sum4, intptr_t stride );
void x264_integral_init4v_sse2 ( uint16_t *sum8, uint16_t *sum4, intptr_t stride );
void x264_integral_init4v_ssse3( uint16_t *sum8, uint16_t *sum4, intptr_t stride );
void x264_integral_init8v_mmx ( uint16_t *sum8, intptr_t stride );
void x264_integral_init8v_sse2( uint16_t *sum8, intptr_t stride );
void x264_mbtree_propagate_cost_sse2( int *dst, uint16_t *propagate_in, uint16_t *intra_costs,
                                      uint16_t *inter_costs, uint16_t *inv_qscales, float *fps_factor, int len );
void x264_mbtree_propagate_cost_avx ( int *dst, uint16_t *propagate_in, uint16_t *intra_costs,
                                      uint16_t *inter_costs, uint16_t *inv_qscales, float *fps_factor, int len );
void x264_mbtree_propagate_cost_fma4( int *dst, uint16_t *propagate_in, uint16_t *intra_costs,
                                      uint16_t *inter_costs, uint16_t *inv_qscales, float *fps_factor, int len );

#define MC_CHROMA(cpu)\
void x264_mc_chroma_##cpu( pixel *dstu, pixel *dstv, intptr_t i_dst, pixel *src, intptr_t i_src,\
                           int dx, int dy, int i_width, int i_height );
MC_CHROMA(mmx2)
MC_CHROMA(sse2)
MC_CHROMA(sse2_misalign)
MC_CHROMA(ssse3)
MC_CHROMA(ssse3_cache64)
MC_CHROMA(avx)
MC_CHROMA(avx_cache64)

#define LOWRES(cpu)\
void x264_frame_init_lowres_core_##cpu( pixel *src0, pixel *dst0, pixel *dsth, pixel *dstv, pixel *dstc,\
                                        intptr_t src_stride, intptr_t dst_stride, int width, int height );
LOWRES(mmx2)
LOWRES(cache32_mmx2)
LOWRES(sse2)
LOWRES(ssse3)
LOWRES(avx)
LOWRES(xop)

#define MC_WEIGHT_WTAB(function, instr, name1, name2, w12version)\
    static void (* x264_mc_##function##_wtab_##instr[6])( pixel *, intptr_t, pixel *, intptr_t, const x264_weight_t *, int ) =\
{\
    x264_mc_##function##_w4_##name1,\
    x264_mc_##function##_w4_##name1,\
    x264_mc_##function##_w8_##name2,\
    x264_mc_##function##_w##w12version##_##instr,\
    x264_mc_##function##_w16_##instr,\
    x264_mc_##function##_w20_##instr,\
};

#if HIGH_BIT_DEPTH
MC_WEIGHT_WTAB(weight,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(offsetadd,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(offsetsub,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(weight,sse2,mmx2,sse2,12)
MC_WEIGHT_WTAB(offsetadd,sse2,mmx2,sse2,16)
MC_WEIGHT_WTAB(offsetsub,sse2,mmx2,sse2,16)

static void x264_weight_cache_mmx2( x264_t *h, x264_weight_t *w )
{
    if( w->i_scale == 1<<w->i_denom )
    {
        if( w->i_offset < 0 )
            w->weightfn = h->mc.offsetsub;
        else
            w->weightfn = h->mc.offsetadd;
        for( int i = 0; i < 8; i++ )
            w->cachea[i] = abs(w->i_offset<<(BIT_DEPTH-8));
        return;
    }
    w->weightfn = h->mc.weight;
    int den1 = 1<<w->i_denom;
    int den2 = w->i_scale<<1;
    int den3 = 1+(w->i_offset<<(BIT_DEPTH-8+1));
    for( int i = 0; i < 8; i++ )
    {
        w->cachea[i] = den1;
        w->cacheb[i] = i&1 ? den3 : den2;
    }
}
#else
MC_WEIGHT_WTAB(weight,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(offsetadd,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(offsetsub,mmx2,mmx2,mmx2,12)
MC_WEIGHT_WTAB(weight,sse2,mmx2,sse2,16)
MC_WEIGHT_WTAB(offsetadd,sse2,mmx2,mmx2,16)
MC_WEIGHT_WTAB(offsetsub,sse2,mmx2,mmx2,16)
MC_WEIGHT_WTAB(weight,ssse3,ssse3,ssse3,16)

static void x264_weight_cache_mmx2( x264_t *h, x264_weight_t *w )
{
    int i;
    int16_t den1;

    if( w->i_scale == 1<<w->i_denom )
    {
        if( w->i_offset < 0 )
            w->weightfn = h->mc.offsetsub;
        else
            w->weightfn = h->mc.offsetadd;
        memset( w->cachea, abs(w->i_offset), sizeof(w->cachea) );
        return;
    }
    w->weightfn = h->mc.weight;
    den1 = 1 << (w->i_denom - 1) | w->i_offset << w->i_denom;
    for( i = 0; i < 8; i++ )
    {
        w->cachea[i] = w->i_scale;
        w->cacheb[i] = den1;
    }
}

static void x264_weight_cache_ssse3( x264_t *h, x264_weight_t *w )
{
    int i, den1;
    if( w->i_scale == 1<<w->i_denom )
    {
        if( w->i_offset < 0 )
            w->weightfn = h->mc.offsetsub;
        else
            w->weightfn = h->mc.offsetadd;

        memset( w->cachea, abs( w->i_offset ), sizeof(w->cachea) );
        return;
    }
    w->weightfn = h->mc.weight;
    den1 = w->i_scale << (8 - w->i_denom);
    for( i = 0; i < 8; i++ )
    {
        w->cachea[i] = den1;
        w->cacheb[i] = w->i_offset;
    }
}
#endif // !HIGH_BIT_DEPTH

const uint8_t x264_hpel_ref0[16] = {0,1,1,1,0,1,1,1,2,3,3,3,0,1,1,1};
const uint8_t x264_hpel_ref1[16] = {0,0,0,0,2,2,3,2,2,2,3,2,2,2,3,2};

#define MC_LUMA(name)\
void x264_mc_luma_##name( pixel *dst,    int i_dst_stride,\
                          pixel *src[4], int i_src_stride,\
                          int mvx, int mvy,\
                          int i_width, int i_height, const x264_weight_t *weight );\

MC_LUMA(mmx2)
MC_LUMA(sse2)
#if !HIGH_BIT_DEPTH
#if ARCH_X86
MC_LUMA(cache32_mmx2)
MC_LUMA(cache64_mmx2)
#endif
MC_LUMA(cache64_sse2)
MC_LUMA(cache64_ssse3)
#endif // !HIGH_BIT_DEPTH
MC_LUMA(sse2_bmi1)

#define GET_REF(name)\
pixel *x264_get_ref_##name( pixel *dst,   int *i_dst_stride,\
                            pixel *src[4], int i_src_stride,\
                            int mvx, int mvy,\
                            int i_width, int i_height, const x264_weight_t *weight );\

GET_REF(mmx2)
GET_REF(sse2)
#if !HIGH_BIT_DEPTH
#if ARCH_X86
GET_REF(cache32_mmx2)
GET_REF(cache64_mmx2)
#endif
GET_REF(sse2_misalign)
GET_REF(cache64_sse2)
GET_REF(cache64_ssse3)
#endif // !HIGH_BIT_DEPTH
GET_REF(sse2_bmi1)

#define HPEL(align, cpu, cpuv, cpuc, cpuh)\
void x264_hpel_filter_v_##cpuv( pixel *dst, pixel *src, int16_t *buf, intptr_t stride, intptr_t width);\
void x264_hpel_filter_c_##cpuc( pixel *dst, int16_t *buf, intptr_t width );\
void x264_hpel_filter_h_##cpuh( pixel *dst, pixel *src, intptr_t width );\
static void x264_hpel_filter_##cpu( pixel *dsth, pixel *dstv, pixel *dstc, pixel *src,\
                                    intptr_t stride, int width, int height, int16_t *buf )\
{\
    intptr_t realign = (intptr_t)src & (align-1);\
    src -= realign;\
    dstv -= realign;\
    dstc -= realign;\
    dsth -= realign;\
    width += realign;\
    while( height-- )\
    {\
        x264_hpel_filter_v_##cpuv( dstv, src, buf+8, stride, width );\
        x264_hpel_filter_c_##cpuc( dstc, buf+8, width );\
        x264_hpel_filter_h_##cpuh( dsth, src, width );\
        dsth += stride;\
        dstv += stride;\
        dstc += stride;\
        src  += stride;\
    }\
    x264_sfence();\
}

HPEL(8, mmx2, mmx2, mmx2, mmx2)
#if HIGH_BIT_DEPTH
HPEL(16, sse2, sse2, sse2, sse2)
#else // !HIGH_BIT_DEPTH
HPEL(16, sse2_amd, mmx2, mmx2, sse2)
#if ARCH_X86_64
void x264_hpel_filter_sse2 ( uint8_t *dsth, uint8_t *dstv, uint8_t *dstc, uint8_t *src, intptr_t stride, int width, int height, int16_t *buf );
void x264_hpel_filter_ssse3( uint8_t *dsth, uint8_t *dstv, uint8_t *dstc, uint8_t *src, intptr_t stride, int width, int height, int16_t *buf );
void x264_hpel_filter_avx  ( uint8_t *dsth, uint8_t *dstv, uint8_t *dstc, uint8_t *src, intptr_t stride, int width, int height, int16_t *buf );
#else
HPEL(16, sse2, sse2, sse2, sse2)
HPEL(16, ssse3, ssse3, ssse3, ssse3)
HPEL(16, avx, avx, avx, avx)
#endif
HPEL(16, sse2_misalign, sse2, sse2_misalign, sse2)
#endif // HIGH_BIT_DEPTH

static void x264_plane_copy_mmx2( pixel *dst, intptr_t i_dst, pixel *src, intptr_t i_src, int w, int h )
{
    int c_w = 16/sizeof(pixel) - 1;
    if( w < 256 ) { // tiny resolutions don't want non-temporal hints. dunno the exact threshold.
        x264_plane_copy_c( dst, i_dst, src, i_src, w, h );
    } else if( !(w&c_w) ) {
        x264_plane_copy_core_mmx2( dst, i_dst, src, i_src, w, h );
    } else if( i_src > 0 ) {
        // have to use plain memcpy on the last line (in memory order) to avoid overreading src
        x264_plane_copy_core_mmx2( dst, i_dst, src, i_src, (w+c_w)&~c_w, h-1 );
        memcpy( dst+i_dst*(h-1), src+i_src*(h-1), w*sizeof(pixel) );
    } else {
        memcpy( dst, src, w*sizeof(pixel) );
        x264_plane_copy_core_mmx2( dst+i_dst, i_dst, src+i_src, i_src, (w+c_w)&~c_w, h-1 );
    }
}

#define PLANE_INTERLEAVE(cpu) \
static void x264_plane_copy_interleave_##cpu( pixel *dst,  intptr_t i_dst,\
                                              pixel *srcu, intptr_t i_srcu,\
                                              pixel *srcv, intptr_t i_srcv, int w, int h )\
{\
    if( !(w&15) ) {\
        x264_plane_copy_interleave_core_##cpu( dst, i_dst, srcu, i_srcu, srcv, i_srcv, w, h );\
    } else if( w < 16 || (i_srcu ^ i_srcv) ) {\
        x264_plane_copy_interleave_c( dst, i_dst, srcu, i_srcu, srcv, i_srcv, w, h );\
    } else if( i_srcu > 0 ) {\
        x264_plane_copy_interleave_core_##cpu( dst, i_dst, srcu, i_srcu, srcv, i_srcv, (w+15)&~15, h-1 );\
        x264_plane_copy_interleave_c( dst+i_dst*(h-1), 0, srcu+i_srcu*(h-1), 0, srcv+i_srcv*(h-1), 0, w, 1 );\
    } else {\
        x264_plane_copy_interleave_c( dst, 0, srcu, 0, srcv, 0, w, 1 );\
        x264_plane_copy_interleave_core_##cpu( dst+i_dst, i_dst, srcu+i_srcu, i_srcu, srcv+i_srcv, i_srcv, (w+15)&~15, h-1 );\
    }\
}

PLANE_INTERLEAVE(mmx2)
PLANE_INTERLEAVE(sse2)
#if HIGH_BIT_DEPTH
PLANE_INTERLEAVE(avx)
#endif

void x264_mc_init_mmx( int cpu, x264_mc_functions_t *pf )
{
    if( !(cpu&X264_CPU_MMX) )
        return;

    pf->load_deinterleave_chroma_fenc = x264_load_deinterleave_chroma_fenc_mmx;
    pf->load_deinterleave_chroma_fdec = x264_load_deinterleave_chroma_fdec_mmx;

    pf->plane_copy_deinterleave = x264_plane_copy_deinterleave_mmx;

    pf->copy_16x16_unaligned = x264_mc_copy_w16_mmx;
    pf->copy[PIXEL_16x16] = x264_mc_copy_w16_mmx;
    pf->copy[PIXEL_8x8]   = x264_mc_copy_w8_mmx;
    pf->copy[PIXEL_4x4]   = x264_mc_copy_w4_mmx;
    pf->memcpy_aligned  = x264_memcpy_aligned_mmx;
    pf->memzero_aligned = x264_memzero_aligned_mmx;
    pf->integral_init4v = x264_integral_init4v_mmx;
    pf->integral_init8v = x264_integral_init8v_mmx;

    if( !(cpu&X264_CPU_MMX2) )
        return;

    pf->prefetch_fenc_420 = x264_prefetch_fenc_420_mmx2;
    pf->prefetch_fenc_422 = x264_prefetch_fenc_422_mmx2;
    pf->prefetch_ref  = x264_prefetch_ref_mmx2;

    pf->plane_copy = x264_plane_copy_mmx2;
    pf->plane_copy_interleave = x264_plane_copy_interleave_mmx2;
    pf->store_interleave_chroma = x264_store_interleave_chroma_mmx2;

    pf->avg[PIXEL_16x16] = x264_pixel_avg_16x16_mmx2;
    pf->avg[PIXEL_16x8]  = x264_pixel_avg_16x8_mmx2;
    pf->avg[PIXEL_8x16]  = x264_pixel_avg_8x16_mmx2;
    pf->avg[PIXEL_8x8]   = x264_pixel_avg_8x8_mmx2;
    pf->avg[PIXEL_8x4]   = x264_pixel_avg_8x4_mmx2;
    pf->avg[PIXEL_4x16]  = x264_pixel_avg_4x16_mmx2;
    pf->avg[PIXEL_4x8]   = x264_pixel_avg_4x8_mmx2;
    pf->avg[PIXEL_4x4]   = x264_pixel_avg_4x4_mmx2;
    pf->avg[PIXEL_4x2]   = x264_pixel_avg_4x2_mmx2;

    pf->mc_luma = x264_mc_luma_mmx2;
    pf->get_ref = x264_get_ref_mmx2;
    pf->mc_chroma = x264_mc_chroma_mmx2;
    pf->hpel_filter = x264_hpel_filter_mmx2;
    pf->weight = x264_mc_weight_wtab_mmx2;
    pf->weight_cache = x264_weight_cache_mmx2;
    pf->offsetadd = x264_mc_offsetadd_wtab_mmx2;
    pf->offsetsub = x264_mc_offsetsub_wtab_mmx2;

    pf->frame_init_lowres_core = x264_frame_init_lowres_core_mmx2;

#if HIGH_BIT_DEPTH
#if ARCH_X86 // all x86_64 cpus with cacheline split issues use sse2 instead
    if( cpu&(X264_CPU_CACHELINE_32|X264_CPU_CACHELINE_64) )
        pf->frame_init_lowres_core = x264_frame_init_lowres_core_cache32_mmx2;
#endif

    if( !(cpu&X264_CPU_SSE2) )
        return;

    pf->frame_init_lowres_core = x264_frame_init_lowres_core_sse2;

    pf->load_deinterleave_chroma_fenc = x264_load_deinterleave_chroma_fenc_sse2;
    pf->load_deinterleave_chroma_fdec = x264_load_deinterleave_chroma_fdec_sse2;

    pf->plane_copy_interleave   = x264_plane_copy_interleave_sse2;
    pf->plane_copy_deinterleave = x264_plane_copy_deinterleave_sse2;

    if( cpu&X264_CPU_SSE2_IS_FAST )
    {
        pf->get_ref = x264_get_ref_sse2;
        pf->mc_luma = x264_mc_luma_sse2;
        pf->hpel_filter = x264_hpel_filter_sse2;
    }

    pf->memcpy_aligned  = x264_memcpy_aligned_sse2;
    pf->memzero_aligned = x264_memzero_aligned_sse2;
    pf->integral_init4v = x264_integral_init4v_sse2;
    pf->integral_init8v = x264_integral_init8v_sse2;
    pf->mbtree_propagate_cost = x264_mbtree_propagate_cost_sse2;
    pf->store_interleave_chroma = x264_store_interleave_chroma_sse2;
    pf->offsetadd = x264_mc_offsetadd_wtab_sse2;
    pf->offsetsub = x264_mc_offsetsub_wtab_sse2;

    if( cpu&X264_CPU_SSE2_IS_SLOW )
        return;

    pf->avg[PIXEL_16x16] = x264_pixel_avg_16x16_sse2;
    pf->avg[PIXEL_16x8]  = x264_pixel_avg_16x8_sse2;
    pf->avg[PIXEL_8x16]  = x264_pixel_avg_8x16_sse2;
    pf->avg[PIXEL_8x8]   = x264_pixel_avg_8x8_sse2;
    pf->avg[PIXEL_8x4]   = x264_pixel_avg_8x4_sse2;
    pf->avg[PIXEL_4x16]  = x264_pixel_avg_4x16_sse2;
    pf->avg[PIXEL_4x8]   = x264_pixel_avg_4x8_sse2;
    pf->avg[PIXEL_4x4]   = x264_pixel_avg_4x4_sse2;
    pf->avg[PIXEL_4x2]   = x264_pixel_avg_4x2_sse2;

    pf->copy[PIXEL_16x16] = x264_mc_copy_w16_aligned_sse2;
    pf->weight = x264_mc_weight_wtab_sse2;

    if( !(cpu&X264_CPU_STACK_MOD4) )
        pf->mc_chroma = x264_mc_chroma_sse2;

    if( !(cpu&X264_CPU_SSSE3) )
        return;

    pf->frame_init_lowres_core = x264_frame_init_lowres_core_ssse3;

    if( (cpu&X264_CPU_SHUFFLE_IS_FAST) && !(cpu&X264_CPU_SLOW_ATOM) )
        pf->integral_init4v = x264_integral_init4v_ssse3;

    if( !(cpu&X264_CPU_AVX) )
        return;

    pf->frame_init_lowres_core = x264_frame_init_lowres_core_avx;
    pf->load_deinterleave_chroma_fenc = x264_load_deinterleave_chroma_fenc_avx;
    pf->load_deinterleave_chroma_fdec = x264_load_deinterleave_chroma_fdec_avx;
    pf->plane_copy_interleave        = x264_plane_copy_interleave_avx;
    pf->plane_copy_deinterleave      = x264_plane_copy_deinterleave_avx;
    pf->store_interleave_chroma      = x264_store_interleave_chroma_avx;

    if( !(cpu&X264_CPU_STACK_MOD4) )
        pf->mc_chroma = x264_mc_chroma_avx;

    if( cpu&X264_CPU_XOP )
        pf->frame_init_lowres_core = x264_frame_init_lowres_core_xop;
#else // !HIGH_BIT_DEPTH

#if ARCH_X86 // all x86_64 cpus with cacheline split issues use sse2 instead
    if( cpu&X264_CPU_CACHELINE_32 )
    {
        pf->mc_luma = x264_mc_luma_cache32_mmx2;
        pf->get_ref = x264_get_ref_cache32_mmx2;
        pf->frame_init_lowres_core = x264_frame_init_lowres_core_cache32_mmx2;
    }
    else if( cpu&X264_CPU_CACHELINE_64 )
    {
        pf->mc_luma = x264_mc_luma_cache64_mmx2;
        pf->get_ref = x264_get_ref_cache64_mmx2;
        pf->frame_init_lowres_core = x264_frame_init_lowres_core_cache32_mmx2;
    }
#endif

    if( !(cpu&X264_CPU_SSE2) )
        return;

    pf->memcpy_aligned = x264_memcpy_aligned_sse2;
    pf->memzero_aligned = x264_memzero_aligned_sse2;
    pf->integral_init4v = x264_integral_init4v_sse2;
    pf->integral_init8v = x264_integral_init8v_sse2;
    pf->hpel_filter = x264_hpel_filter_sse2_amd;
    pf->mbtree_propagate_cost = x264_mbtree_propagate_cost_sse2;

    if( cpu&X264_CPU_SSE2_IS_SLOW )
        return;

    pf->weight = x264_mc_weight_wtab_sse2;
    if( !(cpu&X264_CPU_SLOW_ATOM) )
    {
        pf->offsetadd = x264_mc_offsetadd_wtab_sse2;
        pf->offsetsub = x264_mc_offsetsub_wtab_sse2;
    }

    pf->copy[PIXEL_16x16] = x264_mc_copy_w16_aligned_sse2;
    pf->avg[PIXEL_16x16] = x264_pixel_avg_16x16_sse2;
    pf->avg[PIXEL_16x8]  = x264_pixel_avg_16x8_sse2;
    pf->avg[PIXEL_8x16] = x264_pixel_avg_8x16_sse2;
    pf->avg[PIXEL_8x8]  = x264_pixel_avg_8x8_sse2;
    pf->avg[PIXEL_8x4]  = x264_pixel_avg_8x4_sse2;
    pf->hpel_filter = x264_hpel_filter_sse2;
    if( cpu&X264_CPU_SSE_MISALIGN )
        pf->hpel_filter = x264_hpel_filter_sse2_misalign;
    pf->frame_init_lowres_core = x264_frame_init_lowres_core_sse2;
    if( !(cpu&X264_CPU_STACK_MOD4) )
        pf->mc_chroma = x264_mc_chroma_sse2;

    if( cpu&X264_CPU_SSE2_IS_FAST )
    {
        pf->store_interleave_chroma = x264_store_interleave_chroma_sse2; // FIXME sse2fast? sse2medium?
        pf->load_deinterleave_chroma_fenc = x264_load_deinterleave_chroma_fenc_sse2;
        pf->load_deinterleave_chroma_fdec = x264_load_deinterleave_chroma_fdec_sse2;
        pf->plane_copy_interleave   = x264_plane_copy_interleave_sse2;
        pf->plane_copy_deinterleave = x264_plane_copy_deinterleave_sse2;
        pf->mc_luma = x264_mc_luma_sse2;
        pf->get_ref = x264_get_ref_sse2;
        if( cpu&X264_CPU_CACHELINE_64 )
        {
            pf->mc_luma = x264_mc_luma_cache64_sse2;
            pf->get_ref = x264_get_ref_cache64_sse2;
        }
        if( cpu&X264_CPU_SSE_MISALIGN )
        {
            pf->get_ref = x264_get_ref_sse2_misalign;
            if( !(cpu&X264_CPU_STACK_MOD4) )
                pf->mc_chroma = x264_mc_chroma_sse2_misalign;
        }
    }

    if( !(cpu&X264_CPU_SSSE3) )
        return;

    pf->avg[PIXEL_16x16] = x264_pixel_avg_16x16_ssse3;
    pf->avg[PIXEL_16x8]  = x264_pixel_avg_16x8_ssse3;
    pf->avg[PIXEL_8x16]  = x264_pixel_avg_8x16_ssse3;
    pf->avg[PIXEL_8x8]   = x264_pixel_avg_8x8_ssse3;
    pf->avg[PIXEL_8x4]   = x264_pixel_avg_8x4_ssse3;
    pf->avg[PIXEL_4x16]  = x264_pixel_avg_4x16_ssse3;
    pf->avg[PIXEL_4x8]   = x264_pixel_avg_4x8_ssse3;
    pf->avg[PIXEL_4x4]   = x264_pixel_avg_4x4_ssse3;
    pf->avg[PIXEL_4x2]   = x264_pixel_avg_4x2_ssse3;

    pf->load_deinterleave_chroma_fenc = x264_load_deinterleave_chroma_fenc_ssse3;
    pf->load_deinterleave_chroma_fdec = x264_load_deinterleave_chroma_fdec_ssse3;
    pf->plane_copy_deinterleave = x264_plane_copy_deinterleave_ssse3;

    pf->hpel_filter = x264_hpel_filter_ssse3;
    pf->frame_init_lowres_core = x264_frame_init_lowres_core_ssse3;
    if( !(cpu&X264_CPU_STACK_MOD4) )
        pf->mc_chroma = x264_mc_chroma_ssse3;

    if( cpu&X264_CPU_CACHELINE_64 )
    {
        if( !(cpu&X264_CPU_STACK_MOD4) )
            pf->mc_chroma = x264_mc_chroma_ssse3_cache64;
        pf->mc_luma = x264_mc_luma_cache64_ssse3;
        pf->get_ref = x264_get_ref_cache64_ssse3;

        /* ssse3 weight is slower on Nehalem, so only assign here. */
        pf->weight_cache = x264_weight_cache_ssse3;
        pf->weight = x264_mc_weight_wtab_ssse3;
    }

    if( (cpu&X264_CPU_SHUFFLE_IS_FAST) && !(cpu&X264_CPU_SLOW_ATOM) )
        pf->integral_init4v = x264_integral_init4v_ssse3;

    if( !(cpu&X264_CPU_SSE4) )
        return;

    pf->integral_init4h = x264_integral_init4h_sse4;
    pf->integral_init8h = x264_integral_init8h_sse4;

    if( !(cpu&X264_CPU_AVX) )
        return;

    pf->frame_init_lowres_core = x264_frame_init_lowres_core_avx;
    pf->integral_init8h = x264_integral_init8h_avx;
    pf->hpel_filter = x264_hpel_filter_avx;

    /* ssse3 weight seems to be faster again on Sandy Bridge and Bulldozer. */
    pf->weight_cache = x264_weight_cache_ssse3;
    pf->weight = x264_mc_weight_wtab_ssse3;
    if( !(cpu&X264_CPU_STACK_MOD4) )
        pf->mc_chroma = x264_mc_chroma_avx;

    if( cpu&X264_CPU_XOP )
        pf->frame_init_lowres_core = x264_frame_init_lowres_core_xop;
#endif // HIGH_BIT_DEPTH

    if( !(cpu&X264_CPU_AVX) )
        return;
    pf->mbtree_propagate_cost = x264_mbtree_propagate_cost_avx;

    if( !(cpu&X264_CPU_FMA4) )
        return;
    pf->mbtree_propagate_cost = x264_mbtree_propagate_cost_fma4;
}
