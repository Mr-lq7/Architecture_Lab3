/********************************************************
 * Performance Lab: kenels.c
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "林清泉",     /* 姓名 */
    "SZ170110314",  /* 学号 */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	    for (j = 0; j < dim; j++)
	        dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];


}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    // naive_rotate(dim, src, dst);
    //分块,目的是提高Cache命中率
    //消除循环的低效率,循环展开
    //减少过程调用
    //消除不必要的寄存器引用
    // int t,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;
    int block = 32;
    int i, j, k = 0;
    int ii, jj;
    for (ii = 0; ii < dim; ii += block) {
        // t = ii+0, t1 = ii+1,t2=ii+2,t3=ii+3,t4=ii+4,t5=ii+5,t6=ii+6,t7=ii+7;
        // t8=ii+8,t9=ii+9,t10=ii+10,t11=ii+11,t12=ii+12,t13=ii+13,t14=ii+14,t15=ii+15;
        // t16=ii+16,t17=ii+17,t18=ii+18,t19=ii+19,t20=ii+20,t21=ii+21,t22=ii+22,t23=ii+23;
        // t24=ii+24,t25=ii+25,t26=ii+26,t27=ii+27,t28=ii+28,t29=ii+29,t30=ii+30,t31=ii+31;
        for (jj = 0; jj < dim; jj++) {
            // for (int k = 0; k < block; ++k) {
                dst[RIDX(dim-1-jj, ii+k, dim)] = src[RIDX(ii+k, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+1, dim)] = src[RIDX(ii+k+1, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+2, dim)] = src[RIDX(ii+k+2, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+3, dim)] = src[RIDX(ii+k+3, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+4, dim)] = src[RIDX(ii+k+4, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+5, dim)] = src[RIDX(ii+k+5, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+6, dim)] = src[RIDX(ii+k+6, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+7, dim)] = src[RIDX(ii+k+7, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+8, dim)] = src[RIDX(ii+k+8, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+9, dim)] = src[RIDX(ii+k+9, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+10, dim)] = src[RIDX(ii+k+10, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+11, dim)] = src[RIDX(ii+k+11, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+12, dim)] = src[RIDX(ii+k+12, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+13, dim)] = src[RIDX(ii+k+13, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+14, dim)] = src[RIDX(ii+k+14, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+15, dim)] = src[RIDX(ii+k+15, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+16, dim)] = src[RIDX(ii+k+16, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+17, dim)] = src[RIDX(ii+k+17, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+18, dim)] = src[RIDX(ii+k+18, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+19, dim)] = src[RIDX(ii+k+19, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+20, dim)] = src[RIDX(ii+k+20, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+21, dim)] = src[RIDX(ii+k+21, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+22, dim)] = src[RIDX(ii+k+22, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+23, dim)] = src[RIDX(ii+k+23, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+24, dim)] = src[RIDX(ii+k+24, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+25, dim)] = src[RIDX(ii+k+25, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+26, dim)] = src[RIDX(ii+k+26, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+27, dim)] = src[RIDX(ii+k+27, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+28, dim)] = src[RIDX(ii+k+28, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+29, dim)] = src[RIDX(ii+k+29, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+30, dim)] = src[RIDX(ii+k+30, jj, dim)];
                dst[RIDX(dim-1-jj, ii+k+31, dim)] = src[RIDX(ii+k+31, jj, dim)];

                // dst[RIDX(dim-1-jj, t, dim)] = src[RIDX(t, jj, dim)];
                // dst[RIDX(dim-1-jj, t1, dim)] = src[RIDX(t1, jj, dim)];
                // dst[RIDX(dim-1-jj, t2, dim)] = src[RIDX(t2, jj, dim)];
                // dst[RIDX(dim-1-jj, t3, dim)] = src[RIDX(t3, jj, dim)];
                // dst[RIDX(dim-1-jj, t4, dim)] = src[RIDX(t4, jj, dim)];
                // dst[RIDX(dim-1-jj, t5, dim)] = src[RIDX(t5, jj, dim)];
                // dst[RIDX(dim-1-jj, t6, dim)] = src[RIDX(t6, jj, dim)];
                // dst[RIDX(dim-1-jj, t7, dim)] = src[RIDX(t7, jj, dim)];
                // dst[RIDX(dim-1-jj, t8, dim)] = src[RIDX(t8, jj, dim)];
                // dst[RIDX(dim-1-jj, t9, dim)] = src[RIDX(t9, jj, dim)];
                // dst[RIDX(dim-1-jj, t10, dim)] = src[RIDX(t10, jj, dim)];
                // dst[RIDX(dim-1-jj, t11, dim)] = src[RIDX(t11, jj, dim)];
                // dst[RIDX(dim-1-jj, t12, dim)] = src[RIDX(t12, jj, dim)];
                // dst[RIDX(dim-1-jj, t13, dim)] = src[RIDX(t13, jj, dim)];
                // dst[RIDX(dim-1-jj, t14, dim)] = src[RIDX(t14, jj, dim)];
                // dst[RIDX(dim-1-jj, t15, dim)] = src[RIDX(t15, jj, dim)];
                // dst[RIDX(dim-1-jj, t16, dim)] = src[RIDX(t16, jj, dim)];
                // dst[RIDX(dim-1-jj, t17, dim)] = src[RIDX(t17, jj, dim)];
                // dst[RIDX(dim-1-jj, t18, dim)] = src[RIDX(t18, jj, dim)];
                // dst[RIDX(dim-1-jj, t19, dim)] = src[RIDX(t19, jj, dim)];
                // dst[RIDX(dim-1-jj, t20, dim)] = src[RIDX(t20, jj, dim)];
                // dst[RIDX(dim-1-jj, t21, dim)] = src[RIDX(t21, jj, dim)];
                // dst[RIDX(dim-1-jj, t22, dim)] = src[RIDX(t22, jj, dim)];
                // dst[RIDX(dim-1-jj, t23, dim)] = src[RIDX(t23, jj, dim)];
                // dst[RIDX(dim-1-jj, t24, dim)] = src[RIDX(t24, jj, dim)];
                // dst[RIDX(dim-1-jj, t25, dim)] = src[RIDX(t25, jj, dim)];
                // dst[RIDX(dim-1-jj, t26, dim)] = src[RIDX(t26, jj, dim)];
                // dst[RIDX(dim-1-jj, t27, dim)] = src[RIDX(t27, jj, dim)];
                // dst[RIDX(dim-1-jj, t28, dim)] = src[RIDX(t28, jj, dim)];
                // dst[RIDX(dim-1-jj, t29, dim)] = src[RIDX(t29, jj, dim)];
                // dst[RIDX(dim-1-jj, t30, dim)] = src[RIDX(t30, jj, dim)];
                // dst[RIDX(dim-1-jj, t31, dim)] = src[RIDX(t31, jj, dim)];

                // dst[RIDX(dim-1-jj, ii+k, dim)] = src[RIDX(ii+k, jj, dim)];
                // dst[(dim-1-jj)*dim+ii+k] = src[(ii+k)*dim+jj];
                // *(dst+RIDX(dim-1-jj, ii+k, dim)) = *(src+RIDX(ii+k, jj, dim));
                // *(dst+(dim-1-jj)*dim+ii+k) = *(src+(ii+k)*dim+jj);
            // }
        }
    }
    // for (ii = 0; ii < dim; ii += block) {
    //     for (jj = 0; jj < dim; jj += block) {
    //         for (i = ii; i < ii + block; i++) {
    //             for (j = jj; j < jj + block; j++) {
    //                 // dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
    //                 // dst[RIDX(dim-1-j-1, i, dim)] = src[RIDX(i, j+1, dim)];
    //                 // dst[RIDX(dim-1-j-2, i, dim)] = src[RIDX(i, j+2, dim)];
    //                 // dst[RIDX(dim-1-j-3, i, dim)] = src[RIDX(i, j+3, dim)];
    //                 // dst[RIDX(dim-1-j, i+1, dim)] = src[RIDX(i+1, j, dim)];
    //                 // dst[RIDX(dim-1-j-1, i+1, dim)] = src[RIDX(i+1, j+1, dim)];
    //                 // dst[RIDX(dim-1-j-2, i+1, dim)] = src[RIDX(i+1, j+2, dim)];
    //                 // dst[RIDX(dim-1-j-3, i+1, dim)] = src[RIDX(i+1, j+3, dim)];
    //                 // dst[RIDX(dim-1-j, i+2, dim)] = src[RIDX(i+2, j, dim)];
    //                 // dst[RIDX(dim-1-j-1, i+2, dim)] = src[RIDX(i+2, j+1, dim)];
    //                 // dst[RIDX(dim-1-j-2, i+2, dim)] = src[RIDX(i+2, j+2, dim)];
    //                 // dst[RIDX(dim-1-j-3, i+2, dim)] = src[RIDX(i+2, j+3, dim)];
    //                 // dst[RIDX(dim-1-j, i+3, dim)] = src[RIDX(i+3, j, dim)];
    //                 // dst[RIDX(dim-1-j-1, i+3, dim)] = src[RIDX(i+3, j+1, dim)];
    //                 // dst[RIDX(dim-1-j-2, i+3, dim)] = src[RIDX(i+3, j+2, dim)];
    //                 // dst[RIDX(dim-1-j-3, i+3, dim)] = src[RIDX(i+3, j+3, dim)];

    //                 dst[RIDX(dim-1-i, j, dim)] = src[RIDX(j, i, dim)];
    //                 // dst[RIDX(dim-1-i, j+1, dim)] = src[RIDX(j+1, i, dim)];
    //                 // dst[RIDX(dim-1-i, j+2, dim)] = src[RIDX(j+2, i, dim)];
    //                 // dst[RIDX(dim-1-i, j+3, dim)] = src[RIDX(j+3, i, dim)];
    //                 // dst[RIDX(dim-1-i-1, j, dim)] = src[RIDX(j, i+1, dim)];
    //                 // dst[RIDX(dim-1-i-1, j+1, dim)] = src[RIDX(j+1, i+1, dim)];
    //                 // dst[RIDX(dim-1-i-1, j+2, dim)] = src[RIDX(j+2, i+1, dim)];
    //                 // dst[RIDX(dim-1-i-1, j+3, dim)] = src[RIDX(j+3, i+1, dim)];
    //                 // dst[RIDX(dim-1-i-2, j, dim)] = src[RIDX(j, i+2, dim)];
    //                 // dst[RIDX(dim-1-i-2, j+1, dim)] = src[RIDX(j+1, i+2, dim)];
    //                 // dst[RIDX(dim-1-i-2, j+2, dim)] = src[RIDX(j+2, i+2, dim)];
    //                 // dst[RIDX(dim-1-i-2, j+3, dim)] = src[RIDX(j+3, i+2, dim)];
    //                 // dst[RIDX(dim-1-i-3, j, dim)] = src[RIDX(j, i+3, dim)];
    //                 // dst[RIDX(dim-1-i-3, j+1, dim)] = src[RIDX(j+1, i+3, dim)];
    //                 // dst[RIDX(dim-1-i-3, j+2, dim)] = src[RIDX(j+2, i+3, dim)];
    //                 // dst[RIDX(dim-1-i-3, j+3, dim)] = src[RIDX(j+3, i+3, dim)];

    //                 // dst[(dim-1-i)*dim+j] = src[j*dim+i];


    //             }
    //         }
    //     }
    // }
}

/*********************************************************************
 * register_rotate_functions - 通过为每一个测试函数调用add_rotate_function(),
 * 登记你所有不同版本的旋转代码到评测程序driver中，
 * 当你运行driver程序时，它将测试并且给出每一个已经登记的测试函数的性能。
 *********************************************************************/

void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);   
    /* ... Register additional test functions here */
}


/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static inline int min(int a, int b) { return (a < b ? a : b); }
static inline int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) 
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */
char smooth_descr[] = "smooth: Current working version";
int pre_R[555][555], pre_G[555][555], pre_B[555][555];
void smooth(int dim, pixel *src, pixel *dst) 
{
    // naive_smooth(dim, src, dst);
    // 预处理,计算矩阵的前缀和f[i][j]表示从(0,0)-(i-1,j-1)全部元素的和
    // dp[i][j] = dp[i-1][j] + dp[i][j-1] - dp[i-1][j-1] + data[i-1][j-1]
    int i, j;
    for (i = 1; i <= dim; ++i) {
        for (j = 1; j <= dim; ++j) {
            // pre_R[j][i] = pre_R[j-1][i] + pre_R[j][i-1] - pre_R[j-1][i-1] + src[RIDX(j-1,i-1,dim)].red;
            // pre_G[j][i] = pre_G[j-1][i] + pre_G[j][i-1] - pre_G[j-1][i-1] + src[RIDX(j-1,i-1,dim)].green;
            // pre_B[j][i] = pre_B[j-1][i] + pre_B[j][i-1] - pre_B[j-1][i-1] + src[RIDX(j-1,i-1,dim)].blue;
            
            // pre_R[i][j] = pre_R[i-1][j] + pre_R[i][j-1] - pre_R[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->red;
            // pre_G[i][j] = pre_G[i-1][j] + pre_G[i][j-1] - pre_G[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->green;
            // pre_B[i][j] = pre_B[i-1][j] + pre_B[i][j-1] - pre_B[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->blue;
            
            // *(pre_R+i*555+j) = *(pre_R+(i-1)*555+j) + *(pre_R+i*555+j-1) - *(pre_R+(i-1)*555+j-1) + (src+RIDX(i-1,j-1,dim))->red;
            // *(pre_G+i*555+j) = *(pre_G+(i-1)*555+j) + *(pre_G+i*555+j-1) - *(pre_G+(i-1)*555+j-1) + (src+RIDX(i-1,j-1,dim))->red;
            // *(pre_B+i*555+j) = *(pre_B+(i-1)*555+j) + *(pre_B+i*555+j-1) - *(pre_B+(i-1)*555+j-1) + (src+RIDX(i-1,j-1,dim))->red;
            pre_R[i][j] = pre_R[i-1][j] + pre_R[i][j-1] - pre_R[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->red;
            pre_G[i][j] = pre_G[i-1][j] + pre_G[i][j-1] - pre_G[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->green;
            pre_B[i][j] = pre_B[i-1][j] + pre_B[i][j-1] - pre_B[i-1][j-1] + (src+RIDX(i-1,j-1,dim))->blue;
        }
    }
    //左上
    // dst[0].red = pre_R[2][2] >> 2;
    // dst[0].green = pre_G[2][2] >> 2;
    // dst[0].blue = pre_B[2][2] >> 2;
    (dst+0)->red = pre_R[2][2] >> 2;
    (dst+0)->green = pre_G[2][2] >> 2;
    (dst+0)->blue = pre_B[2][2] >> 2;
    //右上
    // dst[dim-1].red = (pre_R[2][dim] - pre_R[2][dim-2]) >> 2;
    // dst[dim-1].green = (pre_G[2][dim] - pre_G[2][dim-2]) >> 2;
    // dst[dim-1].blue = (pre_B[2][dim] - pre_B[2][dim-2]) >> 2;
    (dst+dim-1)->red = (pre_R[2][dim] - pre_R[2][dim-2]) >> 2;
    (dst+dim-1)->green = (pre_G[2][dim] - pre_G[2][dim-2]) >> 2;
    (dst+dim-1)->blue = (pre_B[2][dim] - pre_B[2][dim-2]) >> 2;
    //左下
    // dst[RIDX(dim-1, 0, dim)].red = (pre_R[dim][2] - pre_R[dim-2][2]) >> 2;
    // dst[RIDX(dim-1, 0, dim)].green = (pre_G[dim][2] - pre_G[dim-2][2]) >> 2;
    // dst[RIDX(dim-1, 0, dim)].blue = (pre_B[dim][2] - pre_B[dim-2][2]) >> 2;
    (dst+RIDX(dim-1, 0, dim))->red = (pre_R[dim][2] - pre_R[dim-2][2]) >> 2;
    (dst+RIDX(dim-1, 0, dim))->green = (pre_G[dim][2] - pre_G[dim-2][2]) >> 2;
    (dst+RIDX(dim-1, 0, dim))->blue = (pre_B[dim][2] - pre_B[dim-2][2]) >> 2;

    //右下
    // dst[RIDX(dim-1, dim-1, dim)].red = (pre_R[dim][dim] - pre_R[dim-2][dim] - pre_R[dim][dim-2] + pre_R[dim-2][dim-2]) >> 2;
    // dst[RIDX(dim-1, dim-1, dim)].green = (pre_G[dim][dim] - pre_G[dim-2][dim] - pre_G[dim][dim-2] + pre_G[dim-2][dim-2]) >> 2;
    // dst[RIDX(dim-1, dim-1, dim)].blue = (pre_B[dim][dim] - pre_B[dim-2][dim] - pre_B[dim][dim-2] + pre_B[dim-2][dim-2]) >> 2;

    (dst+RIDX(dim-1, dim-1, dim))->red = (pre_R[dim][dim] - pre_R[dim-2][dim] - pre_R[dim][dim-2] + pre_R[dim-2][dim-2]) >> 2;
    (dst+RIDX(dim-1, dim-1, dim))->green = (pre_G[dim][dim] - pre_G[dim-2][dim] - pre_G[dim][dim-2] + pre_G[dim-2][dim-2]) >> 2;
    (dst+RIDX(dim-1, dim-1, dim))->blue = (pre_B[dim][dim] - pre_B[dim-2][dim] - pre_B[dim][dim-2] + pre_B[dim-2][dim-2]) >> 2;

    //第一行
    for (i = 1; i < dim-1; ++i) {
        // dst[RIDX(0,i,dim)].red = (pre_R[2][i+2] - pre_R[2][i-1]) / 6;
        // dst[RIDX(0,i,dim)].green = (pre_G[2][i+2] - pre_G[2][i-1]) / 6;
        // dst[RIDX(0,i,dim)].blue = (pre_B[2][i+2] - pre_B[2][i-1]) / 6;
        (dst+RIDX(0,i,dim))->red = (pre_R[2][i+2] - pre_R[2][i-1]) / 6;
        (dst+RIDX(0,i,dim))->green = (pre_G[2][i+2] - pre_G[2][i-1]) / 6;
        (dst+RIDX(0,i,dim))->blue = (pre_B[2][i+2] - pre_B[2][i-1]) / 6;

    }

    //最后一行
    for (i = 1; i < dim-1; ++i) {
        // dst[RIDX(dim-1,i,dim)].red = (pre_R[dim][i+2] - pre_R[dim-2][i+2] - pre_R[dim][i-1] + pre_R[dim-2][i-1]) / 6;
        // dst[RIDX(dim-1,i,dim)].green = (pre_G[dim][i+2] - pre_G[dim-2][i+2] - pre_G[dim][i-1] + pre_G[dim-2][i-1]) / 6;
        // dst[RIDX(dim-1,i,dim)].blue = (pre_B[dim][i+2] - pre_B[dim-2][i+2] - pre_B[dim][i-1] + pre_B[dim-2][i-1]) / 6;
        (dst+RIDX(dim-1,i,dim))->red = (pre_R[dim][i+2] - pre_R[dim-2][i+2] - pre_R[dim][i-1] + pre_R[dim-2][i-1]) / 6;
        (dst+RIDX(dim-1,i,dim))->green = (pre_G[dim][i+2] - pre_G[dim-2][i+2] - pre_G[dim][i-1] + pre_G[dim-2][i-1]) / 6;
        (dst+RIDX(dim-1,i,dim))->blue = (pre_B[dim][i+2] - pre_B[dim-2][i+2] - pre_B[dim][i-1] + pre_B[dim-2][i-1]) / 6;
    }

    //第一列
    for (i = 1; i < dim-1; ++i) {
        // dst[RIDX(i,0,dim)].red = (pre_R[i+2][2] - pre_R[i-1][2]) / 6;
        // dst[RIDX(i,0,dim)].green = (pre_G[i+2][2] - pre_G[i-1][2]) / 6;
        // dst[RIDX(i,0,dim)].blue = (pre_B[i+2][2] - pre_B[i-1][2]) / 6;
        (dst+RIDX(i,0,dim))->red = (pre_R[i+2][2] - pre_R[i-1][2]) / 6;
        (dst+RIDX(i,0,dim))->green = (pre_G[i+2][2] - pre_G[i-1][2]) / 6;
        (dst+RIDX(i,0,dim))->blue = (pre_B[i+2][2] - pre_B[i-1][2]) / 6;
    }

    //最后一列
    for (i = 1; i < dim-1; ++i) {
        dst[RIDX(i,dim-1,dim)].red = (pre_R[i+2][dim] - pre_R[i-1][dim] - pre_R[i+2][dim-2] + pre_R[i-1][dim-2]) / 6;
        dst[RIDX(i,dim-1,dim)].green = (pre_G[i+2][dim] - pre_G[i-1][dim] - pre_G[i+2][dim-2] + pre_G[i-1][dim-2]) / 6;
        dst[RIDX(i,dim-1,dim)].blue = (pre_B[i+2][dim] - pre_B[i-1][dim] - pre_B[i+2][dim-2] + pre_B[i-1][dim-2]) / 6;
        // (dst+RIDX(i,dim-1,dim))->red = (pre_R[i+2][dim] - pre_R[i-1][dim] - pre_R[i+2][dim-2] + pre_R[i-1][dim-2]) / 6;
        // (dst+RIDX(i,dim-1,dim))->green = (pre_G[i+2][dim] - pre_G[i-1][dim] - pre_G[i+2][dim-2] + pre_G[i-1][dim-2]) / 6;
        // (dst+RIDX(i,dim-1,dim))->blue = (pre_B[i+2][dim] - pre_B[i-1][dim] - pre_B[i+2][dim-2] + pre_B[i-1][dim-2]) / 6;
    }
    int tmp = 0;
    for (i = 1; i < dim - 1; ++i) {
        tmp = i*dim;
        for (j = 1; j < dim - 1; ++j) {
            // dst[RIDX(i,j,dim)].red = (pre_R[i+2][j+2]-pre_R[i+2][j-1]-pre_R[i-1][j+2]+pre_R[i-1][j-1]) / 9;
            // dst[RIDX(i,j,dim)].green = (pre_G[i+2][j+2]-pre_G[i+2][j-1]-pre_G[i-1][j+2]+pre_G[i-1][j-1]) / 9;
            // dst[RIDX(i,j,dim)].blue = (pre_B[i+2][j+2]-pre_B[i+2][j-1]-pre_B[i-1][j+2]+pre_B[i-1][j-1]) / 9;

            // (dst+RIDX(i,j,dim))->red = (pre_R[i+2][j+2]-pre_R[i+2][j-1]-pre_R[i-1][j+2]+pre_R[i-1][j-1]) / 9;
            // (dst+RIDX(i,j,dim))->green = (pre_G[i+2][j+2]-pre_G[i+2][j-1]-pre_G[i-1][j+2]+pre_G[i-1][j-1]) / 9;
            // (dst+RIDX(i,j,dim))->blue = (pre_B[i+2][j+2]-pre_B[i+2][j-1]-pre_B[i-1][j+2]+pre_B[i-1][j-1]) / 9;

            // dst[i*dim+j].red = (pre_R[i+2][j+2]-pre_R[i+2][j-1]-pre_R[i-1][j+2]+pre_R[i-1][j-1]) / 9;
            // dst[i*dim+j].green = (pre_G[i+2][j+2]-pre_G[i+2][j-1]-pre_G[i-1][j+2]+pre_G[i-1][j-1]) / 9;
            // dst[i*dim+j].blue = (pre_B[i+2][j+2]-pre_B[i+2][j-1]-pre_B[i-1][j+2]+pre_B[i-1][j-1]) / 9;
            dst[tmp+j].red = (pre_R[i+2][j+2]-pre_R[i+2][j-1]-pre_R[i-1][j+2]+pre_R[i-1][j-1]) / 9;
            dst[tmp+j].green = (pre_G[i+2][j+2]-pre_G[i+2][j-1]-pre_G[i-1][j+2]+pre_G[i-1][j-1]) / 9;
            dst[tmp+j].blue = (pre_B[i+2][j+2]-pre_B[i+2][j-1]-pre_B[i-1][j+2]+pre_B[i-1][j-1]) / 9;

            // dst[RIDX(j,i,dim)].red = (pre_R[j+2][i+2]-pre_R[j+2][i-1]-pre_R[j-1][i+2]+pre_R[j-1][i-1]) / 9;
            // dst[RIDX(j,i,dim)].green = (pre_G[j+2][i+2]-pre_G[j+2][i-1]-pre_G[j-1][i+2]+pre_G[j-1][i-1]) / 9;
            // dst[RIDX(j,i,dim)].blue = (pre_B[j+2][i+2]-pre_B[j+2][i-1]-pre_B[j-1][i+2]+pre_B[j-1][i-1]) / 9;
        }
    }

}


/********************************************************************* 
 * register_smooth_functions - 通过为每一个测试函数调用add_smooth_funtion(),
 * 登记所有不同版本的smooth代码到评测程序driver中。
 * 当你运行driver程序时，它将测试并且给出每一个已经登记的测试函数的性能。
 *********************************************************************/

void register_smooth_functions() {
    add_smooth_function(&smooth, smooth_descr);
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    /* ... Register additional test functions here */
}


