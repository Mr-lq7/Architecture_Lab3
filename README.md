# Architecture

**体系结构实验**

指导书网址：https://hitsz-lab.gitee.io/computer-arch/lab1/part1/

本次实验实际上是csapp的程序性能优化实验perflab

**实验三PartA，程序性能优化，对图像旋转进行优化
=======


**实验三PartB，程序性能优化，对图像平滑进行优化**
=======


1.注意：kernels.c是平滑原始版本的.c文件

2.kernels_V1.c->kernels_v2.c->kernels_V3.c->kernels_V4.c是依次对图像平滑进行优化的文件

（注意：这里的文件都已经对图像旋转进行优化了，所谓的这些版本只是对图像平滑而言的）

最终结果：图像旋转加速比：2.5，图像平滑加速比：4.6

下载perlab包，在linux下

make driver

./driver 即可

如果有需要可以改变config.h中的值用来校正baseline.

注意：

这里默认的环境是ubuntu16.04 32位

如果你的是ubuntu18.04 64位的话

需要下载包，命令为：sudo apt install libc6-dev-i386

声明：实验PPT均为哈工大深圳仇老师制作