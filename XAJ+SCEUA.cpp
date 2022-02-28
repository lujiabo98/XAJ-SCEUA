// XAJ+SCEUA.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 2022/2/28 实现SCE-UA算法率定新安江模型

#include "./SCEUA/SCEUA.h"

//测试SCEUA算法
int main(int argc, char* argv[])
{
	SCEUA optimization(argc, argv);

	optimization.Optimize();

	return 0;
}