// XAJ+SCEUA.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 2022/2/28 实现SCE-UA算法率定新安江模型

#include <iostream>

#include "./SCEUA/SCEUA.h"

//测试SCEUA算法
int main(int argc, char* argv[])
{
	SCEUA optimization(argc, argv);

	clock_t timeBegin = clock();

	optimization.Optimize();

	clock_t timeEnd = clock();

	double duration = static_cast<double>(timeEnd - timeBegin) / CLOCKS_PER_SEC;

	std::cout << std::endl << "优化用时/S: " << duration << std::endl;

	return 0;
}