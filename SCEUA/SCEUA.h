#pragma once

/******************************************************************************
文件名: SCEUA.h
作者: 卢家波
单位：河海大学水文水资源学院
邮箱：lujiabo@hhu.edu.cn
版本：2021.11 创建 V1.0
版权: MIT
引用格式：卢家波，SCEUA算法C++实现. 南京：河海大学，2021.
         LU Jiabo, Shuffled Complex Evolution in C++. Nanjing:Hohai University, 2021.
参考文献：[1]段青云,SCEUA的原始Fortran代码,1992, https://shxy.hhu.edu.cn/2019/0904/c12296a195177/page.htm
		 [2]L. Shawn Matott改编的C++代码,2009, https://github.com/MESH-Model/MESH_Project_Baker_Creek/blob/7e0a7e588213916deb2b6c11589df0d132d9b310/Model/Ostrich/SCEUA.h
		 [3]Van Hoey S改编的Python代码,2011
		 [4]Mostapha Kalami Heris, Shuffled Complex Evolution in MATLAB (URL: https://yarpiz.com/80/ypea110-shuffled-complex-evolution), Yarpiz, 2015.
		 [5]Duan, Q.Y., Gupta, V.K. & Sorooshian, S. Shuffled complex evolution approach for effective and efficient global minimization. J Optim Theory Appl 76, 501–521 (1993). https://doi.org/10.1007/BF00939380.
		 [6]Duan, Q., Sorooshian, S., and Gupta, V. (1992), Effective and efficient global optimization for conceptual rainfall-runoff models, Water Resour. Res., 28( 4), 1015– 1031, https://doi.org/10.1029/91WR02985.
		 [7]Duan, Q., Sorooshian, S., & Gupta, V. K. (1994). Optimal use of the SCE-UA global optimization method for calibrating watershed models. Journal of Hydrology, 158(3-4), 265-284. https://doi.org/10.1016/0022-1694(94)90057-4.
		 [8]王书功. 水文模型参数估计方法及参数估计不确定性研究[M]. 河南:黄河水利出版社,2010.(https://book.douban.com/subject/5377630/)
		 [9]王书功. 水文模型参数估计方法及参数估计不确定性研究[D]. 北京:中国科学院研究生院,2006.(https://jz.docin.com/p-87849994.html)

==========================================================================
考虑最小化问题

待优化参数信息:
		NOPT = 待优化的参数数量，问题维数，目标函数含有的变量数，简记为n
		XNAME(.) = 变量名
		A(.) = 初始参数集
		BL(.) = 参数下限
		BU(.) = 参数上限

SCE-UA算法控制参数：
		NGS = 初始种群中的复形数量，简记为p
		NPG = 每个复形中的点的个数，简记为m
		NPS = 子复形中的点数，简记为q
		ALPHA = CCE步骤4中 a-e 重复的次数，即复形洗牌前每个复形允许的进化步骤数
		BETA = CCE步骤3-5重复的次数
		NPT = 初始种群中的总点数(NPT=NGS*NPG)，简记为s

收敛检查参数：
		MAXN = 优化终止前允许的最大试验次数
		KSTOP = 在终止优化之前，标准值必须改变给定百分比的洗牌循环数
		PCENTO = 给定数量的随机循环中，标准值必须改变的百分比
		PEPS = 给定数量的随机循环中，函数值最低变化率

标志变量:
		INIFLG = 标记是否将初始点包括在种群中
			= false, 不包括
			= true, 包括
		IPRINT = 用于控制每个洗牌循环后的打印输出的标志
			= false, 打印关于种群最佳点的信息
			= true, 打印关于种群每个点的信息
		IDEFLT = 是否采用默认参数值标志
			= false, 采用
			= true, 不采用

优化结果：
		BESTX(.,.) = 每次洗牌循环的最佳点
		BESTF(.) = 每次循环的最佳点函数值

收敛判据：
		ICALL = 模型调用次数
		TIMEOU = 最近KSTOP次数函数值变化率
		GNRNG = 参数范围的归一化几何平均值
******************************************************************************/

#include <vector>
#include <string>


class SCEUA
{
public:
	SCEUA();
	~SCEUA();

	SCEUA(int argc, char* argv[]);
	//优化函数
	void Optimize();  


private:
	//===================封装SCE-UA算法===================
	void scemain(); 

	//===================输入优化参数===================
	void scein();  

	//===================输出优化结果===================
	void sceout();  
	
	//===================SCE-UA算法===================
	void sceua();  

	//2.[生成样本]在参数空间内生成一个含有NPT个点的初始集合		
	void GenerateSample(std::vector<std::vector<double>>& x, std::vector<double>& xf, int& icall);

	//3.[样本点排序]样本点按函数值升序排列
	void RankPoints(std::vector<std::vector<double>>& x, std::vector<double>& xf);

	//4.[划分复形群体]将样本点划分到ngs个复形中，每个复形含npg个点
	void Partition2Complexes(const int& k, const std::vector<std::vector<double>>& x, const std::vector<double>& xf, std::vector<std::vector<double>>& cx, std::vector<double>& cf);

	//5.[复形演化]根据竞争复形演化算法(CCE)，对每个复形进行独立地演化计算
	void cce(std::vector<std::vector<double>>& cx, std::vector<double>& cf, const std::vector<double>& xnstd, int& icall);

	//6.[复形洗牌]将所有复形内的点放回缓冲区D，按函数值升序排列
	void ShuffleComplexes(const int& k, std::vector<std::vector<double>>& x, std::vector<double>& xf, const std::vector<std::vector<double>>& cx, const std::vector<double>& cf);
		
	//7.[收敛判断]如果满足收敛判据，则算法终止，否则回到步骤4
	void CheckConvergence(const std::vector<std::vector<double>>& x, std::vector<double>& xnstd, const std::vector<double>& bound, const int& nloop, int& icall, double& timeou, double& gnrng);

	//===================小工具===================
	//在可行区域内生成一个新点，重载getpnt函数，分别实现均匀分布和正态分布
	void getpnt(std::vector<double>& snew);

	void getpnt(const std::vector<double>& xi, const std::vector<double>& std, std::vector<double>& snew);

	//按函数值的升序对一组点数组和对应的函数值数组进行排序
	void sort(std::vector< std::vector<double> >& x, std::vector<double>& xf);
		
	//检查试验点是否满足所有约束
	void chkcst(const std::vector<double>& xx, bool& ibound);
	
	//===================目标函数===================
	double functn(const int& nopt, const std::vector<double>& x);

	//1.[前处理]在functn中先把自动生成的参数写入到待率定模型的输入文件中
	void PreProcessing(const std::vector<double>& x);
	
	//2.[运行模型]调用待率定模型计算出结果
	void RunModel();
	
	//3.[后处理]计算出待率定模型计算结果与实测结果的指标，将指标值作为functn函数返回值
	double PostProcessing();
	
	std::vector<double> ReadValues(const std::string& fileName);

	double CalculateNSE(const std::vector<double>& simulatedValues, const std::vector<double>& measuredValues) const;
	
	//===================成员变量===================
	
	//每个变量的参数名、初始值、下限、上限
	int m_nopt;  //待优化的参数数量，问题维数，目标函数含有的变量数，简记为n
	std::vector<std::string> m_xname; //变量名
	std::vector<double> m_a;   //参数初始值
	std::vector<double> m_bl;  //参数下限
	std::vector<double> m_bu;  //参数上限
	
	//SCE控制参数	
	int m_ngs;   //种群中的复形数量，简记为p
	int m_npg;   //每个复形中的点的个数，简记为m
	int m_nps;   //子复形即单纯形中的点数，简记为q
	int m_alpha;   //CCE步骤4中 a-e 重复的次数，即复形洗牌前每个复形允许的进化步骤数
	int m_beta;    //CCE步骤3-5重复的次数	
	int m_npt;   //种群中的总点数(NPT=NGS*NPG)，简记为s

	//收敛判据参数
	int m_maxn;    //优化终止前允许的最大试验次数
	int m_kstop;   //在终止优化之前，计算函数值变化率的洗牌循环数
	double m_pcento;  //给定数量的随机循环中，函数值最低变化率
	double m_peps;    //判断所有参数收敛的标准

	//标志变量
	bool m_iniflg;  //标记是否将初始点包括在种群中
	bool m_iprint;  //用于控制每个洗牌循环后的打印输出的标志
	bool m_ideflt;  //采用默认值标志

	//存储结果
	std::vector< std::vector<double> > m_bestx; //每次洗牌循环的最佳点
	std::vector<double> m_bestf; //每次循环的最佳点函数值

	//存储收敛判据
	std::vector<int> m_icall;  //优化目标试验次数
	std::vector<double> m_timeou;  //最近kstop次洗牌循环中，函数值变化率
	std::vector<double> m_gnrng;  //参数范围的归一化几何平均值

	//存储待率定模型实测值与模拟值
	std::vector<double> measuredValues;
	std::vector<double> simulatedValues;

	//存储SCEUA程序命令行参数
	int argc;  //参数个数
	std::vector<std::string> argv;  //参数名
	std::string filePath; //可执行程序所在目录
};


