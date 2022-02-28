/******************************************************************************
文件名: SCEUA.cpp
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

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <numbers> //使用圆周率常数

#include "SCEUA.h"
#include "./XAJ/XAJ.h"


/********************************************************
Optimize()
使用 SCEUA 最小化目标函数
*********************************************************/
void SCEUA::Optimize()
{
    scemain(); //执行SCE程序
}

/**********************************************************
 scemain()
 执行SCE程序，调用scein()、sceua()、sceout()函数
***********************************************************/
void SCEUA::scemain()
{
	scein();   //设置优化参数

	sceua();   //SCE-UA算法

	sceout();  //输出优化结果
}

/************************************************************
 scein()
 该函数读取用于全局优化的SCE方法的输入变量

 输入参数列表：
        StrPath = 文件路径，默认为空，即可执行程序所在路径
************************************************************/
void SCEUA::scein()
{
    //初始化 I/O 变量
    std::ifstream fin(filePath + "scein.txt");
    
    //读取参数个数
    fin >> m_nopt;

    //创建参数名、初始值、上下限数组
    m_xname = std::vector<std::string>(m_nopt);
    m_a = std::vector<double>(m_nopt);
    m_bl = std::vector<double>(m_nopt);
    m_bu = std::vector<double>(m_nopt);

    //依次读取每个变量的参数名、初始值、下限、上限
    for (int i = 0; i < m_nopt; i++)
    {
        fin >> m_xname.at(i)  //读取变量名
            >> m_a.at(i)  //读取参数初始值
            >> m_bl.at(i)  //读取参数下限
            >> m_bu.at(i); //读取参数上限
    }

    //读取默认值标识
    fin >> m_ideflt;

    //如果IDEFLT等于1，接着读入SCE控制参数、收敛判据参数、标志变量
    if (m_ideflt == true)
    {
		fin >> m_ngs
			>> m_npg
			>> m_nps
            >> m_alpha
            >> m_beta

            >> m_maxn
			>> m_kstop
			>> m_pcento
            >> m_peps
 
            >> m_iniflg 
            >> m_iprint;
    }
    //否则将SCE控制参数设置为推荐值即默认值
    else
    {
        m_ngs = 5;

        m_npg = 2 * m_nopt + 1;

        m_nps = m_nopt + 1;

        m_alpha = 1;

        m_beta = 2 * m_nopt + 1;

        m_maxn = 5000;

        m_kstop = 10;

        m_pcento = 0.1;
        
        m_peps = 0.001;


        m_iniflg = true;

        m_iprint = false;
    }

    //计算初始种群中的总点数
    m_npt = m_ngs * m_npg;

	//检查SCE控制参数是否有效

    fin.close();

    //读取实测值
    measuredValues = ReadValues(filePath + "observe.txt");
    
}

/**********************************************************
 sceout()
 输出SCE-UA算法运行结果
 
 输入参数列表：
		StrPath = 文件路径，默认为空，即可执行程序所在路径
 
 局部变量列表：
		fout = 文件输出句柄
***********************************************************/
void SCEUA::sceout()
{   
    //=================1.控制台输出=================
	//打印初始点及其函数值
	std::cout << "===================参数初始值及上下限===================" << std::endl;

    std::cout << "参数名：";

	for (int i = 0; i < m_nopt; i++)
	{
		std::cout << m_xname[i] << "  ";
	}

	std::cout << std::endl;

    std::cout << "参数初始值：";

	for (int i = 0; i < m_nopt; i++)
	{
		std::cout << m_a[i] << "  ";
	}

	std::cout << std::endl;

    std::cout << "参数下限：";

	for (int i = 0; i < m_nopt; i++)
	{
		std::cout << m_bl[i] << "  ";
	}

	std::cout << std::endl;

    std::cout << "参数上限：";

	for (int i = 0; i < m_nopt; i++)
	{
		std::cout << m_bu[i] << "  ";
	}

	std::cout << std::endl;

	std::cout << "初始参数函数值：" << functn(m_nopt, m_a) << std::endl;

	//打印优化结果
    std::cout << "===================SCE-UA搜索的结果===================" << std::endl;

    std::cout << "最优点:"  << std::endl;

	for (int i = 0; i < m_nopt; i++)
	{
		std::cout << m_xname[i] << "  ";
	}

    std::cout << std::endl;

    for (auto i : m_bestx.back())
    {
        std::cout << i << "  ";
    }

	std::cout << std::endl;

    std::cout << "最优函数值:" << m_bestf.back() << std::endl;

    std::cout << m_bestf.size() << "次洗牌演化的最优点函数值：" << std::endl;
	
    for (auto i : m_bestf)
	{
		std::cout << i << "  ";
	}

    //=================2.文件输出=================
    std::ofstream fout(filePath + "sceout.txt");  //创建输出文件
    
	//打印初始点及其函数值
    fout << "===================参数初始值及上下限===================" << std::endl;

    fout << "参数名：";

	for (int i = 0; i < m_nopt; i++)
	{
        fout << m_xname[i] << "  ";
	}

    fout << std::endl;

    fout << "参数初始值：";

	for (int i = 0; i < m_nopt; i++)
	{
        fout << m_a[i] << "  ";
	}

    fout << std::endl;

    fout << "参数下限：";

	for (int i = 0; i < m_nopt; i++)
	{
        fout << m_bl[i] << "  ";
	}

    fout << std::endl;

    fout << "参数上限：";

	for (int i = 0; i < m_nopt; i++)
	{
        fout << m_bu[i] << "  ";
	}

    fout << std::endl;

    fout << "初始参数函数值：" << functn(m_nopt, m_a) << std::endl;

	//打印优化结果
    fout << "===================SCE-UA搜索的结果===================" << std::endl;

    fout << "最优点:" << std::endl;

	for (int i = 0; i < m_nopt; i++)
	{
        fout << m_xname[i] << "  ";
	}

    fout << std::endl;

	for (auto i : m_bestx.back())
	{
        fout << i << "  ";
	}

    fout << std::endl;

    fout << "最优函数值:" << m_bestf.back() << std::endl;

    fout << m_bestf.size() << "次洗牌演化的最优点函数值：" << std::endl;

	for (auto i : m_bestf)
	{
        fout << i << "  ";
	}

    fout.close();  //关闭输出文件
}


/***********************************************************
 sceua()
 该算法对目标函数最小化问题进行全局优化搜索

 局部变量列表：
        X(.,.) = 种群中点的坐标
		XF(.) = X(.,.)的函数值
		CX(.,.) = 复形中点的坐标
		CF(.) = CX(.,.)的函数值
		XNSTD(.) = 种群中参数的标准差
		BOUND(.) = 优化变量的约束
		ICALL = 模型调用次数
	    TIMEOU = 最近kstop次数函数值变化率
	    GNRNG = 参数范围的归一化几何平均值
		NLOOP = 主循环次数	
************************************************************/
void SCEUA::sceua()
{   
	std::cout << "==================================================" << std::endl
		<< "                  进入SCE-UA全局搜索                " << std::endl
		<< "==================================================" << std::endl;

    //===========1.[初始化]定义局部变量及数组，并分配内存===========  
       
    std::vector< std::vector<double> > x(m_npt, std::vector<double>(m_nopt));  //种群中点的坐标
    
    std::vector<double> xf(m_npt);  //种群中点的坐标X(.,.)的函数值
	    
	std::vector< std::vector<double> > cx(m_npg, std::vector<double>(m_nopt)); //复形中点的坐标
    
    std::vector<double> cf(m_npg);  //复形中点的坐标CX(.,.)的函数值
       
    std::vector<double> xnstd(m_nopt, 0.1);   //种群中参数的标准差，初始化为0.1           
    
    std::vector<double> bound(m_nopt);   //优化变量的约束范围    
	
    int icall = 0;  //模型调用次数

    double timeou = 10000.0;  //最近kstop次洗牌循环中函数值变化率，默认为一较大数值

    double gnrng = 10000.0; //参数范围的归一化几何平均值，默认为一较大数值

	int nloop = 0;  //主循环次数   
   
	for (int j = 0; j < m_nopt; j++)
	{
		bound[j] = m_bu[j] - m_bl[j];  //计算待优化参数约束范围
	}
	
    //===========2.[生成样本]在参数空间内生成一个含有NPT个点的初始集合===========
    GenerateSample(x, xf, icall);

    //===========3.[样本点排序]样本点按函数值升序排列===========
    RankPoints(x, xf);

    std::cout << "目标函数调用次数" << "  " <<  "函数值变化率" << "  " << "参数变化范围" << std::endl;

    //===========7.[收敛判断]不满足收敛判据，返回步骤4===========
    while ((icall < m_maxn) && (timeou > m_pcento) && (gnrng > m_peps))
    {
        nloop += 1;  //主循环次数计数

		//对每个复形进行独立演化计算，并放回到缓冲区D中
		for (int k = 0; k < m_ngs; k++)
		{
			//===========4.[划分复形群体]将样本点划分到p个复形中，每个复形含m个点===========
			Partition2Complexes(k, x, xf, cx, cf);

			//===========5.[复形演化]根据CCE对每个复形独立演化计算===========
			cce(cx, cf, xnstd, icall);

			//===========6.[复形洗牌]将所有复形内的点放回缓冲区D===========
			ShuffleComplexes(k, x, xf, cx, cf);
		}

		//===========3.[样本点排序]将所有复形内的点按函数值升序排列===========
		RankPoints(x, xf);

		//===========7.[收敛判断]如果满足收敛判据，则算法终止，否则回到步骤4===========
        CheckConvergence(x, xnstd, bound, nloop, icall, timeou, gnrng);

        std::cout << "      " << icall << "            " << timeou << "        " << gnrng << std::endl;
    }
       
}


/************************************************************
 GenerateSample()
 2.[生成样本]在参数空间内生成一个含有NPT个点的初始集合

 输入参数列表：   
	X(.,.) = 种群中点的坐标
	XF(.) = X(.,.)的函数值
    ICALL = 模型调用次数

 局部变量列表：
    XX(.) = 种群中点X(.,.)中单点的坐标
************************************************************/
void SCEUA::GenerateSample(        
    std::vector< std::vector<double> >& x,
    std::vector<double>& xf,
    int& icall
)
{
	std::vector<double> xx(m_nopt);   //种群中点X(.,.)中单点的坐标

	//生成在参数空间均匀分布的NPT个随机点
	for (int i = 0; i < m_npt; i++)
	{
		//在参数空间内按照均匀分布随机生成点
		getpnt(xx);

		//将生成的点赋值给种群
		x.at(i) = xx;
	}

	//如果INIFLG等于true，设置X(0,.)为初始点A(.)，将初始点包含在初始集合中
	if (m_iniflg == true)
	{
		x.front() = m_a;
	}

    //计算NPT个随机点对应的函数值
    for (int i = 0; i < m_npt; i++)
    {
        xf.at(i) = functn(m_nopt, x.at(i));  //计算生成点的函数值

        icall += 1; //模型调用次数增加1
    }
}

/************************************************************
 RankPoints()
 3.[样本点排序]样本点按函数值升序排列

 输入参数列表：
    X(.,.) = 种群中点的坐标
	XF(.) = X(.,.)的函数值
************************************************************/
void SCEUA::RankPoints(	
	std::vector< std::vector<double> >& x,
    std::vector<double>& xf
)
{
	sort(x, xf);  //根据目标函数值将样本点按从小到大排序

	//记录最好的点

	m_bestx.push_back(x.front());  //最好点的坐标，即第一个点

	m_bestf.push_back(xf.front());  //最好点的函数值，即第一个点的函数值
}


/************************************************************
 Partition2Complexes()
 4.[划分复形群体]将样本点划分到ngs个复形中，每个复形含npg个点

 输入参数列表：
    K = 复形的序号，K = 1,2,···,ngs
    X(.,.) = 种群中点的坐标
	XF(.) = X(.,.)的函数值
	CX(.,.) = 复形中点的坐标
	CF(.) = 复形中点的坐标CX(.,.)的函数值

 局部变量列表：
    INDEX = 按照规则（k + ngs * j）从种群中提取点到复形中的下标
************************************************************/
void SCEUA::Partition2Complexes
(
    const int& k,  
    const std::vector< std::vector<double> >& x,
    const std::vector<double>& xf,
    std::vector< std::vector<double> >& cx,
    std::vector<double>& cf
)
{
    int index = 0;   //种群的下标

    //遍历复形中的每个点，为其从种群中赋值
    for (int j = 0; j < m_npg; j++)
    {
        index = k + m_ngs * j;   //计算下标

        cx.at(j) = x.at(index);   //划分样本点到复形

        cf.at(j) = xf.at(index);  //划分对应的函数值
    }
}

/************************************************************
 cce()
 5.[复形演化]根据CCE对每个复形独立演化计算

 输入参数列表：
	CX(.,.) = 复形中点的坐标
	CF(.) = 复形中点的坐标CX(.,.)的函数值
    XNSTD(.) = 种群中参数的归一化标准差
    ICALL = 模型调用次数

 局部变量列表：
	s(.,.) = 当前单纯形中点的坐标
    sf(.) = 当前单纯形中点的坐标S(.,.)的函数值
    sb(.) = 单纯形的最佳点
    sw(.) = 单纯形的最坏点WO（WORST POINT，简称WO）
    ce(.) = 单纯形排除最坏点（WO）的形心
    snew(.) = 从单纯形生成的新点
    par(.) = 除最差点以外的参数取值集合
    lcs(.) = 索引定位S(.,.)在 CX(.,.)中的位置
	wts = 复形中每个点的权重,weights
	vals = 经过均匀随机扰动的复形中每个点的权重
	valsWithIndices = 权重及下标索引
    gen = 随机数引擎
    u = 产生[0,1)均匀分布随机数
    fw = 最差点的函数值
    fnew = 新点SNEW的函数值    
    ibound = 违反约束指标
            = false,  没有违反
            = true,  违反
************************************************************/
void SCEUA::cce
(
	std::vector< std::vector<double> >& cx,
	std::vector<double>& cf,
    const std::vector<double>& xnstd,  
    int& icall
)
{
    //===========1.[初始化]定义局部变量===========    
    
    //定义数组并分配内存	
	std::vector< std::vector<double> > s(m_nps, std::vector<double>(m_nopt)); //当前单纯形中点的坐标
    
    std::vector<double> sf(m_nps);  //当前单纯形中点的坐标S(.,.)的函数值
    
    std::vector<double> sb(m_nopt); //子复形（即单纯形）的最佳点，其函数值最小
    
    std::vector<double> sw(m_nopt); //子复形（即单纯形）的最差点，其函数值最大
    
    std::vector<double> ce(m_nopt); //子复形（即单纯形）排除最坏点（WO）的形心
    
    std::vector<double> snew(m_nopt); //从子复形（即单纯形）生成的新点
    
    std::vector<double> par(m_nps - 1); //除最差点以外的参数取值集合
    
    std::vector<int> lcs(m_nps);   //索引定位S(.,.)在 CX(.,.)中的位置
    
    std::vector<int> wts(m_npg); //用于存储复形中每个点的权重
    
    std::vector<double> vals;  //用于存储经过随机扰动的复形中每个点的权重
    
    std::vector<std::pair<int, double>> valsWithIndices; //存储下标索引及其权重
	
	auto gen = std::mt19937{ std::random_device{}() }; //随机数引擎采用设备熵值保证随机性
	
	std::uniform_real_distribution<double> u(0.0, 1.0); //产生[0,1)均匀分布随机数

    double fw = 0.0;   //最差点的函数值
    
    double fnew = 0.0;   //新点SNEW的函数值
    
    bool ibound = true;  //违反约束指标，默认违反

    //===========2.[分配权重]根据三角概率分布，为复形内每个点分配概率===========
    //pi=2(m+1-i)/m(m+1),i=1,2,...,m 由于2/m(m+1)相同，故可以省略
    //各下标索引的权重为整数 m+1-i
    //权重赋值，这里的 i=0,1,...,m-1
    //例：npg=5即m=5时，权重依次为5,4,3,2,1
    for (int i = 0; i < m_npg; i++)
    {
        wts.at(i) = m_npg - i;
    }    
        
    //===========6.[迭代]重复步骤 3-5 beta次===========
    for (int ibeta = 0; ibeta < m_beta; ibeta++)
    {
        //===========3.[选择父辈群体]从复形中加权不放回抽取q个点===========
        //并将q个点及其函数值存储于缓冲区B中

        //清空权重数组为排序做准备
        vals.clear();
        valsWithIndices.clear();

        //把权重值作为指数，均匀分布随机数作为底
        //相当于给权重增加了均匀随机扰动
        for (auto iter : wts)
        {
            vals.push_back(std::pow(u(gen), 1.0 / iter));
        }

        //构成<下标索引,权重>对的数组
        for (int i = 0; i < m_npg; i++)
        {
            valsWithIndices.emplace_back(i, vals[i]);
        }

        //按照扰动后的权重值从大到小排序，并扩展到下标索引
        std::sort(valsWithIndices.begin(), valsWithIndices.end(), [](auto x, auto y) {return x.second > y.second; });

        //按照扰动后的权重值从大到小抽样q个点，记录下标索引到lcs数组中
        for (int i = 0; i < m_nps; i++)
        {
            lcs.at(i) = valsWithIndices[i].first;
        }

        //对抽取到的q个点的下标按照从小到大升序排列
        std::sort(lcs.begin(), lcs.end());

        //由抽取到的q个点构成单纯形，记录单纯形中点的坐标及函数值
        for (int i = 0; i < m_nps; i++)
        {
            s.at(i) = cx.at(lcs.at(i));  //从复形中取单纯形中点的坐标

            sf.at(i) = cf.at(lcs.at(i));  //从复形函数值数组中取单纯形中点的函数值
        }


        //===========4.[产生下一代群体]按照下山单纯形法的6步对q个点演化计算,重复alpha次===========
        /* --------------------------------------------------
          确定子复形S的最坏点WO
          计算剩余点的形心CE（CENTROID，简称CE）
          计算步长，最坏点WO和质心CE之间的向量
          确定最坏的函数值FW（The worst function value，简称FW）
        --------------------------------------------------- */

        for (int ialpha = 0; ialpha < m_alpha; ialpha++)
        {

            sb = s.front(); //复制最佳点坐标

            sw = s.back();  //复制最差点坐标

            fw = sf.back();  //最坏的函数值是单纯形S的最后一个点，函数值最大

            //遍历当前单纯形中点的坐标
            for (int j = 0; j < m_nopt; j++)
            {
                for (int i = 0; i < m_nps - 1; i++)
                {
                    par[i] = s[i][j];   //将除最差点以外的参数取值赋值给par数组
                }

                //a.[形心]计算除最差点以外nps-1个点的平均值，即为形心ce
                ce[j] = std::accumulate(par.begin(), par.end(), 0.0) / double(m_nps - 1);
            }

            //计算新点SNEW
            //b.[反射]首先尝试反射步骤，计算最差点的反射点r=2g-uq,
            //这里符号记为snew=2ce-sw
            for (int j = 0; j < m_nopt; j++)
            {
                snew[j] = 2 * ce[j] - sw[j];
            }

            //检查SNEW是否违反约束,即是否在可行域内
            chkcst(snew, ibound);

            /* ---------------------------------------------------
              c.[突变]如果SNEW在界外，
              按正态分布在可行区域内随机选择一个点，
              子复形的最佳点作为种群的均值
              标准差作为STD
            ----------------------------------------------------- */
            if (ibound == true)
            {
                //产生一个正态分布的满足约束的新点
                getpnt(sb, xnstd, snew);
            }

            //计算SNEW点的函数值FNEW
            fnew = functn(m_nopt, snew);
            icall += 1; //模型调用次数增加1

            //将FNEW与最差函数值FW进行比较
            //如果FNEW小于FW，接受新点SNEW并返回
            if (fnew < fw)
            {
                //用新点替换最坏点
                s.back() = snew;

                sf.back() = fnew;
            }
            //d.[收缩]否则计算新点snew=(g+sw)/2和fnew,若fnew<fw,以snew代替最差点sw
            else
            {
                for (int j = 0; j < m_nopt; j++)
                {
                    snew[j] = (ce[j] + sw[j]) / 2.0;  //即为snew = (ce+sw)/2
                }
                //计算新点SNEW的函数值FNEW
                fnew = functn(m_nopt, snew);
                icall += 1; //模型调用次数增加1

                //如果收缩成功，fnew小于fw，接受新点snew，并替换sw
                if (fnew < fw)
                {
                    //用新点替换最坏点
                    s.back() = snew;

                    sf.back() = fnew;
                }
                //e.[突变]此时反射和收缩都失败，
                //则根据正态分布随机生成新点snew,以snew代替sw
                //子复形的最佳点作为种群的平均值并且标准偏差作为STD
                else
                {
                    //产生一个正态分布的新点
                    getpnt(sb, xnstd, snew);

                    //计算新点SNEW的函数值FNEW
                    fnew = functn(m_nopt, snew);
                    icall += 1; //模型调用次数增加1

                    //用新点替换最坏点
                    s.back() = snew;

                    sf.back() = fnew;
                }
            }

            //按照函数值对单纯形中的点从小到大排序
            sort(s, sf);
        }

        //===========5.[以子代取代父代群体]将q个点按照位置放回复形，升序排列===========
        for (int i = 0; i < m_nps; i++)
        {
            cx.at(lcs.at(i)) = s.at(i);

            cf.at(lcs.at(i)) = sf.at(i);
        }

        sort(cx, cf);  //按照函数值升序排列
    }    
}  //cce 结束

/************************************************************
 ShuffleComplexes()
 6.[复形洗牌]将所有复形内的点放回缓冲区D

 输入参数列表：
	K = 复形的序号，K = 1,2,...,ngs
	X(.,.) = 种群中点的坐标
	XF(.) = X(.,.)的函数值
	CX(.,.) = 复形中点的坐标
	CF(.) = 复形中点的坐标CX(.,.)的函数值

 局部变量列表：
	INDEX = 按照规则（k + ngs * j）从复形中提取点到种群中的下标
************************************************************/
void SCEUA::ShuffleComplexes
(
	const int& k,
	std::vector< std::vector<double> >& x,
	std::vector<double>& xf,
	const std::vector< std::vector<double> >& cx,
	const std::vector<double>& cf
)
{
	int index = 0;   //种群的下标

    //遍历复形中的每个点，将其值赋给种群
	for (int j = 0; j < m_npg; j++)
	{
		index = k + m_ngs * j;   //计算下标

		x.at(index) = cx.at(j);   //样本点从复形赋值到种群

		xf.at(index) = cf.at(j);  //函数值从复形赋值到种群
	}
}

/**********************************************************
 CheckConvergence()
 7.[收敛判断]如果满足收敛判据，则算法终止，否则回到步骤4

 输入参数列表：
    X(.,.) = 种群中所有点的坐标
    XNSTD(.) = 种群中参数的归一化标准差
    BOUND(.) = 优化变量的约束
    NLOOP = 主循环次数
    ICALL = 模型调用次数   
    TIMEOU = 最近kstop次数函数值变化率
    GNRNG = 参数范围的归一化几何平均值

 局部变量列表：
    denomi = 分母denominator
    xmax(.) = 各个维度参数值的最大值
    xmin(.) = 各个维度参数值的最小值
    xmean(.) = 各个维度参数值的平均值
    par(.) = 某个参数的所有取值集合
    delta = 极小值，防止log真值为0
    parsum = 某个参数值总和
    parsum2 = 某个参数的平方和
    gsum = 所有参数的几何平均值
***********************************************************/
void SCEUA::CheckConvergence
(
	const std::vector< std::vector<double> >& x, //种群中所有点的坐标
	std::vector<double>& xnstd,  //种群中参数的归一化标准差
	const std::vector<double>& bound,  //优化变量的约束    
    const int& nloop,
    int& icall,
    double& timeou,
    double& gnrng  
)
{
    //==========1.检查是否超过了函数评估次数的最大值==========
	if (icall > m_maxn)
	{
		//搜索终止
		std::cout << "经过" << nloop << "次洗牌演化，"
            << "优化搜索已经终止，因为超过了最大试验次数" << m_maxn << "次的限制" << std::endl;
	}

    //==========2.检查无函数值改进的连续循环次数==========

    double denomi = 0.0;    //分母denominator

    if (nloop >= m_kstop)
    {
        //计算分母denominator，最近kstop次数的函数值绝对值的平均值
		denomi = std::fabs(std::accumulate(m_bestf.end() - m_kstop, m_bestf.end(), 0.0)) / double(m_kstop);
		
        //最近kstop次数函数值变化率
        timeou = fabs(m_bestf.at(nloop - 1) - m_bestf.at(nloop - m_kstop)) / denomi;
        
        //如果变化率小于pcento，则已经收敛
		if (timeou < m_pcento)
		{
			std::cout << "最佳点在最近" << m_kstop << 
                "次循环中函数值变化率小于阈值" << m_pcento << std::endl;

            std::cout << "经过" << nloop << "次洗牌演化，"
                << "基于目标函数标准已经实现了收敛！" << std::endl;
		}
    }
    
    //==========3.检查种群是否聚集到一个足够小的空间==========

	std::vector<double> xmax(m_nopt);  //各个维度参数值的最大值
	std::vector<double> xmin(m_nopt);  //各个维度参数值的最小值
	std::vector<double> xmean(m_nopt); //各个维度参数值的平均值
	std::vector<double> par(m_npt); //某个参数的所有取值集合

	const double delta = 1.0e-20;   //极小值，防止log真值为0

	//计算参数值的最大值、最小值和标准差
	double parsum = 0.0; //某个参数值总和
	double parsum2 = 0.0; //某个参数的平方和
	double gsum = 0.0;  //所有参数的几何平均值

	for (int k = 0; k < m_nopt; k++)
	{

		for (int i = 0; i < m_npt; i++)
		{
			par[i] = x[i][k];  //存储某个参数的所有取值
		}

		//计算参数值的最大值、最小值、总和、平均值、标准差、归一化标准差
		xmax[k] = *std::max_element(par.begin(), par.end()); //最大值

		xmin[k] = *std::min_element(par.begin(), par.end()); //最小值

		parsum = std::accumulate(par.begin(), par.end(), 0.0); //总和

		xmean[k] = parsum / double(m_npt);  //平均值

        parsum2 = std::inner_product(par.begin(), par.end(), par.begin(), 0.0);  //计算平方和sum(xi^2)

		//根据方差的性质3：Var(x) = E(x^2) - E(x)^2，期望平方内减外
		xnstd[k] = sqrt(parsum2 / double(m_npt) - pow(xmean[k], 2)); //标准差
		xnstd[k] = xnstd[k] / bound[k];  //归一化标准差

		//对每个参数的几何平均值求和
		gsum += log(delta + (xmax[k] - xmin[k]) / bound[k]);
	}

	gnrng = exp(gsum / double(m_nopt));   //种群参数变化范围的相对值

	if (gnrng < m_peps)
	{
        std::cout << "经过" << nloop << "次洗牌演化，"
            << "种群已经收敛到一个预先指定的小参数空间" << std::endl;
	}
  
    //==========4.存储每次洗牌循环的收敛判据==========

	m_icall.push_back(icall);  //记录优化目标试验次数

	m_timeou.push_back(timeou);  //记录最近kstop次洗牌循环中，函数值变化率

	m_gnrng.push_back(gnrng);  //记录参数范围的归一化几何平均值
}


/**********************************************************
 getpnt()
 在可行区域内根据均匀分布生成一个新点

 输入参数列表：
	SNEW(.) = 产生的新点

 局部变量列表：
	 gen = 随机引擎
	 ibound = 违反约束指标
				= false,  没有违反
				= true,  违反
	 u = 均匀分布
***********************************************************/
void SCEUA::getpnt(std::vector<double>& snew) 
{
	
    //随机数引擎采用设备熵值保证随机性
	auto gen = std::mt19937{ std::random_device{}() };

    //违反约束指标，默认违反
    bool ibound = true; 

    while (ibound == true)
    {
        //根据要求产生新点
		for (int j = 0; j < m_nopt; j++)
		{
			//产生均匀分布随机数
			std::uniform_real_distribution<double> u(m_bl[j], m_bu[j]);

			snew[j] = u(gen);
		}

        //显性和隐性检查
        chkcst(snew, ibound);
    }
}

/**********************************************************
 getpnt()
 在可行区域内根据正态分布生成一个新点

 输入参数列表：
    SNEW(.) = 产生的新点
    std(.) = 概率分布的标准差
    xi(.) = 焦点

 局部变量列表：
     gen = 随机引擎
     ibound = 违反约束指标
                = false,  没有违反
                = true,  违反
	 n = 正态分布
***********************************************************/
void SCEUA::getpnt
(
    const std::vector<double>& xi,   //焦点
	const std::vector<double>& std,  //概率分布的标准差
	std::vector<double>& snew  //产生的新点
)
{	
    //随机数引擎采用设备熵值保证随机性
	auto gen = std::mt19937{ std::random_device{}() };

	//显性和隐性检查
	bool ibound = true; //违反约束指标，默认违反

	while (ibound == true)
	{
		//根据要求产生新点
		for (int j = 0; j < m_nopt; j++)
		{

			//产生正态分布随机数
            std::normal_distribution<double> n{ xi[j], std[j] };

			snew[j] = n(gen);

		}

		//显性和隐性检查
		chkcst(snew, ibound);
	}
}



/************************************************************
 sort()
 按函数值的升序对一组点数组和对应的函数值数组进行排序
 
 输入参数列表：
	 X(.,.) = 种群中所有点的坐标
	 XF(.) = X(.,.)的函数值

 局部变量列表：
    XXF(.) = 存储坐标和函数值对的数组
*************************************************************/
void SCEUA::sort
(
    std::vector< std::vector<double> >& x, 
    std::vector<double>& xf  
)
{
    //将点的坐标与函数值对应上，以便于扩展排序
    std::vector< std::pair<std::vector<double>, double> > xxf;
	
    //构造向量xxf，存储坐标和函数值对
    for (size_t i = 0; i < xf.size(); i++)
    {
        xxf.emplace_back(x[i], xf[i]);
    }

    //匿名函数定义比较规则，用于pair数据结构按照函数值升序排列
    std::sort(xxf.begin(), xxf.end(), [](auto a, auto b) { return a.second < b.second; });
    
    //将排序后的xxf赋值给点数组x和函数值数组xf
    for (size_t i = 0; i < xxf.size(); i++)
    {
        x[i] = xxf[i].first;  //取出点坐标

        xf[i] = xxf[i].second;   //取出函数值
    }
}

/************************************************************
 chkcst()
 检查试验点是否满足所有约束

 输入参数列表：
	XX(.) = X 中单点的坐标
	IBOUND = 违反约束指标
			= false,  没有违反
			= true,  违反
 
 局部变量列表：
        i = 数组 x, bl 和 bu 的第i个变量
*************************************************************/
void SCEUA::chkcst
(
    const std::vector<double>& xx, //X 中单点的坐标
    bool& ibound  //违反约束指标
)
{
	ibound = true;   //默认违反约束

	//检查是否违反了显式边界约束
	for (int i = 0; i < m_nopt; i++)
	{
        if ((xx[i] < m_bl[i]) || (xx[i] > m_bu[i]))
		{
            // 至少违反了一个约束
            return;  //终止函数运行
		}
	}
   
	// 检查是否违反了隐式约束
    // 此处用于填写隐式约束
	// 此函数没有隐式约束
	


	// 通过上述检查，表明没有违反约束
	ibound = false;   
}

/************************************************************
 functn()
 计算目标函数值

 输入参数列表：
    NOPT = 待优化的参数数量，问题维数，目标函数含有的变量数，简记为n
    X(.) = 单点的坐标，即各参数构成的数组
*************************************************************/
double SCEUA::functn(const int& nopt, const std::vector<double>& x)
{
    

    //1.[前处理]在functn中先把自动生成的参数写入到新安江模型的输入文件中
    PreProcessing(x);

    //2.[运行模型]再调用新安江模型计算出结果
    RunModel();

    //3.[后处理]最后计算出NSE，将1-NSE作为functn函数返回值
    double returnVal = PostProcessing();

    return returnVal;
}

void SCEUA::PreProcessing(const std::vector<double>& x)
{
	//初始化 I/O 变量
    std::ifstream fin(filePath + "parameter.tpl");

    std::ofstream fout(filePath + "parameter.txt");

    //读取参数
    std::string parameter = "";

    while (fin >> parameter)
    {
        //查找是否在参数名数组中
		auto iter = std::find(m_xname.begin(), m_xname.end(), parameter);

		if (iter != m_xname.end())
		{
			//在参数名数组中，将该参数名对应值写入parameter.txt中
            //获取该参数名索引
            int index = static_cast<int>( iter - m_xname.begin() );
            fout << x.at(index) << std::endl;
		}
		else
		{
            //不在参数名数组中，即不是待率定参数，直接写入模型输入文件即可
            fout << parameter << std::endl;
		}
    }

    fout.close();
    fin.close();
}

void SCEUA::RunModel()
{
    XAJ(filePath);
}

double SCEUA::PostProcessing()
{

    simulatedValues = ReadValues(filePath + "Q.txt");

    double NSE = CalculateNSE(simulatedValues, measuredValues);

    return 1 - NSE;
}

std::vector<double> SCEUA::ReadValues(const std::string& fileName)
{
	//初始化 I/O 变量
	std::ifstream fin(fileName);

	std::vector<double> values;

    double val = 0.0;
	while (fin >> val)
	{
		values.push_back(val);
	}

    return values;
}

double SCEUA::CalculateNSE(const std::vector<double>& simulatedValues, const std::vector<double>& measuredValues) const
{
    double measuredValuesSum = std::accumulate(measuredValues.begin(), measuredValues.end(), 0.0);
    
    double measuredValuesAvg = measuredValuesSum / measuredValues.size();

    double numerator = 0.0;

    double denominator = 0.0;

    for (double val : measuredValues)
    {
        denominator += pow(val - measuredValuesAvg, 2);
    }
    
    for (int index = 0; index < measuredValues.size(); ++index)
    {
        numerator += pow(simulatedValues.at(index) - measuredValues.at(index), 2);
    }

    double NSE = 1 - numerator / denominator;

    return NSE;
}

//构造函数
SCEUA::SCEUA()
{
    m_nopt = 5;
    m_xname = std::vector<std::string>{ "x1", "x2", "x3", "x4", "x5" };
    m_a = std::vector<double>{ 0, 0, 0, 0, 0 };
    m_bl = std::vector<double>{ -2, -2, -2, -2, -2 };
    m_bu = std::vector<double>{ 2, 2, 2, 2, 2 };

	m_ngs = 5;
	m_npg = 2 * m_nopt + 1;
	m_nps = m_nopt + 1;
	m_alpha = 1;
	m_beta = 2 * m_nopt + 1;
    m_npt = m_ngs * m_npg;

	m_maxn = 3000;
	m_kstop = 10;
	m_pcento = 0.1;
	m_peps = 0.001;

	m_iniflg = true;
	m_iprint = false;
    m_ideflt = false;

    filePath = "";
}

//析构函数
SCEUA::~SCEUA()
{
}

SCEUA::SCEUA(int argc_, char* argv_[]): SCEUA()
{

    argc = argc_;

    for (int i = 0; i < argc_; ++i)
    {
        argv.emplace_back( static_cast<std::string>(argv_[i]) );
    }

	std::string Path = argv[0];  //获取可执行文件所在路径
	filePath = Path.substr(0, Path.rfind("\\") + 1);   //获取可执行文件所在目录
}
