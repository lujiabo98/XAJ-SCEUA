#pragma once

/******************************************************************************
�ļ���: SCEUA.h
����: ¬�Ҳ�
��λ���Ӻ���ѧˮ��ˮ��ԴѧԺ
���䣺lujiabo@hhu.edu.cn
�汾��2021.11 ���� V1.0
��Ȩ: MIT
���ø�ʽ��¬�Ҳ���SCEUA�㷨C++ʵ��. �Ͼ����Ӻ���ѧ��2021.
         LU Jiabo, Shuffled Complex Evolution in C++. Nanjing:Hohai University, 2021.
�ο����ף�[1]������,SCEUA��ԭʼFortran����,1992, https://shxy.hhu.edu.cn/2019/0904/c12296a195177/page.htm
		 [2]L. Shawn Matott�ı��C++����,2009, https://github.com/MESH-Model/MESH_Project_Baker_Creek/blob/7e0a7e588213916deb2b6c11589df0d132d9b310/Model/Ostrich/SCEUA.h
		 [3]Van Hoey S�ı��Python����,2011
		 [4]Mostapha Kalami Heris, Shuffled Complex Evolution in MATLAB (URL: https://yarpiz.com/80/ypea110-shuffled-complex-evolution), Yarpiz, 2015.
		 [5]Duan, Q.Y., Gupta, V.K. & Sorooshian, S. Shuffled complex evolution approach for effective and efficient global minimization. J Optim Theory Appl 76, 501�C521 (1993). https://doi.org/10.1007/BF00939380.
		 [6]Duan, Q., Sorooshian, S., and Gupta, V. (1992), Effective and efficient global optimization for conceptual rainfall-runoff models, Water Resour. Res., 28( 4), 1015�C 1031, https://doi.org/10.1029/91WR02985.
		 [7]Duan, Q., Sorooshian, S., & Gupta, V. K. (1994). Optimal use of the SCE-UA global optimization method for calibrating watershed models. Journal of Hydrology, 158(3-4), 265-284. https://doi.org/10.1016/0022-1694(94)90057-4.
		 [8]���鹦. ˮ��ģ�Ͳ������Ʒ������������Ʋ�ȷ�����о�[M]. ����:�ƺ�ˮ��������,2010.(https://book.douban.com/subject/5377630/)
		 [9]���鹦. ˮ��ģ�Ͳ������Ʒ������������Ʋ�ȷ�����о�[D]. ����:�й���ѧԺ�о���Ժ,2006.(https://jz.docin.com/p-87849994.html)

==========================================================================
������С������

���Ż�������Ϣ:
		NOPT = ���Ż��Ĳ�������������ά����Ŀ�꺯�����еı����������Ϊn
		XNAME(.) = ������
		A(.) = ��ʼ������
		BL(.) = ��������
		BU(.) = ��������

SCE-UA�㷨���Ʋ�����
		NGS = ��ʼ��Ⱥ�еĸ������������Ϊp
		NPG = ÿ�������еĵ�ĸ��������Ϊm
		NPS = �Ӹ����еĵ��������Ϊq
		ALPHA = CCE����4�� a-e �ظ��Ĵ�����������ϴ��ǰÿ����������Ľ���������
		BETA = CCE����3-5�ظ��Ĵ���
		NPT = ��ʼ��Ⱥ�е��ܵ���(NPT=NGS*NPG)�����Ϊs

������������
		MAXN = �Ż���ֹǰ���������������
		KSTOP = ����ֹ�Ż�֮ǰ����׼ֵ����ı�����ٷֱȵ�ϴ��ѭ����
		PCENTO = �������������ѭ���У���׼ֵ����ı�İٷֱ�
		PEPS = �������������ѭ���У�����ֵ��ͱ仯��

��־����:
		INIFLG = ����Ƿ񽫳�ʼ���������Ⱥ��
			= false, ������
			= true, ����
		IPRINT = ���ڿ���ÿ��ϴ��ѭ����Ĵ�ӡ����ı�־
			= false, ��ӡ������Ⱥ��ѵ����Ϣ
			= true, ��ӡ������Ⱥÿ�������Ϣ
		IDEFLT = �Ƿ����Ĭ�ϲ���ֵ��־
			= false, ����
			= true, ������

�Ż������
		BESTX(.,.) = ÿ��ϴ��ѭ������ѵ�
		BESTF(.) = ÿ��ѭ������ѵ㺯��ֵ

�����оݣ�
		ICALL = ģ�͵��ô���
		TIMEOU = ���KSTOP��������ֵ�仯��
		GNRNG = ������Χ�Ĺ�һ������ƽ��ֵ
******************************************************************************/

#include <vector>
#include <string>


class SCEUA
{
public:
	SCEUA();
	~SCEUA();

	SCEUA(int argc, char* argv[]);
	//�Ż�����
	void Optimize();  


private:
	//===================��װSCE-UA�㷨===================
	void scemain(); 

	//===================�����Ż�����===================
	void scein();  

	//===================����Ż����===================
	void sceout();  
	
	//===================SCE-UA�㷨===================
	void sceua();  

	//2.[��������]�ڲ����ռ�������һ������NPT����ĳ�ʼ����		
	void GenerateSample(std::vector<std::vector<double>>& x, std::vector<double>& xf, int& icall);

	//3.[����������]�����㰴����ֵ��������
	void RankPoints(std::vector<std::vector<double>>& x, std::vector<double>& xf);

	//4.[���ָ���Ⱥ��]�������㻮�ֵ�ngs�������У�ÿ�����κ�npg����
	void Partition2Complexes(const int& k, const std::vector<std::vector<double>>& x, const std::vector<double>& xf, std::vector<std::vector<double>>& cx, std::vector<double>& cf);

	//5.[�����ݻ�]���ݾ��������ݻ��㷨(CCE)����ÿ�����ν��ж������ݻ�����
	void cce(std::vector<std::vector<double>>& cx, std::vector<double>& cf, const std::vector<double>& xnstd, int& icall);

	//6.[����ϴ��]�����и����ڵĵ�Żػ�����D��������ֵ��������
	void ShuffleComplexes(const int& k, std::vector<std::vector<double>>& x, std::vector<double>& xf, const std::vector<std::vector<double>>& cx, const std::vector<double>& cf);
		
	//7.[�����ж�]������������оݣ����㷨��ֹ������ص�����4
	void CheckConvergence(const std::vector<std::vector<double>>& x, std::vector<double>& xnstd, const std::vector<double>& bound, const int& nloop, int& icall, double& timeou, double& gnrng);

	//===================С����===================
	//�ڿ�������������һ���µ㣬����getpnt�������ֱ�ʵ�־��ȷֲ�����̬�ֲ�
	void getpnt(std::vector<double>& snew);

	void getpnt(const std::vector<double>& xi, const std::vector<double>& std, std::vector<double>& snew);

	//������ֵ�������һ�������Ͷ�Ӧ�ĺ���ֵ�����������
	void sort(std::vector< std::vector<double> >& x, std::vector<double>& xf);
		
	//���������Ƿ���������Լ��
	void chkcst(const std::vector<double>& xx, bool& ibound);
	
	//===================Ŀ�꺯��===================
	double functn(const int& nopt, const std::vector<double>& x);

	//1.[ǰ����]��functn���Ȱ��Զ����ɵĲ���д�뵽���ʶ�ģ�͵������ļ���
	void PreProcessing(const std::vector<double>& x);
	
	//2.[����ģ��]���ô��ʶ�ģ�ͼ�������
	void RunModel();
	
	//3.[����]��������ʶ�ģ�ͼ�������ʵ������ָ�꣬��ָ��ֵ��Ϊfunctn��������ֵ
	double PostProcessing();
	
	std::vector<double> ReadValues(const std::string& fileName);

	double CalculateNSE(const std::vector<double>& simulatedValues, const std::vector<double>& measuredValues) const;
	
	//===================��Ա����===================
	
	//ÿ�������Ĳ���������ʼֵ�����ޡ�����
	int m_nopt;  //���Ż��Ĳ�������������ά����Ŀ�꺯�����еı����������Ϊn
	std::vector<std::string> m_xname; //������
	std::vector<double> m_a;   //������ʼֵ
	std::vector<double> m_bl;  //��������
	std::vector<double> m_bu;  //��������
	
	//SCE���Ʋ���	
	int m_ngs;   //��Ⱥ�еĸ������������Ϊp
	int m_npg;   //ÿ�������еĵ�ĸ��������Ϊm
	int m_nps;   //�Ӹ��μ��������еĵ��������Ϊq
	int m_alpha;   //CCE����4�� a-e �ظ��Ĵ�����������ϴ��ǰÿ����������Ľ���������
	int m_beta;    //CCE����3-5�ظ��Ĵ���	
	int m_npt;   //��Ⱥ�е��ܵ���(NPT=NGS*NPG)�����Ϊs

	//�����оݲ���
	int m_maxn;    //�Ż���ֹǰ���������������
	int m_kstop;   //����ֹ�Ż�֮ǰ�����㺯��ֵ�仯�ʵ�ϴ��ѭ����
	double m_pcento;  //�������������ѭ���У�����ֵ��ͱ仯��
	double m_peps;    //�ж����в��������ı�׼

	//��־����
	bool m_iniflg;  //����Ƿ񽫳�ʼ���������Ⱥ��
	bool m_iprint;  //���ڿ���ÿ��ϴ��ѭ����Ĵ�ӡ����ı�־
	bool m_ideflt;  //����Ĭ��ֵ��־

	//�洢���
	std::vector< std::vector<double> > m_bestx; //ÿ��ϴ��ѭ������ѵ�
	std::vector<double> m_bestf; //ÿ��ѭ������ѵ㺯��ֵ

	//�洢�����о�
	std::vector<int> m_icall;  //�Ż�Ŀ���������
	std::vector<double> m_timeou;  //���kstop��ϴ��ѭ���У�����ֵ�仯��
	std::vector<double> m_gnrng;  //������Χ�Ĺ�һ������ƽ��ֵ

	//�洢���ʶ�ģ��ʵ��ֵ��ģ��ֵ
	std::vector<double> measuredValues;
	std::vector<double> simulatedValues;

	//�洢SCEUA���������в���
	int argc;  //��������
	std::vector<std::string> argv;  //������
	std::string filePath; //��ִ�г�������Ŀ¼
};


