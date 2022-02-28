#pragma once
#include "Data.h"

//ˮԴ����
class Source
{
public:

	void SetParmameter(const Parameter* parameter);   //���ò���

	void SetState(const State* state);   //����״̬

	void UpdateState(State* state);    //����״̬

	void calculate();   //������ˮԴ

	Source(double sm = 20.0, double ex = 1.5, double kg = 0.35, double ki = 0.35,
		double r = 0.0, double rs = 0.0, double ri = 0.0, double rg = 0.0,
		double pe = 0.0, double fr = 0.0, double s0 = 0.0, double s = 0.0,
		double m = 0.0, double kid = 0.0, double kgd = 0.0,
		double smm = 0.0, double smmf = 0.0, double smf = 0.0, double au = 0.0,
		double rsd = 0.0, double rid = 0.0, double rgd = 0.0, double fr0 = 0.0,
		int n = 1, double q = 0.0, double kidd = 0.0, double kgdd = 0.0, double dt_ = 0.0);

	~Source();

protected:
private:
	//========ģ�Ͳ���========//

	double SM;   //�������ˮ��ˮ����/mm ������

	double EX;   //�������ˮ��ˮ�������Σ������У�1.0~1.5

	double KG;   //�������ˮ��ˮ��Ե���ˮ���ճ���ϵ��������

	double KI;   //�������ˮ��ˮ������������ճ���ϵ��������

	//========ģ��״̬========//

	double R;	 //�ܾ�������mm

	double RS;   //���澶����mm

	double RI;   //��������mm

	double RG;   //���¾�����mm

	double PE;	 //��������mm

	double FR;	 //��ʱ�β����������

	double S0;   //��ʱ�γ�������ˮ������mm

	double S;    //��ʱ�ε�����ˮ������mm


	double M;    //һ�컮�ֵļ���ʱ����

	double KID;   //�������ˮ��ˮ����������ļ���ʱ�γ���ϵ��������

	double KGD;    //�������ˮ��ˮ��Ե���ˮ�ļ���ʱ�γ���ϵ��������

	double SMM;    //ȫ���򵥵���������ˮ��ˮ������mm

	double SMMF;  //������������һ�������ˮ��ˮ������mm

	double SMF;   //��������ϵ�ƽ������ˮ��ˮ�����mm

	double AU;   //��Ӧƽ����ˮ��������ˮ�S0ֵ��Ӧ�������꣬mm

	double RSD;  //���㲽�����澶����mm

	double RID;  //���㲽����������mm

	double RGD;  //���㲽�����¾�����mm

	double FR0;   //��һʱ�β����������

	int N;       //N Ϊ����ʱ�ηֶ�����ÿһ��Ϊ���㲽��

	double Q;     //Q ��ÿ�����㲽���ڵľ�������mm

	double KIDD;   //�������ˮ��ˮ����������ļ��㲽������ϵ��������

	double KGDD;   //�������ˮ��ˮ��Ե���ˮ�ļ��㲽������ϵ��������

	double dt;   //ģ�ͼ���ʱ�γ�,h
};