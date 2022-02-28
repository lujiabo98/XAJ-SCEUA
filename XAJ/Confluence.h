#pragma once
#include "Data.h"

//��Ԫ��������������µػ����ͺ�������������������ˮ�ⷨ
class Confluence
{
public:
	void SetParmameter(const Parameter* parameter);   //���ò���

	void SetState(const State* state);   //����״̬

	void UpdateState(State* state);    //����״̬

	void calculate();   //�µػ����ͺ�������

	Confluence(double cs = 0.1, double ci = 0.6, double cg = 0.95, double cr = 0.2, double im = 0.01,
		double qs = 0.0, double qi = 0.0, double qg = 0.0, double qt = 0.0, double qu = 0.0,
		double rs = 0.0, double ri = 0.0, double rg = 0.0, double rim = 0.0,
		double qi0 = 0.0, double qg0 = 0.0, double qu0 = 0.0, double f = 0.0,
		double u = 0.0, double m = 0.0, double csd = 0.0, 
		double cid = 0.0, double cgd = 0.0, double crd = 0.0, double dt_ = 24);

	~Confluence();
protected:
private:
	//========ģ�Ͳ���========//

	double CS;   //���澶������ϵ��������

	double CI;   //����������ϵ��������

	double CG;   //����ˮ����ϵ��������

	double CR;   //������ˮ����ϵ��������

	double IM;   //��͸ˮ���ռȫ��������ı�����������

	//========ģ��״̬========//

	double QS;   //��Ԫ������澶����m3/s

	double QI;   //��Ԫ������������m3/s

	double QG;   //��Ԫ������¾�����m3/s

	double QT;   //��Ԫ�������������
				 //�����뵥Ԫ����ĵ��澶�����������͵��¾���֮�ͣ���m3/s

	double QU;   //��Ԫ�������������m3/s

	double RS;   //���澶������mm

	double RI;   //��������������mm

	double RG;   //���¾�������mm

	double RIM;   //��͸ˮ����ϵĲ�������mm

	double QI0;  //QI(t-1)��ǰһʱ����������m3/s

	double QG0;  //QG(t-1)��ǰһʱ�̵��¾�����m3/s

	double QU0;  //QU(t-1)��ǰһʱ�̵�Ԫ�������������m3/s

	double F;   //��Ԫ���������km2

	double U;   //��λת��ϵ��

	double M;   //һ�컮�ֵļ���ʱ����

	double CSD;   //����ʱ���ڵ��澶����ˮ�������ϵ��

	double CID;   //����ʱ������������ˮ�������ϵ��

	double CGD;   //����ʱ���ڵ���ˮ��ˮ�������ϵ��

	double CRD;   //����ʱ���ں�����ˮ����ϵ��

	double dt;   //ģ�ͼ���ʱ�γ�,h
};