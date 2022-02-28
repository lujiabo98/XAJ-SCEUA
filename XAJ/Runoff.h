#pragma once
#include "Data.h"

//��������ģ��
class Runoff
{
public:

	void SetParmameter(const Parameter* parameter);   //���ò���

	void SetState(const State* state);   //����״̬

	void UpdateState(State* state);    //����״̬

	void calculate();	//�������㲢����������ˮ��

	Runoff(double wm = 120.0, double b = 0.2, double im = 0.01, 
		double wum = 20.0, double wlm = 80.0, double wdm = 20.0,
		double r = 0.0, double rim = 0.0, double w = 0.0, 
		double wu = 10.0, double wl = 30.0, double wd = 20.0,
		double wmm = 0.0, double a = 0.0, 
		double eu = 0.0, double el = 0.0, double ed = 0.0,
		double ep = 0.0, double pe =0.0, double p =0.0);

	~Runoff();

protected:
private:
	//========ģ�Ͳ���========//

	double WM;   //����ƽ������ˮ����/mm�������У�120~200

	double B;    //����ˮ��ˮ�������߷��Σ������У�0.1~0.4

	double IM;   //��͸ˮ���ռȫ��������ı�����������

	double WUM;   //�ϲ�����ˮ����/mm�����У�10~50

	double WLM;   //�²�����ˮ����/mm�����У�60~90

	double WDM;   //�������ˮ����/mm������WM - WUM - WLM�������ڲ���

	//========ģ��״̬========//

	double R;	//�ܾ�������mm

	double RIM;   //��͸ˮ����ϵĲ�������mm

	double W;	//����ƽ����ʼ������ˮ����mm

	double WU;   //�ϲ�����ˮ������mm

	double WL;   //�²�����ˮ������mm

	double WD;   //�������ˮ������mm

	double WMM;   //��������ˮ�������ֵ��mm

	double A;    //��ʼ������ˮ�����ֵ��mm

	double EU;   //�ϲ���ɢ������mm

	double EL;   //�²���ɢ������mm

	double ED;   //�����ɢ������mm

	double EP;   //��������������mm

	double PE;	//��������mm


	//========�ⲿ����========//

	double P;    //��������mm

};