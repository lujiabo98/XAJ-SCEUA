#pragma once
#include "Data.h"

//������ɢ��ģ��
class Evapotranspiration
{
public:
	void SetParmameter(const Parameter * parameter);   //���ò���

	void SetState(const State* state);   //����״̬

	void UpdateState(State* state);    //����״̬

	void calculate();    //����������ɢ����

	Evapotranspiration(double kc = 0.8, double lm = 80.0, double c = 0.15,
		double wu = 10, double wl = 30, double ep = 0.0,
		double e = 0.0, double eu = 0.0, double el = 0.0, double ed = 0.0, 
		double p = 0.0, double em = 0.0);	//���캯��

	~Evapotranspiration();	 //��������

protected:
private:
	//========ģ�Ͳ���========//

	double KC;   //������ɢ������ϵ��������

	double LM;   //�²�����ˮ����/mm�����У�60~90

	double C;    //�����ɢ������ϵ���������У�0.10~0.20

	//========ģ��״̬========//

	double WU;   //�ϲ�����ˮ������mm

	double WL;   //�²�����ˮ������mm

	double EP;   //��Ԫ��������������mm

	double E;    //�ܵ���ɢ������mm

	double EU;   //�ϲ���ɢ������mm

	double EL;   //�²���ɢ������mm

	double ED;   //�����ɢ������mm


	//========�ⲿ����========//

	double P;    //��������mm

	double EM;   //ˮ����������mm

};