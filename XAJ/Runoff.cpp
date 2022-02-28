#include "Runoff.h"

void Runoff::SetParmameter(const Parameter* parameter)
{
	WM = parameter->m_WM;

	B = parameter->m_B;

	IM = parameter->m_IM;

	WUM = parameter->m_UM;

	WLM = parameter->m_LM;

}

void Runoff::SetState(const State* state)
{
	WU = state->m_WU;

	WL = state->m_WL;

	WD = state->m_WD;

	W = state->m_W;

	P = state->m_P;

	EU = state->m_EU;

	EL = state->m_EL;

	ED = state->m_ED;

	EP = state->m_EP;
}

void Runoff::UpdateState(State* state)
{
	state->m_WU = WU;

	state->m_WL = WL;

	state->m_WD = WD;

	state->m_W = W;

	state->m_R = R;

	state->m_PE = PE;

	state->m_RIM = RIM;
}

void Runoff::calculate()
{
	//========�������========//
	WMM = (1 + B) / (1 - IM) * WM;   //��������ˮ�������ֵ��mm

	A = WMM * (1 - pow(1 - W / WM, 1 / (1 + B)));   //��ʼ������ˮ�����ֵ��mm

	PE = P - EP;

	if (PE <= 1e-5)    //������Ϊ������С��1e-5ʱ��ΪС�ڵ���0
	{
		R = 0.0;

		RIM = 0.0;   //���㲻͸ˮ����ϵĲ�����
	}
	else
	{
		if (A + PE <= WMM)
		{
			R = PE + W - WM + WM * pow(1 - (A + PE) / WMM, B + 1);
		}
		else
		{
			R = PE - (WM - W);
		}

		RIM = PE * IM;   //���㲻͸ˮ����ϵĲ�����
	}

	//========������һʱ�γ�������ˮ��========//
	WU = WU + P - EU - R;

	WL = WL - EL;

	WD = WD - ED;

	if (WD < 0)
	{
		WD = 0;   //��ֹ�������ˮ����С��0
	}

	//��������ˮ����������
	if (WU > WUM)
	{
		WL = WL + WU - WUM;

		WU = WUM;
	}

	if (WL > WLM)
	{
		WD = WD + WL - WLM;

		WL = WLM;
	}

	WDM = WM - WUM - WLM;  //�����������ˮ����

	if (WD > WDM)
	{
		WD = WDM;
	}

	//����������ˮ��
	W = WU + WL + WD;
}


Runoff::Runoff(double wm, double b, double im, 
	double wum, double wlm, double wdm,
	double r, double rim, double w,
	double wu, double wl, double wd,
	double wmm, double a, 
	double eu, double el, double ed, 
	double ep, double pe, double p)
{
	WM = wm;
	B = b;
	IM = im;
	WUM = wum;
	WLM = wlm;
	WDM = wdm;

	R = r;
	RIM = rim;
	W = w;
	WU = wu;
	WL = wl;
	WD = wd;
	WMM = wmm;
	A = a;
	EU = eu;
	EL = el;
	ED = ed;
	EP = ep;
	PE = pe;

	P = p;

}

Runoff::~Runoff()
{

}