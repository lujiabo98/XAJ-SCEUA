#include "Confluence.h"

void Confluence::SetParmameter(const Parameter* parameter)
{
	CS = parameter->m_CS;

	CI = parameter->m_CI;

	CG = parameter->m_CG;

	CR = parameter->m_CR;

	IM = parameter->m_IM;
}

void Confluence::SetState(const State* state)
{
	RIM = state->m_RIM;

	RS = state->m_RS;

	RI = state->m_RI;

	RG = state->m_RG;

	QS = state->m_QS;

	QI = state->m_QI;

	QG = state->m_QG;

	QU0 = state->m_QU;

	F = state->m_F;

	dt = state->m_dt;
}

void Confluence::UpdateState(State* state)
{
	state->m_QS = QS;

	state->m_QI = QI;

	state->m_QG = QG;

	state->m_QU = QU;

	state->m_QU0 = QU0;
}

void Confluence::calculate()
{
	//��͸ˮ����ϵĲ�������̯����Ԫ������

	RS = RS * (1 - IM);

	RI = RI * (1 - IM);

	RG = RG * (1 - IM);

	//����ϵ������

	M = 24.0 / dt;   //һ�컮�ֵļ���ʱ����

	CSD = pow(CS, 1.0 / M);   //����ʱ���ڵ��澶����ˮ�������ϵ��

	CID = pow(CI, 1.0 / M);   //����ʱ������������ˮ�������ϵ��

	CGD = pow(CG, 1.0 / M);   //����ʱ���ڵ���ˮ��ˮ�������ϵ��

	CRD = pow(CR, 1.0 / M);   //����ʱ���ں�����ˮ����ϵ��

	//�µػ���
	U = F / 3.6 / 24;   //��λת��ϵ��

	QS = CSD * QS + (1 - CSD) * (RS + RIM) * U;   //���澶��������澶��ˮ�⣬
	                                              //��������(CSD)����Ϊ���澶���Ժ�����������QS

	QI = CID * QI + (1 - CID) * RI * U;    //����������������ˮ�⣬
										   //��������(CID)����Ϊ�������Ժ�����������QI

	QG = CGD * QG + (1 - CGD) * RG * U;    //���¾����������ˮ��ˮ�⣬��������ˮ
										   //��ˮ�������(CGD)����Ϊ����ˮ�Ժ�����������QG

	QT = QS + QI + QG;

	//������������������ˮ�ⷨ���ҽ�����Ԫ�����������200km2ʱ�ż����������
	if (F < 200)
	{
		QU = QT;   //��Ԫ������������Һӵ��϶̣���ˮ���˶��ĵ�������ͨ����С
		           //�����ֵ������úϲ��ڵ��澶���͵��¾�����һ���������������ͨ�����Ժ���
	}
	else
	{
		QU = CRD * QU + (1 - CRD) * QT;   //����ˮ�ⷨ
		                                  //ֻ���ڵ�Ԫ��������ϴ����������������临��ʱ
		                                  //�ſ��ǵ�Ԫ����ڵĺ�������
	}
}

Confluence::Confluence(double cs, double ci, double cg, double cr, double im,
	double qs, double qi, double qg, double qt, double qu,
	double rs, double ri, double rg, double rim,
	double qi0, double qg0, double qu0, double f,
	double u, double m, double csd, double cid, double cgd, double crd, double dt_)
{
	CS = cs;

	CI = ci;

	CG = cg;
	
	CR = cr;

	IM = im;

	QS = qs;

	QI = qi;

	QG = qg;

	QT = qt;

	QU = qu;

	RS = rs;

	RI = ri;

	RG = rg;

	RIM = rim;

	QI0 = qi0;

	QG0 = qg0;

	QU0 = qu0;

	F = f;

	U = u;

	M = m;

	CSD = csd;

	CID = cid;

	CGD = cgd;

	CRD = crd;

	dt = dt_;
}

Confluence::~Confluence()
{

}