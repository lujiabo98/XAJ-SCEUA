#include "Source.h"

void Source::SetParmameter(const Parameter* parameter)
{
	SM = parameter->m_SM;

	EX = parameter->m_EX;

	KG = parameter->m_KG;

	KI = parameter->m_KI;
}

void Source::SetState(const State* state)
{
	R = state->m_R;

	PE = state->m_PE;

	FR = state->m_FR;

	S0 = state->m_S0;

	dt = state->m_dt;
}

void Source::UpdateState(State* state)
{
	state->m_RS = RS;

	state->m_RI = RI;

	state->m_RG = RG;

	state->m_S0 = S;

	state->m_FR = FR;

}

void Source::calculate()
{
	//����ϵ������
	M = 24.0 / dt;   //һ�컮�ֵļ���ʱ����

	KID = (1 - pow(1 - (KI + KG), 1.0 / M)) / (1 + KG / KI);   //�������ˮ��ˮ����������ļ���ʱ�γ���ϵ��������

	KGD = KID * KG / KI;   //�������ˮ��ˮ��Ե���ˮ�ļ���ʱ�γ���ϵ��������

	//����ˮԴ[4]
	if (PE <= 1e-5)
	{
		//������С�ڵ���0ʱ,������Ϊ������С��1e-5ʱ��ΪС�ڵ���0
		RS = 0;

		RI = KID * S0 * FR;   /*��������С�ڵ���0ʱ����������ˮ��ˮ���е�ˮ
							  �������������Ϊ��һʱ�εı�������*/
		RG = KGD * S0 * FR;

		S = S0 * (1 - KID - KGD);   //������һʱ�γ�������ˮ����

	}
	else
	{
		//����������0ʱ
		SMM = SM * (1 + EX);   //ȫ���򵥵���������ˮ��ˮ������mm

		FR0 = FR;   //��һʱ�β����������

		FR = R / PE;   //���㱾ʱ�β����������

		if (FR > 1)    //���FR����С���������������1���������ǿ����Ϊ1
		{
			FR = 1;
		}

		S = S0 * FR0 / FR;

		N = int(PE / 5.0) + 1;   //N Ϊ����ʱ�ηֶ�����ÿһ��Ϊ���㲽��

		Q = PE / N;    //Q ��ÿ�����㲽���ڵľ�������mm
		                      //R ���ܾ�������PE�ǵ�Ԫ��������ȡ�
				              //��R��5mm�֣�ʵ�ʲ������ǰ�PE��5mm�֣���Ϊ�˼�С��������仯�Ĳ����
				              //����ˮ��ˮ��ֻ�����ڲ�������ϣ���׿�Ϊ�������FR����ȻFR����ʱ��仯�ġ�
				              //������R����ˮ�⼴�ڲ�������ϲ���PE�ľ����
				              //FR = f / F = R / PE 

		KIDD = (1 - pow(1 - (KID + KGD), 1.0 / N)) / (1 + KGD / KID);  //�������ˮ��ˮ����������ļ��㲽������ϵ��������

		KGDD = KIDD * KGD / KID;  //�������ˮ��ˮ��Ե���ˮ�ļ��㲽������ϵ��������

		//�Ѹ�ʱ�ε�RS��RI��RG��0�����ں����ۼӼ��㲽���ڵ�RSD��RID��RGD
		RS = 0.0;

		RI = 0.0;

		RG = 0.0;

		//���������������һ�������ˮ��ˮ���� SMMF
		if (EX == 0.0)
		{
			SMMF = SMM;   //EX����0ʱ����������ˮ��ˮ�����ֲ�����
		}
		else
		{
			//�ٶ�SMMF��������FR��ȫ���������������ˮ��ˮ����SMM��Ϊ�����߷ֲ�
			//��SMMFӦ������ʽ����
			SMMF = (1 - pow(1 - FR, 1.0 / EX)) * SMM;
		}

		SMF = SMMF / (1.0 + EX);

		//��ÿ������ʱ�ε�������R�ֳ�N�Σ�����������㲽���ڵ�RSD��RID��RGD�����ۼӵõ�RS��RI��RG
		for (int i = 1; i <= N; i++)
		{
			if (S > SMF)
			{
				S = SMF;
			}

			AU = SMMF * (1 - pow(1 - S / SMF, 1.0 / (1 + EX)));

			if (Q + AU <= 0)
			{
				RSD = 0;

				RID = 0;

				RGD = 0;

				S = 0;
			}
			else
			{
				if (Q + AU >= SMMF)
				{
					//���㲽���ڵľ�������������ˮ��ˮ���ʹ������ˮ��ˮ��������������󵥵�����ˮ��ˮ��
					RSD = (Q + S - SMF) * FR;

					RID = SMF * KIDD * FR;

					RGD = SMF * KGDD * FR;

					S = SMF * (1 - KIDD - KGDD);
				}
				else
				{
					//����ˮ��ˮ��δ���������������󵥵�����ˮ��ˮ��
					RSD = (S + Q - SMF + SMF * pow(1 - (Q + AU) / SMMF, 1 + EX)) * FR;

					RID = (S * FR + Q * FR - RSD) * KIDD;

					RGD = (S * FR + Q * FR - RSD) * KGDD;

					S = S + Q - (RSD + RID + RGD) / FR;
				}
			}

			//�ۼӼ���ʱ���ڵĵ��澶�����������͵��¾���
			RS = RS + RSD;

			RI = RI + RID;

			RG = RG + RGD;
		}
	}
}


Source::Source(double sm, double ex, double kg, double ki,
	double r, double rs, double ri, double rg,
	double pe, double fr, double s0, double s,
	double m, double kid, double kgd,
	double smm, double smmf, double smf, double au,
	double rsd, double rid, double rgd, double fr0,
	int n, double q, double kidd, double kgdd, double dt_)
{
	SM = sm;

	EX = ex;

	KG = kg;

	KI = ki;

	R = r;

	RS = rs;

	RI = ri;

	RG = rg;

	PE = pe;

	FR = fr;

	S0 = s0;

	S = s;

	M = m;

	KID = kid;

	KGD = kgd;

	SMM = smm;

	SMMF = smmf;

	SMF = smf;

	AU = au;

	RSD = rsd;

	RID = rid;

	RGD = rgd;

	FR0 = fr0;

	N = n;

	Q = q;

	KIDD = kidd;

	KGDD = kgdd;

	dt = dt_;
}

Source::~Source()
{

}