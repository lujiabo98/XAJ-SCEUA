#include "Watershed.h"
#include <fstream>

void Watershed::ReadFromFile(std::string StrPath)
{
	std::ifstream fin(StrPath + "watershed.txt");

	fin >> m_name
		>> m_area
		>> m_numRainfallStation
		>> m_numEvaporationStation
		>> m_numSubWatershed;

	//��ȡ��Ԫ�������
	m_areaSubWatershed = new double[m_numSubWatershed];

	for (int i = 0; i < m_numSubWatershed; i++)
	{
		fin >> m_areaSubWatershed[i];
	}

	//��ȡ����վ�������
	//������������*����վ�еĶ�ά��̬����������ʼ��Ϊ0
	m_rateRainfallStation = 
		std::vector< std::vector<double> >(m_numSubWatershed, std::vector<double>(m_numRainfallStation, 0.0));

	for (int i = 0; i < m_numSubWatershed; i++)
	{
		for (int j = 0; j < m_numRainfallStation; j++)
		{
			fin >> m_rateRainfallStation[i][j];
		}
	}

	//��ȡ����վ�������
    //������������*����վ�еĶ�ά��̬����������ʼ��Ϊ0
	m_rateEvaporationStation = 
		std::vector< std::vector<double> >(m_numSubWatershed, std::vector<double>(m_numEvaporationStation, 0.0));

	for (int i = 0; i < m_numSubWatershed; i++)
	{
		for (int j = 0; j < m_numEvaporationStation; j++)
		{
			fin >> m_rateEvaporationStation[i][j];
		}
	}

	//��ȡ����վ������
	m_nameRainfallStation = new std::string[m_numRainfallStation];

	for (int j = 0; j < m_numRainfallStation; j++)
	{
		fin >> m_nameRainfallStation[j];
	}

	//��ȡ����վ������
	m_nameEvaporationStation = new std::string[m_numEvaporationStation];

	for (int j = 0; j < m_numEvaporationStation; j++)
	{
		fin >> m_nameEvaporationStation[j];
	}

	fin.close();
}

void Watershed::SetValues(std::string name, double area,
	int numRainfallStation, int numEvaporationStation, int numSubWatershed,
	double* areaSubWatershed, std::vector< std::vector<double> > rateRainfallStation,
	std::vector< std::vector<double> > rateEvaporationStation,
	std::string* nameRainfallStation, std::string* nameEvaporationStation,
	std::vector< std::vector<double> > P, std::vector< std::vector<double> > EM)
{
	m_name = name;

	m_area = area;

	m_numRainfallStation = numRainfallStation;

	m_numEvaporationStation = numEvaporationStation;

	m_numSubWatershed = numSubWatershed;

	m_areaSubWatershed = areaSubWatershed;

	m_rateRainfallStation = rateRainfallStation;

	m_rateEvaporationStation = rateEvaporationStation;

	m_nameRainfallStation = nameRainfallStation;

	m_nameEvaporationStation = nameEvaporationStation;

	m_P = P;

	m_EM = EM;
}

void Watershed::calculate(const IO * io)
{
	int nrows = io->nrows;   //��¼��������ʱ����

	int ncols = m_numSubWatershed;    //��Ԫ�������

	//�������Ԫ������ʱ�ν�������mm
	m_P = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0)); 

	for (int r = 0; r < nrows; r++)
	{
		for (int c = 0; c < ncols; c++)
		{
			for (int i = 0; i < m_numRainfallStation; i++)
			{
				m_P[r][c] += io->m_P[r][i] * m_rateRainfallStation[c][i];   //���������㵥Ԫ��������
			}
		}
	}

	//�������Ԫ������ʱ��ˮ����������mm
	m_EM = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0));

	for (int r = 0; r < nrows; r++)
	{
		for (int c = 0; c < ncols; c++)
		{
			for (int i = 0; i < m_numEvaporationStation; i++)
			{
				m_EM[r][c] += io->m_EM[r][i] * m_rateEvaporationStation[c][i];   //���������㵥Ԫ����ˮ��������
			}
		}
	}
}

double Watershed::GetP(int nt, int nw)
{
	return m_P[nt][nw];
}

double Watershed::GetEM(int nt, int nw)
{
	return m_EM[nt][nw];
}

double Watershed::GetF(int nw)
{
	return m_areaSubWatershed[nw];
}

int Watershed::GetnW()
{
	return m_numSubWatershed;
}

Watershed::Watershed()
{
	m_name = "Ĭ������";

	m_area = 0.0;

	m_numRainfallStation = 0;

	m_numEvaporationStation = 0;

	m_numSubWatershed = 0;

	m_areaSubWatershed = nullptr;

	m_rateRainfallStation = 
		std::vector< std::vector<double> >(m_numSubWatershed, std::vector<double>(m_numRainfallStation, 0.0));

	m_rateEvaporationStation = 
		std::vector< std::vector<double> >(m_numSubWatershed, std::vector<double>(m_numEvaporationStation, 0.0));

	m_nameRainfallStation = nullptr;

	m_nameEvaporationStation = nullptr;

	m_P = std::vector< std::vector<double> >(0, std::vector<double>(m_numSubWatershed, 0.0));

	m_EM = std::vector< std::vector<double> >(0, std::vector<double>(m_numSubWatershed, 0.0));
}

Watershed::~Watershed()
{
	delete[] m_areaSubWatershed;

}


void IO::ReadFromFile(std::string StrPath)
{
	std::ifstream fin(StrPath + "P.txt");

	//���뽵�����ݵ�������������������ʾʱ������������ʾ����վ��
	fin >> nrows >> ncols;

	//������̬��ά���ݣ����ڴ洢������վ��������mm
	m_P = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0));

	//��ȡ������վ���������ݣ�mm
	for (int r = 0; r < nrows; r++)
	{
		for (int c = 0; c < ncols; c++)
		{
			fin >> m_P[r][c];
		}
	}

	fin.close();

	fin.open(StrPath + "EM.txt");

	//�����������ݵ�������������������ʾʱ������������ʾ����վ��
	fin >> nrows >> ncols;

	//������̬��ά���ݣ����ڴ洢������վ��������mm
	m_EM = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0));

	//��ȡ������վ���������ݣ�mm
	for (int r = 0; r < nrows; r++)
	{
		for (int c = 0; c < ncols; c++)
		{
			fin >> m_EM[r][c];
		}
	}

	fin.close();
}

void IO::WriteToFile(std::string StrPath)
{
	//��Q.txt���������ڶ����������̣�û�и��ļ����½�
	std::ofstream fout(StrPath + "Q.txt");

	//�����������
	for (int i = 0; i < nrows; i++)
	{
		fout << m_Q[i] << std::endl;
	}

	fout.close();
}

IO::IO()
{
	m_Q = nullptr;

	nrows = 0;

	ncols = 0;

	m_P = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0));

	m_EM = std::vector< std::vector<double> >(nrows, std::vector<double>(ncols, 0.0));
}

IO::~IO()
{
	delete[] m_Q;
}