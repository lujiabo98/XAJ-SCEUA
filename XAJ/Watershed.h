#pragma once
#include <string>
#include <vector>

class IO    //���ı��е��뽵����������ݣ�����������̵��ı���
{
public:

	void ReadFromFile(std::string StrPath = "");   //���ı��е������������������

	void WriteToFile(std::string StrPath = "");    //���������ڶ����������̵��ı���

	IO();

	~IO();

public:

	std::vector< std::vector<double> > m_P;  //������վ��ʱ�ν�������mm

	std::vector< std::vector<double> > m_EM;  //������վ��ʱ��ˮ����������mm

	double* m_Q;    //������ڶ����������̣�m3/s

	int nrows;     //��������

	int ncols;     //��������

};

//����ֿ�
class Watershed
{
public:
	void ReadFromFile(std::string StrPath = "");
	
	void SetValues(std::string name, double area, 
		int numRainfallStation, int numEvaporationStation, int numSubWatershed,
		double* areaSubWatershed, std::vector< std::vector<double> > rateRainfallStation,
		std::vector< std::vector<double> > rateEvaporationStation,
		std::string* nameRainfallStation, std::string* nameEvaporationStation,
		std::vector< std::vector<double> > P, std::vector< std::vector<double> > EM);

	void calculate(const IO * io);   //�������Ԫ����������ˮ��������

	double GetP(int nt, int nw);   //�õ�ntʱ�̣���nw��Ԫ����Ľ�������mm

	double GetEM(int nt, int nw);   //�õ�ntʱ�̣���nw��Ԫ�����ˮ����������mm

	double GetF(int nw);    //�õ���nw��Ԫ����������km2

	int GetnW();   //�õ���Ԫ�������

	Watershed();

	~Watershed();

protected:
private:
	std::string m_name;	//��������

	double m_area;	//������������km2

	int m_numRainfallStation;		//����վ����

	int m_numEvaporationStation;    //����վ����

	int m_numSubWatershed;	//��Ԫ�������

	double* m_areaSubWatershed;	//��Ԫ���������km2

	std::vector< std::vector<double> > m_rateRainfallStation;  //����Ԫ�����Ӧ������վ��������̩ɭ����η���

	std::vector< std::vector<double> > m_rateEvaporationStation;  //����Ԫ�����Ӧ������վ��������̩ɭ����η���

	std::string* m_nameRainfallStation;	//����վ������

	std::string* m_nameEvaporationStation;  //����վ������

	std::vector< std::vector<double> > m_P;  //����Ԫ������ʱ�ν�������mm

	std::vector< std::vector<double> > m_EM;  //����Ԫ������ʱ��ˮ����������mm

};

