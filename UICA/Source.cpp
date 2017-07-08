#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
#include<iomanip>
#include<sstream>
using namespace std;

#define NumofLight 10        //�O���Ӽ�
#define NumofInstruction 16 //�G�׫��O�Ӽ�
#define NumofTd 5          //�ؼ��I�Ӽ�
#define OtherLight 0        //�B�~�������Ӷ�
//float OtherLightT[OtherLight]={0};     //�B�~�������ӫ�

#define BalanceRatio 5 //�h��Կ����šA��O�t�Z5%��

float the_best = 224.22;

int RunTime = 50;  //�]�h�֦�
#define ImpossiblePowerConsume 9999    //�ӫפ��X��A�������i�઺���j��q���ӡA�Ϩ�e���Q�ӫצX��̨��N

#define NumofCountry 100000    //��NumofCountry�ոѧ@�����(��a��)
#define NumofEmpire 10     //�Ұ��
#define VelocityLimit 3
#define InitialW 0.4             //�D�ʫY�ơA�Ѥj��p
#define FinalW 0.4
#define InitialCOne 2          //������C�Y�ơA�Ѥj��p
#define FinalCOne 2
#define InitialCTwo  2       //�s�骺C�Y�ơA�Ѥp��j
#define FinalCTwo 2


float I[NumofInstruction][2];   //(�ӯ�,�G��)
float Td[NumofTd];
float KParameter[NumofTd][NumofLight + OtherLight];

class Country  //�C�Ӳɤl���򥻳]�w�A�ܰʽե����O����
{
private:
	float PBest[NumofLight + 1],  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,����̧C�ӯ�)
		GBest[NumofLight + 1];  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,�s��̧C�ӯ�)
	int Position[NumofLight];  //�O������A��U���O�����ե����O <----> Position   //�䤤�O���Ƥ]�Osolution�����סA�ĴX���O��ܲĴX��
	int Velocity[NumofLight];             //�ե����O�W�ɩΤU�� <----> Velocity
	float W;//�����D�ʥΪ�weighting�A��ȶV�j�A�V�A�X����j�M�A��ĳ��0.9~0.4
	float C1, C2;//�Y�ơA�վ㰾�V�U��A�ΰ��V�s��
public:
	Country(int* PPointer, float Power, int* VPointer);
	void SetW(float NewW){ W = NewW; }
	void SetC(float NewC1, float NewC2){ C1 = NewC1; C2 = NewC2; }
	float* GetPBest(){ return &PBest[0]; }
	float GetFitness(){ return PBest[NumofLight]; }
	void SetGBest(float* GBestPointer)
	{
		for (int i = 0; i<(NumofLight + 1); i++)
		{
			GBest[i] = *(GBestPointer + i);                                                                 //cout<<"Gbest["<<i<<"]:"<<GBest[i]<<",";
		}                                                                                               //system("pause");              
	}
	void UpdatePosition()
	{
		for (int i = 0; i<NumofLight; i++)//i�O����
		{                                                                                               //cout<<"��P"<<i<<":"<<Position[i]<<",";
			Position[i] = Position[i] + Velocity[i];                                                            //cout<<"V:"<<Velocity[i]<<endl;
			if (Position[i]>NumofInstruction - 1)             //�ե����O���W�U��
				Position[i] = NumofInstruction - 1;
			if (Position[i]<0)
				Position[i] = 0;                                                                           //cout<<"��P"<<i<<":"<<Position[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void UpdateVelocity()
	{
		for (int i = 0; i<NumofLight; i++)//i�O����
		{                                                                                               //cout<<"��V"<<i<<":"<<Velocity[i]<<",";
			Velocity[i] = W*Velocity[i] + C1*(PBest[i] - Position[i]) + C2*(GBest[i] - Position[i]);         //cout<<"��V"<<i<<":"<<Velocity[i]<<endl;
			if (Velocity[i]>VelocityLimit)                   //�t�פW�U��
				Velocity[i] = VelocityLimit;
			if (Velocity[i]<(-VelocityLimit))
				Velocity[i] = (-VelocityLimit);                                                            //cout<<"��V"<<i<<":"<<Velocity[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void CheckPosition()
	{
		float Tpso = 0;  //�ھڥثe�ɤl����m�Һ�X�Ӫ�T(�ӫ�)��
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<NumofLight; j++)
			{
				Tpso += KParameter[i][j] * I[Position[j]][1];
			}
			if (Tpso<Td[i])
			{
				break;
			}
			Tpso = 0;
			if (i == (NumofTd - 1))
			{
				ComputingFitness();
			}
		}
	}
	void ComputingFitness()
	{
		float Fitness = 0;
		for (int i = 0; i<NumofLight; i++)
		{                                                                                            //cout<<i<<":"<<Position[i]<<",";
			Fitness += I[Position[i]][0];
		}                                                                                            //cout<<endl;system("pause");
		if (Fitness<PBest[NumofLight] || PBest[NumofLight] <= 0)
		{                                                                                            //cout<<Fitness<<"<"<<PBest[NumofLight]<<",";
			for (int j = 0; j<NumofLight; j++)
			{
				PBest[j] = Position[j];                                                               ///cout<<PBest[j];
			}
			PBest[NumofLight] = Fitness;                                                              //cout<<"New PBest:"<<PBest[NumofLight]<<endl;
		}
	}
};

Country::Country(int* PPointer, float Power, int* VPointer)
{
	for (int i = 0; i<NumofLight; i++)
	{
		Position[i] = *(PPointer + i);                                                //cout<<Position[i]<<",";                                                         
		Velocity[i] = *(VPointer + i);
		PBest[i] = Position[i];
	}                                                                             //cout<<endl<<endl;
	PBest[NumofLight] = Power;                                                       //cout<<PBest[NumofLight]<<endl;
	W = InitialW;
	C1 = InitialCOne;
	C2 = InitialCTwo;
}
//Country���O����

int RandomNumber(int MinValue, int MaxValue);
int CheckT(int MatrixRow, int i[]);
float Empire();
void quickSort(Country* EmpireCountry[], int left, int right);

int main()
{
	char FileInput[500];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofInstruction
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				stringstream(FileInput) >> I[i][j];
				//sscanf(FileInput, "%f", &I[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofTd
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			stringstream(FileInput) >> Td[i];
			//sscanf(FileInput, "%f", &Td[i]);        //�r����Ʀr
		}
	}
	fr.close();

	//KParameter
	string InputProcess;
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput));//�x�}���B�e
		//fr.getline(FileInput, sizeof(FileInput), ',');//�x�}�e
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput));
			int j = 0;
			/**
			for (int j = 0; j<NumofLight; j++)
			{
			fr.getline(FileInput, sizeof(FileInput), ',');
			stringstream(FileInput) >> KParameter[i][j];
			//sscanf(FileInput, "%f", &KParameter[i][j]);        //�r����Ʀr
			cout << KParameter[i][j] << ", ";
			}*/
			InputProcess = FileInput;
			istringstream ss(InputProcess);
			string token;
			while (std::getline(ss, token, ','))
			{
				stringstream(token) >> KParameter[i][j];
				//cout << KParameter[i][j] << ", ";
				j++;
			}
			//cout << endl;
		}
		//system("pause");
	}
	fr.close();
	clock_t start, stop;
	start = clock(); //�}�l�ɶ�
	float AverageResult = 0;
	float EnergyResume;//�˹�Ѹ�
	int NoAnswer = 0;//�L���`����
	int BestTimes = 0;
	int SecTimes = 0;
	int ThirdTimes = 0;

	int NotBest = 0;//�D�̨θѪ�����
	for (int i = 0; i<RunTime; i++)
	{
		EnergyResume = Empire();
		if (EnergyResume == ImpossiblePowerConsume)
		{
			NoAnswer++;
		}
		AverageResult += EnergyResume;
		if (EnergyResume < the_best+0.001)
		{
			BestTimes++;
		}
	}
	cout << "�̨θѭӼ�" << BestTimes << endl;
	//cout << endl << "�D�̨θѪ����Ƭ�" << NotBest << endl;
	stop = clock(); //�����ɶ�
	cout << "����" << RunTime << "��" << endl;
	cout << "�����ӯ�:" << (AverageResult / RunTime) << "��" << endl;
	cout << "�X�p�O��:" << double(stop - start) / CLOCKS_PER_SEC << "��" << endl << endl;
	cout << "�L�Ѧ���:" << NoAnswer << "��" << endl;
	system("pause");
	return 0;
}

float Empire()
{
	float ENDFlag = ImpossiblePowerConsume;//�ˬd�Ѿl���ɡA�O�_�i�J���Ū��A
	int CirCleOrNot = 0;//�O�_�i�J�`��

	int remain = -1;
	int remain_count = 0;
	int flag = 0;

	Country* AllCountry[NumofCountry];
	int RandomVelocity[NumofLight];
	clock_t start, stop;
	start = clock(); //�}�l�ɶ�

	for (int CountryNo = 0; CountryNo<NumofCountry; CountryNo++)
	{
		int CountryConstruct[NumofLight];

		for (int h = 0; h<NumofLight; h++)
		{
			CountryConstruct[h] = RandomNumber(0, NumofInstruction - 1);                           //cout<<CountryConstruct[h]<<",";
		}

		//��l��m���ͧ���

		for (int MatrixRow = 0; MatrixRow<NumofTd; MatrixRow++)
		{
			int T = CheckT(MatrixRow, CountryConstruct);
			if (T == 0)                                                     //���X��APower�]ImpossiblePowerConsume�A���i��o�ͪ����p�A�@�����
			{                                                                              //cout<<"ImpossiblePower"<<","; 
				for (int h = 0; h<NumofLight; h++)
				{
					RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
				}
				AllCountry[CountryNo] = new Country(CountryConstruct, ImpossiblePowerConsume, RandomVelocity);
				break;
			}
			if (MatrixRow == NumofTd - 1)//����Row���ŦX����A�}�l�p��ӯ�
			{
				float Power = 0;
				for (int k = 0; k<NumofLight; k++)                                                 //�p��ӯ�
				{
					Power += I[CountryConstruct[k]][0];
				}                                                                             //cout<<Power<<",";
				for (int h = 0; h<NumofLight; h++)                                                 //�H����l�t��
				{
					RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
				}
				AllCountry[CountryNo] = new Country(CountryConstruct, Power, RandomVelocity);                //��a����
			}
		}
	}
	//�̷Ӱ�a�����ƧǡA�����V�p�A��O�V�j	
	quickSort(AllCountry, 0, NumofCountry - 1);    //(data,left,right)

	Country* EmpireCountry[NumofEmpire];
	Country* ColonyCountry[NumofEmpire][NumofCountry - NumofEmpire];   //[���ݫҰ��p���s��][�Q�ޥ���s��]�A�@�ӫҰ��p���A�֦����Ǵޥ���
	float EmpireCost[NumofEmpire];
	float EmpirePower[NumofEmpire];
	float TotalCost = 0;
	int EmpireNumofColony[NumofEmpire];  //�@�ӫҰ��p���ү�֦����ޥ���ƥ�
	for (int i = 0; i<NumofEmpire; i++)
	{
		EmpireCountry[i] = AllCountry[i];
		EmpireCost[i] = (EmpireCountry[i]->GetFitness()) - (AllCountry[NumofEmpire]->GetFitness());  //���W�Ʀ���:Cn=����:cn-�̤j����ci:max{ci}
		TotalCost += EmpireCost[i];
	}
	for (int i = 0; i<NumofEmpire; i++)
	{
		EmpirePower[i] = EmpireCost[i] / TotalCost;
		EmpireNumofColony[i] = EmpirePower[i] * (NumofCountry - NumofEmpire);
	}
	int Remain = NumofCountry - NumofEmpire;      //�p��O�_������a�����ޥ��P�Ϊv���Y
	for (int i = 0; i<NumofEmpire; i++)
	{
		Remain = Remain - EmpireNumofColony[i];
	}
	for (int i = 0, j = 0; i<Remain; i++)         //�����ɪ���a�ơA�n�ɶi
	{
		EmpireNumofColony[j]++;
		if (j >= NumofEmpire)
		{
			j = 0;
		}
		j++;
	}
	//for(int i=0;i<NumofEmpire;i++)
	//cout<<EmpireNumofColony[i]<<",";cout<<endl;
	int ColonyNo[NumofEmpire];                  //�[�J�Ұ��p�������ǽs���A�]�Ω�p��ثe���h�ִޥ���,[�ĴX�Ұ��p��]
	for (int j = 0; j<NumofEmpire; j++)
	{
		ColonyNo[j] = 0;
	}
	for (int i = NumofEmpire; i<NumofCountry; i++)
	{
		int EmpireSelection = RandomNumber(0, (NumofEmpire - 1));
		if (ColonyNo[EmpireSelection]<EmpireNumofColony[EmpireSelection])   //�ˬd��ܪ��Ұ��p���O�_�|���B��
		{
			ColonyCountry[EmpireSelection][ColonyNo[EmpireSelection]] = AllCountry[i];
			ColonyNo[EmpireSelection] = ColonyNo[EmpireSelection] + 1;
		}
		else                                                              //�Y�O�B�����ܡA�ݭn���s���
		{
			i--;
		}
	}
	//��l������
	float BalanceCost = 0;
	int Counter = 0;
	int Collapse;
	//�ĤG���}�l
	while (1)         //��Ȧs��̤@�Ұ��p���A�z�L�̫᪺if����break���X
	{
		for (int j = 0; j<NumofEmpire; j++)                 //�̫Ұ��p�������B�z
		{
			for (int k = 0; k<ColonyNo[j]; k++)            //�Ұ��p�������ޥ���̧ǳB�z
			{
				ColonyCountry[j][k]->SetGBest(EmpireCountry[j]->GetPBest());
				ColonyCountry[j][k]->UpdateVelocity();
				ColonyCountry[j][k]->UpdatePosition();
				ColonyCountry[j][k]->CheckPosition();
			}
		}

		for (int j = 0; j<NumofEmpire; j++)                 //�̫Ұ��p�������B�z
		{
			for (int k = 0; k<ColonyNo[j]; k++)            //�Ұ��p�������ޥ���̧ǳB�z
			{
				if ((ColonyCountry[j][k]->GetFitness()) < EmpireCountry[j]->GetFitness())   //�ޥ���«��A�����Ұ�֤�
				{
					Country* SwitchingPointer = ColonyCountry[j][k];
					ColonyCountry[j][k] = EmpireCountry[j];
					EmpireCountry[j] = SwitchingPointer;
				}
			}
		}
		//�j�j�Ұ�m�̮z�p�ꪺ�ޥ���
		int WeakestEmpireNo = -1;
		int AwfulEmpireNo = -1;
		float TotalColonyCost = 0;
		//////////////////////////

		for (int j = 0; j<NumofEmpire; j++)
		{
			//�p��U��j�j�{��
			if (ColonyNo[j] != 0)
			{
				TotalColonyCost = 0;
				for (int k = 0; k<ColonyNo[j]; k++)
				{
					TotalColonyCost += ColonyCountry[j][k]->GetFitness();
				}
				TotalColonyCost = (TotalColonyCost) / ColonyNo[j];
			}
			EmpireCost[j] = ((EmpireCountry[j]->GetFitness() * 7) + (TotalColonyCost * 3)) / 10;//(�Ұ��v�O*7+�����ޥ����v�O*3)/10
		}

		//�Կ��Ԯɪ����ŧ����}��
		int No = RandomNumber(0, NumofEmpire - 1);
		while (ColonyNo[No] == 0)
		{
			No = RandomNumber(0, NumofEmpire - 1);
		}
		EmpireCost[No] -= BalanceCost;
		BalanceCost = 0;
		//
		for (int j = 0; j<NumofEmpire; j++)
		{
			//
			if ((AwfulEmpireNo == -1 || EmpireCost[j] < EmpireCost[AwfulEmpireNo]))   //�M��̱j�j�Ұ�(�����̤p)
			{
				AwfulEmpireNo = j;
			}
			int CheckExist = 0;//�Ѿl�ޥ��a�Ӽ�

			{
				CheckExist += ColonyNo[j];
			}
			if ((WeakestEmpireNo == -1 || EmpireCost[j] > EmpireCost[WeakestEmpireNo]) && (CheckExist != 0))   //�M��{�s�̮z�p�Ұ�(�����̤j)
			{
				WeakestEmpireNo = j;
			}
		}

		///////////////////////////
		////////////////////////////

		//�M��̮z�p�ꪺ�̮z�p�ޥ��a
		Country* WeakestColony = NULL;
		int WeakestColonyNo;
		for (int j = 0; j<ColonyNo[WeakestEmpireNo]; j++)
		{
			if (WeakestColony == NULL || (ColonyCountry[WeakestEmpireNo][j]->GetFitness()) >(WeakestColony->GetFitness()))  //�M��̮z�p�ޥ���(�����̤j)
			{
				WeakestColony = ColonyCountry[WeakestEmpireNo][j];
				WeakestColonyNo = j;
			}
		}

		//�̮z�p�Ұꪺ�̮z�p�ޥ���Q�m
		//�Q�ַm
		int Marauder = -1;
		TotalCost = 0;
		int MarauderPassibility;
		float WeakestCost = EmpireCost[WeakestEmpireNo];
		for (int j = 0; j<NumofEmpire; j++)
		{
			int CheckExist = 0;//�Ѿl�ޥ��a�Ӽ�

			{
				CheckExist += ColonyNo[j];
			}
			if (CheckExist != 0)
			{
				EmpireCost[j] = EmpireCost[j] - WeakestCost;      //���W�Ʀ���
				TotalCost += EmpireCost[j];
			}
		}
		while (Marauder == WeakestEmpireNo || Marauder == -1)
		{
			MarauderPassibility = RandomNumber(1, 100);
			for (int j = 0; j<NumofEmpire; j++)
			{
				int CheckExist = 0;//�Ѿl�ޥ��a�Ӽ�

				CheckExist += ColonyNo[j];

				if (CheckExist != 0)
				{
					EmpirePower[j] = EmpireCost[j] / TotalCost;                      //�v�O�p��
					//�O���O�o�ӫҰ�m��?
					if (MarauderPassibility <= (EmpirePower[j] * 100))
					{
						Marauder = j;
						break;
					}
					else
					{
						MarauderPassibility -= EmpirePower[j] * 100;
					}
					//
				}
			}
			//cout<<WeakestEmpireNo<<","<<ColonyNo[WeakestEmpireNo]<<";"<<Marauder<<","<<ColonyNo[Marauder]<<endl;
		}
		//�Q�m�F

		ColonyCountry[Marauder][ColonyNo[Marauder]] = ColonyCountry[WeakestEmpireNo][WeakestColonyNo];
		//�z�p��ޥ��a��z
		if (WeakestColonyNo != (ColonyNo[WeakestEmpireNo] - 1))//�̮z�p���ޥ��a�A���O�s���̫᭱���ޥ��a��
		{
			ColonyCountry[WeakestEmpireNo][WeakestColonyNo] = ColonyCountry[WeakestEmpireNo][ColonyNo[WeakestEmpireNo] - 1];
			ColonyCountry[WeakestEmpireNo][ColonyNo[WeakestEmpireNo] - 1] = NULL;
		}
		else
		{
			ColonyCountry[WeakestEmpireNo][WeakestColonyNo] = NULL;
		}
		ColonyNo[WeakestEmpireNo]--;//�ޥ��a�֤@��
		//
		ColonyNo[Marauder]++;//�ޥ��a�h�@��


		/////////////////////////////
		Collapse = 0;
		for (int j = 0; j<NumofEmpire; j++)
		{

			if (ColonyNo[j] == 0)
			{
				Collapse++;
			}
		}
		if (Collapse >= (NumofEmpire - 1))  //�Ѿl�@��
		{
			break;
		}

		else if (Collapse >= (NumofEmpire - 2))  //�Ѿl2��A��O������A�i�J���Ū��A�A�j�z���Y���g����2��
		{
			if (CirCleOrNot == 0)
			{
				//cout<<"�j�p��:"<<EmpireCost[AwfulEmpireNo]<<endl;
				//cout<<"�z�p��:"<<EmpireCost[WeakestEmpireNo]<<endl;
				ENDFlag = EmpireCost[AwfulEmpireNo];//�j��z�ꤣ�_��աA�ҥH�S�w��a�����j�ꪺ�g���]�O���
				CirCleOrNot = 1;
			}
			else if (CirCleOrNot == 1)
			{
				CirCleOrNot = 2;
			}
			else if (CirCleOrNot == 2 && ENDFlag == EmpireCost[AwfulEmpireNo])
			{
				break;
			}
			else
			{
				CirCleOrNot = 0;
			}
		}
		//�h��Կ�
		float Balance = 0;
		for (int i = 0; i<NumofEmpire - 1; i++)
		{
			if (ColonyNo[i] != 0)
			{
				Balance += EmpireCost[i];
			}
		}
		if (-(Balance / (NumofEmpire - Collapse)) <= -(EmpireCost[AwfulEmpireNo] * (1 - BalanceRatio / 100)))
		{
			BalanceCost = EmpireCost[AwfulEmpireNo] * BalanceRatio / 100;
		}

		//�[�ݾԧ�
		Counter++;
		if (Counter >= 100)
		{
			flag = 1;
		}
		if (Counter == 1000)
		{
			cout << "���椤,�Ѿl" << NumofEmpire - Collapse << "�Ұ�,�̱j:" << EmpireCost[AwfulEmpireNo] << ",�̮z:" << EmpireCost[WeakestEmpireNo] << endl;
			Counter = 0;
		}

		if (flag == 1)
		{
			if (remain == -1)
			{
				remain = NumofEmpire - Collapse;
			}
			else if (remain == NumofEmpire - Collapse)
			{
				remain_count++;
			}
			else
			{
				remain = NumofEmpire - Collapse;
				remain_count = 0;
			}
			if (remain_count == 5)
			{
				break;
			}
		}
	}

	int BetterEmpireNo = -1;
	for (int i = 0; i<NumofEmpire; i++)
	{
		float BetterEmpireFitness = ImpossiblePowerConsume;
		if (ColonyNo[i] != 0)
		{
			if (EmpireCountry[i]->GetFitness() < BetterEmpireFitness || BetterEmpireNo == -1)
			{
				BetterEmpireNo = i;
			}
		}
	}
	//cout<<EmpireCountry[BetterEmpireNo]->GetFitness()<<endl;
	/********************/
	cout << "���O:";
	for (int j = 0; j<NumofLight; j++)
	{
		cout << *(EmpireCountry[BetterEmpireNo]->GetPBest() + j) << ",";
	}
	cout << "�ӯ�:" << EmpireCountry[BetterEmpireNo]->GetFitness()<<", ";

	stop = clock(); //�����ɶ�
	cout << "�O��:" << double(stop - start) / CLOCKS_PER_SEC << "��" << endl;
	//system("pause");
	/********************/

	return EmpireCountry[BetterEmpireNo]->GetFitness();
}

//�T�{�G�צX��P�_
int CheckT(int MatrixRow, int i[])     //MatrixRow���n�ˬd�O�_�X�檺row�Ai���������O���x�}
{
	float T = 0;
	int j;
	for (int k = 0; k<(NumofLight + OtherLight); k++)
	{
		if (k >= NumofLight)
		{
			//	T+=KParameter[MatrixRow][k]*OtherLightT[k-NumofLight];
		}
		else
			T += KParameter[MatrixRow][k] * I[i[k]][1];//[��ܿO�����O,1��ܫG��]
	}
	if (T<Td[MatrixRow])
		return 0;
	else
		return 1;
}

void quickSort(Country* EmpireCountry[], int left, int right)
{
	int i = left, j = right;
	Country* SortingTmp;
	float Pivot = EmpireCountry[(left + right) / 2]->GetFitness();

	// partition 
	while (i <= j) {
		while ((EmpireCountry[i]->GetFitness() < Pivot))
			i++;
		while (EmpireCountry[j]->GetFitness() > Pivot)
			j--;
		if (i <= j) {
			SortingTmp = EmpireCountry[i];
			EmpireCountry[i] = EmpireCountry[j];
			EmpireCountry[j] = SortingTmp;
			i++;
			j--;
		}
	};
	// recursion 
	if (left < j)
		quickSort(EmpireCountry, left, j);
	if (i < right)
		quickSort(EmpireCountry, i, right);
}

/****************************************Random Number**********************************************
�����Ʊ檺�̤j�ȻP�̤p�ȡA�N�ಣ�ͤ���o��ӭȤ����H���ơA�o���H���üƤ]�i�൥��̤j�ȡA�ε���̤p�ȡC
��üƭȬO�̷Ӯɶ��Ҳ��͡A�ҥH���ƩI�s�ɡA�|�������P���üƭ�(���)�C
***************************************************************************************************/
int RandomNumber(int MinValue, int MaxValue)   //Both MinValue and MaxValue are included
{
	int R = (rand() % (MaxValue - MinValue + 1)) + MinValue;
	return R;
}