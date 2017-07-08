#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;



//Ū�ɨ��o
int NumofLight;       //�O���Ӽ�
int NumofInstruction; //�G�׫��O�Ӽ�
int NumofTd;          //�ؼ��I�Ӽ�

int RunTime = 1;  //�]�h�֦�
float RealBest = -1;

#define OtherLight 0        //�B�~����-�Ӷ�
//float OtherLightT[OtherLight]={0};        //�B�~�������G��

int NumofParticle;        //��NumofParticle�ոѧ@��PSO�����(�ɤl��)
int NumofRedistribute;

float RepeatValue = -1;
float DifferentRange;
int NumofGroup;

//�t��k�w�q
#define ImpossibleResume 99999
#define InitialW 0.9             //�D�ʫY�ơA�Ѥj��p
#define FinalW 0.4
#define InitialCOne 2.5          //������C�Y�ơA�Ѥj��p�A�]���j�M�����������
#define FinalCOne 1
#define InitialCTwo  1       //�s�骺C�Y�ơA�Ѥp��j�A�]���j�M���������s��
#define FinalCTwo 2.5
int CircleofPSO;
int VelocityLimit;
int NumofGroupCircle;
double ResumeTime;

void Parameter()
{
	CircleofPSO = 50;//NumofParticle*(0.02);      //PSO���檺�`�����ơA�ɤl�ƪ�2%

	VelocityLimit = NumofInstruction/3;//NumofInstruction*(0.5);  //���ʳt�פW��

	NumofGroupCircle = 10;//NumofGroup;    //�U�B�z������y�һݡA�����B�z���Ӽ�

	NumofRedistribute = NumofParticle / 2;  //�n�洫���ɤl��
}
//


int CheckT(float T, int MatrixRow, int i, int j, int k);

void PSO();

void* new2d(int h, int w, int size);
int RandomNumber(int MinValue, int MaxValue);


float** I;
float* Td;
float** KParameter;

//�B�~���i��������B�z(����)
/*
void OtherLightProcess()
{
for(int i=0;i<NumofTd;i++)
{
Td[i]-=KParameter[i][NumofLight+OtherLight-1]*OtherLightT[OtherLight];
}
}
*/
class Particle  //�C�Ӳɤl���򥻳]�w�A�ܰʽե����O����
{
private:
	float* PBest=new float[NumofLight + 1],  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,����̧C�ӯ�)
		*GBest = new float[NumofLight + 1];  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,�s��̧C�ӯ�)
	int* Position = new int[NumofLight];  //�O������A��U���O�����ե����O <----> Position   //�䤤�O���Ƥ]�Osolution�����סA�ĴX���O��ܲĴX��
	int* Velocity = new int[NumofLight];             //�ե����O�W�ɩΤU�� <----> Velocity
	float W;//�����D�ʥΪ�weighting�A��ȶV�j�A�V�A�X����j�M�A��ĳ��0.9~0.4
	float C1, C2;//�Y�ơA�վ㰾�V�U��A�ΰ��V�s��
public:
	Particle(float* PPointer, float Power, int* VPointer);
	Particle(int RangeofLightwithInstruction[], float Power, int* VPointer);
	void SetW(float NewW, float NewC1, float NewC2)
	{
		W = NewW;
		C1 = NewC1;
		C2 = NewC2;
	}
	float* GetPBest(){ return &PBest[0]; }
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
		{
			Position[i] = Position[i] + Velocity[i];
			if (Position[i]>NumofInstruction - 1)             //�ե����O���W�U��
				Position[i] = NumofInstruction - 1;
			if (Position[i]<0)
				Position[i] = 0;
		}
	}
	void UpdateVelocity()
	{
		for (int i = 0; i<NumofLight; i++)//i�O����
		{
			Velocity[i] = (W*Velocity[i] + C1*((float)RandomNumber(0, 100) / 100)*(PBest[i] - Position[i]) + C2*((float)RandomNumber(0, 100) / 100)*(GBest[i] - Position[i]));
			if (Velocity[i]>VelocityLimit)                   //�t�פW�U��
				Velocity[i] = VelocityLimit;
			if (Velocity[i]<(-VelocityLimit))
				Velocity[i] = (-VelocityLimit);
		}
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
		{
			Fitness += I[Position[i]][0];
		}
		if (Fitness<PBest[NumofLight] || PBest[NumofLight] <= 0)
		{
			for (int j = 0; j<NumofLight; j++)
			{
				PBest[j] = Position[j];
			}
			PBest[NumofLight] = Fitness;
		}
	}

};

Particle::Particle(float* PPointer, float Power, int* VPointer)
{
	for (int i = 0; i<NumofLight; i++)
	{
		Position[i] = *(PPointer + i);
		Velocity[i] = *(VPointer + i);
		PBest[i] = Position[i];
	}
	PBest[NumofLight] = Power;
	W = InitialW;
	C1 = InitialCOne;
	C2 = InitialCTwo;
}
Particle::Particle(int RangeofLightwithInstruction[], float Power, int* VPointer)
{
	for (int i = 0; i<NumofLight; i++)
	{
		Velocity[i] = *(VPointer + i);
	}
	for (int k = 0; k<NumofLight; k++)
	{
		Position[k] = RangeofLightwithInstruction[k];
		PBest[k] = RangeofLightwithInstruction[k];
	}
	PBest[NumofLight] = Power;
	W = InitialW;
	C1 = InitialCOne;
	C2 = InitialCTwo;
}
//PSO���O����



int main()
{
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));

	//Number of Particles
	fr.open("PSOParameter.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open NumofParticle file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%d", &NumofParticle);
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%f", &DifferentRange);
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%d", &NumofGroup);
	}
	fr.close();

	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open LightInstruction file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofInstruction
		sscanf_s(FileInput, "%d", &NumofInstruction);

		I = (float **)new2d(NumofInstruction, 2, sizeof(float));

		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf_s(FileInput, "%f", &I[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open Td file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofTd
		sscanf_s(FileInput, "%d", &NumofTd);

		Td = new float[NumofTd];

		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf_s(FileInput, "%f", &Td[i]);        //�r����Ʀr
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open KParameter file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}��
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}�e
		sscanf_s(FileInput, "%d", &NumofLight);

		KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));

		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf_s(FileInput, "%f", &KParameter[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();
	//���Ū������

	Parameter();

	//OtherLightProcess();//�B�~�����B�z

	/*******PSO���������G**********/
	PSO();

	/***����Ŷ�***/
	delete[] I;
	delete[] Td;
	delete[] KParameter;

	//system("pause");
	return 0;
}

//�T�{�G�צX��P�_
int CheckT(int MatrixRow, int i[])     //MatrixRow���n�ˬd�O�_�X�檺row�Ai���������O���x�}
{
	float T = 0;
	for (int k = 0; k<NumofLight; k++)
	{
		T += KParameter[MatrixRow][k] * I[i[k]][1];//[��ܿO�����O,1��ܫG��]
	}
	if (T<Td[MatrixRow])
		return 0;
	else
		return 1;
}
/**********************************************************************************************/
void PSO()
{
	clock_t PSOstart, PSOstop;
	float* TotalGBest=new float[NumofLight + 1];
	TotalGBest[NumofLight] = ImpossibleResume;

	float** GBest = (float**)new2d(NumofGroup, NumofLight + 1, sizeof(float));
	float *PBestPoint, Power = 0;
	for (int i = 0; i<NumofGroup; i++)
		GBest[i][NumofLight] = ImpossibleResume;
	int ParticleNo = 0;                  //�ؼЬO���NumofParticle�Ӹѧ@����l�ȡA�ҥH�]�@�ӺX�ХΥH�p��
	Particle*** Group=(Particle***)new2d(NumofGroup,NumofParticle,sizeof(Particle*));
	float T = 0;

	/*********���������t�ɨϥ�************/
	Particle*** TempGroup=(Particle***)new2d(NumofGroup,NumofParticle,sizeof(Particle*));
	int* TempGroupNo=new int[NumofGroup];  //�p��洫�ɤl�Ϊ��Ȧs�s�O�_�w��
	for (int i = 0; i<NumofGroup; i++)
	{
		TempGroupNo[i] = 0;
	}
	/**************************************/
	PSOstart = clock();
	/********************************�H�����ͪ�l�ȡA���ר�O�_�X��*****************************************/
	int* RandomVelocity=new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)                                       //����j�������зj�M"��������"�Aj�s�F�Ҧ��O�����O��T
		{
			int* j=new int[NumofLight];
			//�ɤl��l��m�H�����ͤ覡3
			for (int h = 0; h<NumofLight; h++)
			{
				j[h] = RandomNumber(0, NumofInstruction - 1);
			}

			//
			//��l��m���ͧ���

			for (int MatrixRow = 0; MatrixRow<NumofTd; MatrixRow++)
			{
				T = CheckT(MatrixRow, j);
				if (T == 0)                                                     //���X��APower�]ImpossibleResume�A���i��o�ͪ����p�A�@�����
				{
					for (int h = 0; h<NumofLight; h++)
					{
						RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
					}
					Group[GroupNo][ParticleNo++] = new Particle(j, ImpossibleResume, &RandomVelocity[0]);
					break;
				}
				else
				{
					T = 0;
				}
				if (MatrixRow == NumofTd - 1)//����Row���ŦX����A�}�l�p��ӯ�
				{
					for (int k = 0; k<NumofLight; k++)                                                 //�p��ӯ�
					{
						Power += I[j[k]][0];
					}
					for (int h = 0; h<NumofLight; h++)                                                 //�H����l�t��
					{
						RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
					}
					Group[GroupNo][ParticleNo++] = new Particle(j, Power, &RandomVelocity[0]);                //���Ͳɤl
					Power = 0;
				}
			}
		}
		for (int i = 0; i<NumofParticle; i++)
		{
			PBestPoint = Group[GroupNo][i]->GetPBest();
			if (*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
			{
				for (int j = 0; j<NumofLight + 1; j++)
				{
					GBest[GroupNo][j] = *(PBestPoint + j);
				}
			}
		}
	}

	//��Ȳ��ͧ���
	//�P�ɡA�Ĥ@��PSO�����A�ѤUCircleofPSO-1��
	//�s�W����B�z�A�b�̥~�h
	/*******************************************************************************************************/
	for (int GroupCircle = 0; GroupCircle<NumofGroupCircle; GroupCircle++)
	{
		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			for (int i = 1; i<CircleofPSO; i++)
			{
				float NewW = InitialW - (InitialW - FinalW)*((i + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
				float NewC1 = InitialCOne - (InitialCOne - FinalCOne)*((i + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
				float NewC2 = InitialCTwo + (FinalCTwo - InitialCTwo)*((i + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
				for (int i = 0; i<NumofParticle; i++)
				{
					Group[GroupNo][i]->SetW(NewW, NewC1, NewC2);
					Group[GroupNo][i]->SetGBest(&GBest[GroupNo][0]);
					Group[GroupNo][i]->UpdateVelocity();
					Group[GroupNo][i]->UpdatePosition();
					Group[GroupNo][i]->CheckPosition();
				}
				for (int i = 0; i<NumofParticle; i++)
				{
					PBestPoint = Group[GroupNo][i]->GetPBest();
					if ((*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume) && *(PBestPoint + NumofLight) != ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
					{
						for (int j = 0; j<NumofLight + 1; j++)
						{
							GBest[GroupNo][j] = *(PBestPoint + j);
						}
					}
				}
				if (GBest[GroupNo][NumofLight]<TotalGBest[NumofLight] || TotalGBest[NumofLight] == ImpossibleResume)
				{
					for (int p = 0; p<NumofLight + 1; p++)
						TotalGBest[p] = GBest[GroupNo][p];
				}
			}
		}


		//�����覡
		/****************************************************
		//��ҥ��]����
		int RedistributeList[NumofGroup];
		int RedistributeParticleList[NumofGroup][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;

		//�إ߸s��list
		for(int i=0;i<NumofGroup;i++)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		for(int j=0;j<i;j++)
		{
		if(RedistributeList[j]==RandomGroup)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		j=-1;                                                             //for loop���Y�A�ˬd�e�|�[1�A�ҥH�]-1�A���L�q0�}�l�ˬd
		}
		}
		RedistributeList[i]=RandomGroup;
		}

		//��ܭn�洫���ɤl�A�إ߹�Ӫ�
		for(int k=0;k<NumofGroup;k++)
		{
		for(int i=0;i<NumofRedistribute;i++)
		{
		RandomParticle=RandomNumber(0,NumofParticle-1);
		for(int j=0;j<i;j++)
		{
		if(RedistributeParticleList[k][j]==RandomParticle)
		{
		RandomParticle=RandomNumber(0,NumofParticle-1);
		j=-1;                                                //for loop���Y�A�ˬd�e�|�[1�A�ҥH�]-1�A���L�q0�}�l�ˬd
		}
		}
		RedistributeParticleList[k][i]=RandomParticle;
		}
		}

		//�}�l�洫
		for(int i=0;i<NumofGroup;i++)  //�ƥ�
		{
		for(int j=0;j<NumofParticle;j++)
		{
		TempGroup[i][j]=Group[i][j];
		}
		}
		for(int i=0;i<NumofGroup;i++) //i�O�Ҧb���s,k�O�ĴX���@�s�������洫
		{
		for(int j=0;j<NumofRedistribute;j++)
		{
		Group[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempGroup[i][ RedistributeParticleList[i][j] ];
		}
		}
		/************��������************
		int RandomGroup;
		for(int i=0;i<NumofGroup;i++)
		{
		for(int j=0;j<NumofParticle;j++)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		while(TempGroupNo[RandomGroup]==NumofParticle)//�Ȧs�s�w��
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		}
		TempGroup[RandomGroup][TempGroupNo[RandomGroup]++] = Group[i][j];
		}
		}
		for(int i=0;i<NumofGroup;i++)
		{
		for(int j=0;j<NumofParticle;j++)
		{
		Group[i][j]=TempGroup[i][j];
		}
		TempGroupNo[i]=0;
		}
		/*******************************/
		/********************************/
		//////////////�����������s�A���s�զX
		int RandomGroup;
		int RandomParticle1, RandomParticle2;
		Particle *TempMemory;
		for (int j = 0; j<NumofGroup; j++)
		{
			for (int i = 0; i<NumofRedistribute; i++)
			{
				RandomGroup = RandomNumber(0, NumofGroup - 1);
				while (RandomGroup == j)
				{
					RandomGroup = RandomNumber(0, NumofGroup - 1);
				}
				RandomParticle1 = RandomNumber(0, NumofParticle - 1);
				RandomParticle2 = RandomNumber(0, NumofParticle - 1);
				//
				TempMemory = Group[j][RandomParticle1];
				Group[j][RandomParticle1] = Group[RandomGroup][RandomParticle2];
				Group[RandomGroup][RandomParticle2] = TempMemory;
			}
		}
		///////////*///��������

		/*�M�䭫���᪺�U�s��gbest*/
		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
		    for(int i=0;i<NumofParticle;i++)
		    {
				GBest[GroupNo][NumofLight] = ImpossibleResume;
		    }
		    for(int i=0;i<NumofParticle;i++)
		    {
				PBestPoint = Group[GroupNo][i]->GetPBest();
				if ((*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume) && *(PBestPoint + NumofLight) != ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
		        {
		            for(int j=0;j<NumofLight+1;j++)
		            {
						GBest[GroupNo][j] = *(PBestPoint + j);
		            }
		        }
		    }
		}

		/****�M�䧹��****/

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			/*�U�s���s��X��gbest�g�^*/
			for(int i=0;i<NumofParticle;i++)
			{
				Group[GroupNo][i]->SetGBest(&GBest[GroupNo][0]);
			}
			/******/
			/*�N��X��Gbest�g�^���U�s*
			for (int i = 0; i<NumofLight + 1; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			/*********/
		}

		//�{���D�X���Ѥwí�w
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
			RepeatValue = TotalGBest[NumofLight];
		else if (RepeatValue - TotalGBest[NumofLight]<TotalGBest[NumofLight] * DifferentRange && RepeatValue != -1)
			break;
	}

	PSOstop = clock(); //�����ɶ�

	/***PSO���浲���A�[�ݵ��G***
	if (TotalGBest[NumofLight] != ImpossibleResume)
	{
		//cout<<"PSO�j�M"<<endl
		cout << "���O:";
		for (int i = 0; i<NumofLight; i++)
		{
			cout << TotalGBest[i] << ",";
		}
		cout << "�ӯ�:" << TotalGBest[NumofLight] << endl;
	}
	else
	{
		cout << "�L��" << endl;
	}

	/****�g��****/
	if (TotalGBest[NumofLight] <= RealBest || RealBest == -1)
	{
		fstream fw;
		string filename = "The_Best_Solution.txt";
		fw.open(filename, ios::out);
		if (TotalGBest[NumofLight] != ImpossibleResume)
		{
			for (int i = 0; i < NumofLight; i++)
			{
				fw << TotalGBest[i] << endl;
			}
			fw << TotalGBest[NumofLight] << endl;
		}
		else
		{
			fw << "�L��" << endl;
		}
		fw << double(PSOstop - PSOstart) / CLOCKS_PER_SEC << endl;
		fw.close();
	}

	/*�U�I�ӫ�*/
	fstream fw;
	string filename = "RealIllumination.txt";
	fw.open(filename, ios::out);
	float RealT = 0;
	for (int MatrixRow = 0; MatrixRow < NumofTd; MatrixRow++)
	{
		for (int k = 0; k < NumofLight; k++)
		{
			RealT += KParameter[MatrixRow][k] * I[(int)TotalGBest[k]][1];//[��ܿO�����O,1��ܫG��]
		}
		fw << RealT << endl;
		RealT = 0;
	}
	fw.close();
	/****�g�ɧ���****/

	cout << "�O��:" << double(PSOstop - PSOstart) / CLOCKS_PER_SEC << "��" << endl << endl;

	/*****************************/
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)  //����Ŷ�
	{
		for (int i = 0; i<NumofParticle; i++)
		{
			delete Group[GroupNo][i];
		}
	} /**/

	//return TotalGBest[NumofLight];
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

/**********************************************
�ʺA����2��array
/**********************************************/
void* new2d(int h, int w, int size)
{
	register int i;
	void **p;

	p = (void**)new char[h*sizeof(void*) + h*w*size];
	for (i = 0; i < h; i++)
	{
		p[i] = ((char *)(p + h)) + i*w*size;
	}

	return p;
}