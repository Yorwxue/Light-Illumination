
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <curand_kernel.h>

#include <stdio.h>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//�ϥΪ̩w�q
int NumofLight;        //�O���Ӽ�
int NumofInstruction; //�G�׫��O�Ӽ�
int NumofTd;          //�ؼ��I�Ӽ�

#define MaxNumofLight 12        //�O���Ӽ�
#define MaxNumofInstruction 16 //�G�׫��O�Ӽ�
#define MaxNumofTd 5          //�ؼ��I�Ӽ�

int NumofGroup;  //���s�ƥءA�W�[���ѼƯണ�ɮį�

float RepeatValue = -1;
float DifferentRange;
int NumofParticle;       //��NumofParticle�ոѧ@��PSO�����(�ɤl��) �A�C��particle�j�p��200
int NumofRedistribute;



#define OtherLight 0        //�B�~�������Ӷ�

#define ThreadX 8  ///�\��y��X��thread
#define ThreadY 8  ///�\��y��Y��thread

#define BlockX  128 ///�\��y��X��block 
#define BlockY  128 ///�\��y��Y��block

dim3 ThreadPerBlocks(ThreadX, ThreadY);//threads�bblock�����\��覡�A�y��X�PY
dim3 NumofBlocks(BlockX / ThreadPerBlocks.x, BlockY / ThreadPerBlocks.y);//block���\��覡�A�y��X�PY
//float OtherLightT[OtherLight]={0};        //�B�~�������ӫ�


//���ե�
int RunTime = 1;  //����h�֦�
float RealBest = -1;
#define in_cir_stable 2
#define in_cir_threshold 0.5


//�{���۩w�q
#define ImpossibleResume 99999
#define InitialW 0.9             //�D�ʫY�ơA�Ѥj��p
#define FinalW 0.4
#define InitialCOne 2.5          //������C�Y�ơA�Ѥj��p�A�]���j�M�����������
#define FinalCOne 1
#define InitialCTwo  1       //�s�骺C�Y�ơA�Ѥp��j�A�]���j�M���������s��
#define FinalCTwo 2.5

#define CircleofPSO 50
#define VelocityLimit 5
#define NumofGroupCircle 10



float Host_I[MaxNumofInstruction][2];   //(�ӯ�,�G��)
float* Host_Td;
float Host_KParameter[MaxNumofTd][MaxNumofLight + OtherLight];

int CheckT(float T, int MatrixRow, int i, int j, int k);
void PSO();
int RandomNumber(int MinValue, int MaxValue);
void* new2d(int h, int w, int size);
float similarity(int* A, int* B, int Length);
__device__ float device_similarity(int* A, int* B, int Length);
__device__ int device_RandomNumber(unsigned int thread_id, int MinValue, int MaxValue);

class Particle  //�C�Ӳɤl���򥻳]�w�A�ܰʽե����O����
{
private:
	int PBest[MaxNumofLight];  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,����̧C�ӯ�)
	float PBestFitness;
	int GBest[MaxNumofLight];
	float GBestFitness;
	int Position[MaxNumofLight];  //�O������A��U���O�����ե����O <----> Position   //�䤤�O���Ƥ]�Osolution�����סA�ĴX���O��ܲĴX��
	int Velocity[MaxNumofLight];             //�ե����O�W�ɩΤU�� <----> Velocity

public:
	Particle(int* PPointer, int* VPointer);
	Particle(){ PBestFitness = ImpossibleResume; }
	__host__  int* GetPBest(){ return PBest; }
	__host__  float GetPBestFitness(){ return PBestFitness; }
	__host__  void CheckPosition()
	{
		float Trow = 0;  //�ھڥثe�ɤl����m�Һ�X�Ӫ�T(�ӫ�)��
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight); j++)
			{
				Trow += Host_KParameter[i][j] * Host_I[Position[j]][1];
			}
			if (Trow<Host_Td[i])
			{
				break;
			}
			Trow = 0;
			if (i == (NumofTd - 1))
			{
				ComputingFitness();
			}
		}
	}
	__host__  void ComputingFitness()
	{
		float Fitness = 0;
		for (int i = 0; i<NumofLight; i++)
		{
			Fitness += Host_I[Position[i]][0];
		}

		if (Fitness<PBestFitness || PBestFitness == ImpossibleResume)
		{
			for (int j = 0; j<NumofLight; j++)
			{
				PBest[j] = Position[j];
			}
			PBestFitness = Fitness;
		}
	}

	__device__  int* d_GetPBest(){ return PBest; }
	__device__  float d_GetPBestFitness(){ return PBestFitness; }
	__device__  void DeviceFunction(int* GBestPointer, float NewGBestFitness, float C0, float C1, float C2, float* I, float* Td, float* KParameter, int NumofLight, int NumofInstruction, int NumofTd, int ParticleNoinGroup)
	{
		//SetGBest
		GBestFitness = NewGBestFitness;

		//SetGBest & UpdateVelocity & Position 
		for (int i = 0; i<NumofLight; i++)//i�O����
		{
			//SetGBest
			GBest[i] = *(GBestPointer + i);

			//UpdateVelocity
			Velocity[i] = (C0*Velocity[i] + C1*(float)device_RandomNumber(ParticleNoinGroup, 0, 100) / 100 * (PBest[i] - Position[i]) + C2*(float)device_RandomNumber(ParticleNoinGroup,0, 100) / 100 * (GBest[i] - Position[i]));
			if (Velocity[i]<(-VelocityLimit))                //�t�פW�U��
			{
				Velocity[i] = (-VelocityLimit);
			}
			else if (Velocity[i]>VelocityLimit)                   //�t�פW�U��
			{
				Velocity[i] = VelocityLimit;
			}

			//UpdatePosition 
			Position[i] = Position[i] + Velocity[i];
			if (Position[i]<0)                                //�ե����O���W�U��
			{
				Position[i] = 0;
			}
			else if (Position[i]>(NumofInstruction - 1))             //�ե����O���W�U��
			{
				Position[i] = (NumofInstruction - 1);
			}
		}
		//CheckPosition
		float Trow = 0;  //�ھڥثe�ɤl����m�Һ�X�Ӫ�T(�ӫ�)��     

		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				Trow += *(KParameter + i*(MaxNumofLight + OtherLight) + j) * *(I + 2 * Position[j] + 1);//KParameter*I;  
			}
			if (Trow<*(Td + i))
			{
				break;
			}
			Trow = 0;
			if (i == (NumofTd - 1))
			{
				//ComputingFitness
				float Fitness = 0;
				for (int i = 0; i<NumofLight; i++)
				{
					Fitness += *(I + 2 * Position[i] + 0);
				}

				if (Fitness<PBestFitness)
				{
					for (int j = 0; j<NumofLight; j++)
					{
						PBest[j] = Position[j];
					}
					PBestFitness = Fitness;
				}
				//ComputingFitnessEND
			}
		}
		//CheckPositionEND
	}

};

Particle::Particle(int* PPointer, int* VPointer)
{
	for (int i = 0; i<NumofLight; i++)
	{
		Position[i] = *(PPointer + i);
		Velocity[i] = *(VPointer + i);
		PBest[i] = Position[i];
	}
	PBestFitness = ImpossibleResume;
}//PSO���O����


__global__ void PSOKernel(Particle* Group, int* GBest, float* GBestFitness, int GroupCircle, float* I, float* Td, float* KParameter, int NumofLight, int NumofInstruction, int NumofTd, int NumofParticle)
{
	int CoordinatesX = blockIdx.x * blockDim.x + threadIdx.x;
	int CoordinatesY = blockIdx.y * blockDim.y + threadIdx.y;
	int ParticleNoinGroup = CoordinatesY*(ThreadX*BlockX) + CoordinatesX;//Y*(row size)+(X)
	if (ParticleNoinGroup < NumofParticle && CoordinatesX < (ThreadX*BlockX) && CoordinatesY < (ThreadY*BlockY))
	{
		float last_sim = -1;
		int stable_count = 0;
		for (int circle = 1; circle < CircleofPSO; circle++)
		{
			float C0 = InitialW - (InitialW - FinalW)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
			float C1 = InitialCOne - (InitialCOne - FinalCOne)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
			float C2 = InitialCTwo + (FinalCTwo - InitialCTwo)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));

			(Group + ParticleNoinGroup)->DeviceFunction(GBest, *(GBestFitness), C0, C1, C2, I, Td, KParameter, NumofLight, NumofInstruction, NumofTd, ParticleNoinGroup);

			if (((Group + ParticleNoinGroup)->d_GetPBestFitness()) < *(GBestFitness))
			{
				for (int j = 0; j < NumofLight; j++)
				{
					*(GBest + j) = *(((Group + ParticleNoinGroup)->d_GetPBest()) + j);
				}
				*(GBestFitness) = (Group + ParticleNoinGroup)->d_GetPBestFitness();
			}
			
			//similarity
			float sim = 0;
			sim = device_similarity((Group + ParticleNoinGroup)->d_GetPBest(), (Group + ParticleNoinGroup)->d_GetPBest(), NumofLight);

			float temp_check = (sim - last_sim) / last_sim;
			if ((temp_check >= 0 && temp_check <= in_cir_threshold) || (temp_check < 0 && temp_check >= -in_cir_threshold))
			{
				stable_count++;
			}
			if (stable_count >= in_cir_stable)
			{
				break;
			}
			//circle--;
			last_sim = sim;
		}
	}
}


int main()
{
	/***********************Ū���ɮ�**********************/
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));

	//Number of Particles
	fr.open("PSOParameter.txt", ios::in); //C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\PSOParameter.txt
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open PSOParameter file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%d", &NumofParticle);
		NumofRedistribute = NumofParticle / 2;
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%f", &DifferentRange);
		fr.getline(FileInput, sizeof(FileInput), ',');
		sscanf_s(FileInput, "%d", &NumofGroup);
	}
	fr.close();

	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\LightInstruction.txt
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofInstruction
		sscanf(FileInput, "%d", &NumofInstruction);
		//Host_I = (float**)new2d(MaxNumofInstruction, 2, sizeof(float));
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_I[i][j]);        //�r����Ʀr
				//cout << Host_I[i][j] << ",";
			}
			//cout<<endl;
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\Td.txt
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofTd
		sscanf(FileInput, "%d", &NumofTd);
		Host_Td = new float[NumofTd];
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf(FileInput, "%f", &Host_Td[i]);        //�r����Ʀr
			//cout<<Host_Td[i]<<endl;
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\KParameter.txt
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}��
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}�e
		sscanf(FileInput, "%d", &NumofLight);
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_KParameter[i][j]);        //�r����Ʀr
				//cout<<Host_KParameter[i][j]<<",";
			}
			//cout<<endl;
		}
	}
	fr.close();
	//���Ū������

	for (int i = 0; i<RunTime; i++)
	{
		PSO();
	}
	//system("pause");
	return 0;
}




void PSO()
{
	fstream fw;
	string filename;

	clock_t StartTime, EndTime;
	int* TotalGBest = new int[NumofLight];
	float TotalGBestFitness = ImpossibleResume;

	int** GBest = (int**)new2d(NumofGroup, MaxNumofLight, sizeof(int));
	float* GBestFitness = new float[NumofGroup];

	for (int i = 0; i<NumofGroup; i++)
	{
		GBestFitness[i] = ImpossibleResume;
	}

	int ProduceParticleNo = 0;

	Particle* GroupPointer;
	Particle** GroupInfo = (Particle**)new2d(NumofGroup, NumofParticle, sizeof(Particle));

	float total_ave_sim = 0;

	////////for CUDA//////////////////////

	//GPU�ǳƪŶ�
	float* I;   //(�ӯ�,�G��)
	float* Td;
	float* KParameter;

	size_t SizeofI = MaxNumofInstruction * 2 * sizeof(float);
	cudaMalloc(&I, SizeofI);

	size_t SizeofTd = NumofTd*sizeof(float);
	cudaMalloc(&Td, SizeofTd);

	size_t SizeofKParameter = MaxNumofTd*(MaxNumofLight + OtherLight)*sizeof(float);
	cudaMalloc(&KParameter, SizeofKParameter);

	//�N��Ʋ���GPU
	cudaMemcpy(I, Host_I, SizeofI, cudaMemcpyHostToDevice);
	cudaMemcpy(Td, Host_Td, SizeofTd, cudaMemcpyHostToDevice);
	cudaMemcpy(KParameter, Host_KParameter, SizeofKParameter, cudaMemcpyHostToDevice);
	//���ʧ���

	int* d_GBest;
	Particle* d_Group;
	float* d_GBestFitness;

	size_t SizeofGroup = NumofParticle*sizeof(Particle);
	size_t SizeofGBestFitness = sizeof(float);
	size_t SizeofGBest = NumofLight*sizeof(int);


	cudaMalloc(&d_Group, SizeofGroup);
	cudaMalloc(&d_GBestFitness, SizeofGBestFitness);
	cudaMalloc(&d_GBest, SizeofGBest);

	///////////////////////////////////////////////////////
	StartTime = clock(); //�}�l�ɶ�


	/********************************�H�����ͪ�l��*****************************************/
	int* RandomVelocity = new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ProduceParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)
		{
			int* j = new int[NumofLight];
			//�ɤl��l��m�H�����ͤ覡
			for (int h = 0; h<NumofLight; h++)
			{
				j[h] = RandomNumber(0, NumofInstruction - 1);
				RandomVelocity[h] = (RandomNumber(0, 1) ? -1 : 1)*RandomNumber(0, VelocityLimit);
			}
			GroupPointer = new Particle(j, RandomVelocity);
			GroupPointer->CheckPosition();
			if ((GroupPointer->GetPBestFitness())<GBestFitness[GroupNo])
			{
				GBestFitness[GroupNo] = GroupPointer->GetPBestFitness();
				for (int j = 0; j<NumofLight; j++)
				{
					GBest[GroupNo][j] = *(GroupPointer->GetPBest() + j);
				}
			}
			GroupInfo[GroupNo][ProduceParticleNo] = *GroupPointer;  //�h�����
			delete GroupPointer;                  //����Ŷ�
			ProduceParticleNo++;
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
			//////////////CUDA
			cudaMemcpy(d_Group, GroupInfo[GroupNo], SizeofGroup, cudaMemcpyHostToDevice);
			cudaMemcpy(d_GBestFitness, &GBestFitness[GroupNo], SizeofGBestFitness, cudaMemcpyHostToDevice);
			cudaMemcpy(d_GBest, GBest[GroupNo], SizeofGBest, cudaMemcpyHostToDevice);


			PSOKernel << <NumofBlocks, ThreadPerBlocks >> >(d_Group, d_GBest, d_GBestFitness, GroupCircle, I, Td, KParameter, NumofLight, NumofInstruction, NumofTd, NumofParticle);


			cudaMemcpy(GroupInfo[GroupNo], d_Group, SizeofGroup, cudaMemcpyDeviceToHost);
			cudaMemcpy(&GBestFitness[GroupNo], d_GBestFitness, SizeofGBestFitness, cudaMemcpyDeviceToHost);
			cudaMemcpy(GBest[GroupNo], d_GBest, SizeofGBest, cudaMemcpyDeviceToHost);
			///////CODA over///////

			//similarity
			float sim_ave = 0;
			for (int i = 0; i < NumofParticle; i++)
			{
				sim_ave += similarity(GBest[GroupNo], GroupInfo[GroupNo][i].GetPBest(), NumofLight);
			}
			sim_ave = sim_ave / NumofParticle;

			//cout << "Group" << GroupNo + 1 << ", �{�b�������ۦ���:" << sim_ave << endl;
			total_ave_sim += sim_ave;
		}


		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			if (GBestFitness[GroupNo]<TotalGBestFitness || TotalGBestFitness == ImpossibleResume)
			{
				for (int p = 0; p<NumofLight; p++)
				{
					TotalGBest[p] = GBest[GroupNo][p];
				}
				TotalGBestFitness = GBestFitness[GroupNo];
			}
		}

		//�{���D�X���Ѥwí�w
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
			RepeatValue = TotalGBest[NumofLight];
		else if (RepeatValue - TotalGBest[NumofLight] < TotalGBest[NumofLight] * DifferentRange && RepeatValue != -1)
		{
			cout << "�~�j�馸��: " << GroupCircle << endl;
			break;
		}

		//�����������s�A���s�զX
		int RandomCPU;
		int RandomParticle1, RandomParticle2;
		Particle TempMemory;
		for (int j = 0; j<NumofGroup; j++)
		{
			for (int i = 0; i<NumofRedistribute; i++)
			{
				RandomCPU = RandomNumber(0, NumofGroup - 1);
				while (RandomCPU == j)
				{
					RandomCPU = RandomNumber(0, NumofGroup - 1);
				}
				RandomParticle1 = RandomNumber(0, NumofParticle - 1);
				RandomParticle2 = RandomNumber(0, NumofParticle - 1);
				//
				TempMemory = GroupInfo[j][RandomParticle1];
				GroupInfo[j][RandomParticle1] = GroupInfo[RandomCPU][RandomParticle2];
				GroupInfo[RandomCPU][RandomParticle2] = TempMemory;
			}
		}
		//��������//

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			for (int i = 0; i<NumofLight; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			GBestFitness[GroupNo] = TotalGBestFitness;
		}
	}
	EndTime = clock();//�����ɶ�


	/***PSO���浲���A�[�ݵ��G***
	if (TotalGBestFitness != ImpossibleResume)
	{
	//cout<<"PSO�j�M"<<endl
	cout << "���O:";                                                                                 //fw<<"PSO�j�M"<<endl<<"���O:";
	for (int i = 0; i<NumofLight; i++)
	{
	cout << TotalGBest[i] << ",";                                                                            //fw<<GBest[i];
	}
	cout << "�ӯ�:" << TotalGBestFitness << endl;                                                  //fw<<endl<<"�ӯ�:"<<GBest[NumofLight]<<endl;
	}
	else
	{
	cout << "�L��" << endl;
	}
	cout << "�O��:" << double(EndTime - StartTime) / CLOCKS_PER_SEC << "��" << endl << endl;
	/*****************************/

	cudaFree(d_GBest);
	cudaFree(d_Group);
	cudaFree(d_GBestFitness);

	cudaFree(I);
	cudaFree(Td);
	cudaFree(KParameter);

	/*****write file*****/
	filename = "The_Best_Solution.txt";//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\The_Best_Solution.txt
	fw.open(filename, ios::out);//�}���ɮ�
	if (!fw){//�p�G�}���ɮץ��ѡAfw��0�F���\�Afw���D0
		cout << "Fail to open file: " << filename << endl;
	}
	//cout << "�O����O:" << endl;
	if (TotalGBestFitness != ImpossibleResume)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			//cout << TotalGBest[i] << ", ";
			fw << TotalGBest[i] << endl;
		}
		//cout << endl;
		//cout << "�ӯ�:" << TotalGBestFitness << endl;
		fw << TotalGBestFitness << endl;
	}
	else
	{
		//cout << "�L��" << endl;
		fw << "�L��" << endl;
	}
	//cout << "�Ӯ�: " << double(EndTime - StartTime) / CLOCKS_PER_SEC << endl;
	fw << double(EndTime - StartTime) / CLOCKS_PER_SEC << endl;
	fw.close();//�����ɮ�

	//cout << "���饭�����X����50�ɪ��ۦ���: " << total_ave_sim / (NumofGroupCircle*NumofGroup) << endl;

	//����Ŷ�             //GroupInfo �D�ʺA�}�C�A���Φۦ�R��  //�{�b�O�ʺA�F
	//for (int i = 0; i < NumofGroup; i++)
	//{
	//for (int j = 0; j < NumofParticle; j++)
	//{
	delete GroupInfo;
	//}
	//}
	//
	//system("pause");
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

__device__ int device_RandomNumber(unsigned int thread_id,int MinValue, int MaxValue)   //Both MinValue and MaxValue are included
{
	unsigned int seed = thread_id;
	curandState s;
	// seed a random number generator 
	curand_init(seed, 0, 0, &s);

	int R = ((int)curand_uniform(&s) % (MaxValue - MinValue + 1)) + MinValue;
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

/***********************************************
�p��ۦ���
/**********************************************/
float similarity(int* A, int* B, int Length)
{
	//for case 1
	float variation = 0, innerProduct = 0;
	float Similarity;

	//Euclidean distance(�V�p�V��)
	for (int i = 0; i < Length; i++)
	{
		innerProduct += pow(A[i], 2);
		variation += pow(A[i] - B[i], 2);
	}
	Similarity = sqrt(variation) / sqrt(innerProduct);

	//cout <<"variation:"<< fixed << setprecision(5) << Similarity << ", ";
	return Similarity;
}


__device__ float device_similarity(int* A, int* B, int Length)
{
	//for case 1
	float variation = 0, innerProduct = 0;
	float Similarity;

	//Euclidean distance(�V�p�V��)
	for (int i = 0; i < Length; i++)
	{
		innerProduct += A[i] * A[i];
		variation += (A[i] - B[i]) * (A[i] - B[i]);
	}
	Similarity = sqrt(variation) / sqrt(innerProduct);

	//cout <<"variation:"<< fixed << setprecision(5) << Similarity << ", ";
	return Similarity;
}