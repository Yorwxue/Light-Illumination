
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//�ϥΪ̩w�q
#define NumofLight 12        //�O���Ӽ�
#define NumofInstruction 16 //�G�׫��O�Ӽ�
#define NumofTd 5          //�ؼ��I�Ӽ�

#define NumofGroup 10  //���s�ƥءA�W�[���ѼƯണ�ɮį�


#define NumofParticle 100       //��NumofParticle�ոѧ@��PSO�����(�ɤl��) �A�C��particle�j�p��200
#define NumofRedistribute 50



#define OtherLight 0        //�B�~�������Ӷ�

#define ThreadX 4  ///�\��y��X��thread
#define ThreadY 4  ///�\��y��Y��thread

#define BlockX  100 ///�\��y��X��block 
#define BlockY  100 ///�\��y��Y��block

dim3 ThreadPerBlocks(ThreadX, ThreadY);//threads�bblock�����\��覡�A�y��X�PY
dim3 NumofBlocks(BlockX / ThreadPerBlocks.x, BlockY / ThreadPerBlocks.y);//block���\��覡�A�y��X�PY
//float OtherLightT[OtherLight]={0};        //�B�~�������ӫ�


//���ե�
int RunTime = 1;  //�]�h�֦�
float RealBest = -1;//246.34;


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

float RepeatValue = -1;
float DifferentRange = 0.0005;

//__shared__ int CircleofPSO;  
//__device__ int VelocityLimit;
//__shared__ int NumofGroupCircle;

/*
void Parameter()
{
CircleofPSO=NumofParticle*(0.02);      //PSO���檺�`�����ơA�ɤl�ƪ�2%
if(CircleofPSO<10)   //�������ƹL�ַ|�y���L��
CircleofPSO=10;

VelocityLimit=NumofInstruction*(0.5);  //���ʳt�פW���A���ʪŶ���50%
if(VelocityLimit<3)   //���ʭ�����Y�A�|�ɭP���鰱���
VelocityLimit=3;

NumofGroupCircle=NumofGroup;    //�U�B�z������y�һݡA�����B�z���Ӽ�
}
*/
//
//__shared__ float I[NumofInstruction][2];   //(�ӯ�,�G��)    //__shared__���O�A�bdeviceŪ���Ȯɷ|�y�����~
//__shared__ float Td[NumofTd];
//__shared__ float KParameter[NumofTd][NumofLight+OtherLight];


float** Host_I;   //(�ӯ�,�G��)
float* Host_Td;
float** Host_KParameter;

int CheckT(float T, int MatrixRow, int i, int j, int k);
float PSO();
int RandomNumber(int MinValue, int MaxValue);
void* new2d(int h, int w, int size);

class Particle  //�C�Ӳɤl���򥻳]�w�A�ܰʽե����O����
{
public:
	int* PBest;  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,����̧C�ӯ�)
	float PBestFitness;
	int GBest[NumofLight];
	float GBestFitness;
	int Position[NumofLight];  //�O������A��U���O�����ե����O <----> Position   //�䤤�O���Ƥ]�Osolution�����סA�ĴX���O��ܲĴX��
	int Velocity[NumofLight];             //�ե����O�W�ɩΤU�� <----> Velocity

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
	__host__  void SetGBest(int* GBestPointer, float NewGBestFitness)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			GBest[i] = *(GBestPointer + i);
		}
		GBestFitness = NewGBestFitness;
	}

	__device__  int* d_GetPBest(){ return PBest; }
	__device__  float d_GetPBestFitness(){ return PBestFitness; }
	/***********
	__device__  void SetGBest(int* GBestPointer,float NewGBestFitness)
	{
	for(int i=0;i<NumofLight;i++)
	{
	GBest[i]=*(GBestPointer+i);
	}
	GBestFitness=NewGBestFitness;
	}
	__device__  void UpdatePosition()
	{
	for(int i=0;i<NumofLight;i++)//i�O����
	{
	Position[i]=Position[i]+Velocity[i];
	if(Position[i]<0)                                //�ե����O���W�U��
	{
	Position[i]=0;
	}
	else if(Position[i]>(NumofInstruction-1))             //�ե����O���W�U��
	{
	Position[i]=(NumofInstruction-1);
	}
	}
	}
	__device__  void UpdateVelocity(float C0,float C1,float C2)
	{
	for(int i=0;i<NumofLight;i++)//i�O����
	{
	Velocity[i]=(C0*Velocity[i]+C1*(PBest[i]-Position[i])+C2*(GBest[i]-Position[i]));
	if(Velocity[i]<(-VelocityLimit))                //�t�פW�U��
	{
	Velocity[i]=(-VelocityLimit);
	}
	if(Velocity[i]>VelocityLimit)                   //�t�פW�U��
	{
	Velocity[i]=VelocityLimit;
	}
	}
	}
	__device__  void d_CheckPosition(float* I,float* Td,float* KParameter)  //�s�b���D
	{
	float Trow=0;  //�ھڥثe�ɤl����m�Һ�X�Ӫ�T(�ӫ�)��

	for(int i=0;i<NumofTd;i++)
	{
	for(int j=0;j<(NumofLight+OtherLight);j++)
	{
	Trow+=*(KParameter+i*(NumofLight+OtherLight)+j) * *(I+2*Position[j]+1);//KParameter*I;
	}
	if(Trow<*(Td+i))
	{
	break;
	}
	Trow=0;
	if(i==(NumofTd-1))
	{
	d_ComputingFitness(I);
	}
	}
	}
	__device__  void d_ComputingFitness(float* I)
	{
	float Fitness=0;
	for(int i=0;i<NumofLight;i++)
	{
	Fitness += *(I+2*Position[i]+0);
	}

	if(Fitness<PBestFitness)
	{
	for(int j=0;j<NumofLight;j++)
	{
	PBest[j] = Position[j];
	}
	PBestFitness=Fitness;
	}
	}
	********/
	__device__  void DeviceFunction(int* GBestPointer, float NewGBestFitness, float C0, float C1, float C2, float* I, float* Td, float* KParameter)
	{
		//SetGBest
		GBestFitness = NewGBestFitness;

		//SetGBest & UpdateVelocity & Position 
		for (int i = 0; i<NumofLight; i++)//i�O����
		{
			//SetGBest
			GBest[i] = *(GBestPointer + i);

			//UpdateVelocity
			Velocity[i] = (C0*Velocity[i] + C1*(PBest[i] - Position[i]) + C2*(GBest[i] - Position[i]));
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
				Trow += *(KParameter + i*(NumofLight + OtherLight) + j) * *(I + 2 * Position[j] + 1);//KParameter*I;  
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
	//test
	__device__  void testSetPBest(int Value)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			PBest[i] = Value;
		}
	}
};

Particle::Particle(int* PPointer, int* VPointer)
{
	PBest = new int[NumofLight];
	for (int i = 0; i<NumofLight; i++)
	{
		Position[i] = *(PPointer + i);
		Velocity[i] = *(VPointer + i);
		PBest[i] = Position[i];
	}
	PBestFitness = ImpossibleResume;
}//PSO���O����


__global__ void PSOKernel(Particle* Group, int* GBest, float* GBestFitness, int GroupCircle, float* I, float* Td, float* KParameter)
{
	int CoordinatesX = blockIdx.x * blockDim.x + threadIdx.x;
	int CoordinatesY = blockIdx.y * blockDim.y + threadIdx.y;
	int ParticleNoinGroup = CoordinatesY*(ThreadX*BlockX) + CoordinatesX;//Y*(row size)+(X)
	if (ParticleNoinGroup < NumofParticle && CoordinatesX < (ThreadX*BlockX) && CoordinatesY < (ThreadY*BlockY))
		for (int circle = 1; circle < CircleofPSO; circle++)
		{
		float C0 = InitialW - (InitialW - FinalW)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
		float C1 = InitialCOne - (InitialCOne - FinalCOne)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
		float C2 = InitialCTwo + (FinalCTwo - InitialCTwo)*((circle + GroupCircle) / (NumofGroupCircle + CircleofPSO - 1));
		/****
		//��sGBest
		(Group+ParticleNoinGroup)->SetGBest(GBest,*(GBestFitness));
		//��s�t��
		(Group+ParticleNoinGroup)->UpdateVelocity(C0,C1,C2);
		//��s��m
		(Group+ParticleNoinGroup)->UpdatePosition();
		//�ˬd�X��
		(Group+ParticleNoinGroup)->d_CheckPosition(I,Td,KParameter);
		/*****/
		(Group + ParticleNoinGroup)->DeviceFunction(GBest, *(GBestFitness), C0, C1, C2, I, Td, KParameter);
		/*****/
		//(Group+ParticleNoinGroup)->DelayFunction(); 
		/****
		__syncthreads();

		for(int e=0;e<NumofParticle;e++)
		{
		if(((Group+e)->d_GetPBestFitness())< *(GBestFitness) )
		{
		for(int j=0;j<NumofLight;j++)
		{
		*(GBest+j)=*(((Group+e)->d_GetPBest())+j);
		}
		*(GBestFitness)=(Group+e)->d_GetPBestFitness();
		}
		}
		/****/
		if (((Group + ParticleNoinGroup)->d_GetPBestFitness()) < *(GBestFitness))
		{
			for (int j = 0; j < NumofLight; j++)
			{
				*(GBest + j) = *(((Group + ParticleNoinGroup)->d_GetPBest()) + j);
			}
			*(GBestFitness) = (Group + ParticleNoinGroup)->d_GetPBestFitness();
		}
		/****/
		}
	if (ParticleNoinGroup == 0)  //for test
	{
		(Group + ParticleNoinGroup)->testSetPBest(-1);
	}
}


int main()
{
	/***********************Ū���ɮ�**********************/
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�o�ONumofInstruction
		Host_I = (float **)new2d(NumofInstruction, 2, sizeof(float));
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_I[i][j]);        //�r����Ʀr
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
		Host_Td = new float[NumofTd];
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf(FileInput, "%f", &Host_Td[i]);        //�r����Ʀr
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}��
		fr.getline(FileInput, sizeof(FileInput), ',');//�x�}�e
		Host_KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_KParameter[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();
	//���Ū������

	float PSOResult;
	int NumberofBest = 0;
	int Numberof5 = 0;
	float ResultTime = 0;
	int Others = 0;
	for (int i = 0; i<RunTime; i++)
	{
		PSOResult = PSO();
		if (PSOResult<RealBest + 0.001)
			NumberofBest++;
		else if (PSOResult<(RealBest + RealBest*0.05) && PSOResult>RealBest)
			Numberof5++;
		else
			Others++;
		//if(i!=0) 
		//ResultTime+=PSOResult;
	}
	/*********************���ק�PSO�禡���^�ǭ�************************/
	cout << "Best:" << NumberofBest << "��" << endl;
	cout << "5%:" << Numberof5 << "��" << endl;
	cout << "Others:" << Others << "��" << endl;
	//cout<<"�����O��"<< ResultTime/RunTime<<"��"<<endl;

	/***����Ŷ�***/
	delete[] Host_I;
	delete[] Host_Td;
	delete[] Host_KParameter;

	system("pause");
	return 0;
}




float PSO()
{
	fstream fw;
	string filename;

	clock_t StartTime, EndTime;
	int* TotalGBest=new int[NumofLight];
	float TotalGBestFitness = ImpossibleResume;

	int** GBest=(int**)new2d(NumofGroup,NumofLight,sizeof(int));
	float* GBestFitness=new float[NumofGroup];
	for (int i = 0; i<NumofGroup; i++)
	{
		GBestFitness[i] = ImpossibleResume;
	}

	int ProduceParticleNo = 0;

	Particle* GroupPointer;
	Particle** GroupInfo=(Particle**)new2d(NumofGroup,NumofParticle,sizeof(Particle));

	/*********���������t�ɨϥ�************
	Particle TempGroup[NumofGroup][NumofParticle];
	int TempGroupNo[NumofGroup];  //�p��洫�ɤl�Ϊ��Ȧs�s�O�_�w��
	for (int i = 0; i<NumofGroup; i++)
	{
	TempGroupNo[i] = 0;
	}
	/**************************************/
	////////for CUDA//////////////////////

	//GPU�ǳƪŶ�
	float* I;   //(�ӯ�,�G��)
	float* Td;
	float* KParameter;

	size_t SizeofI = NumofInstruction * 2 * sizeof(float);
	cudaMalloc(&I, SizeofI);

	size_t SizeofTd = NumofTd*sizeof(float);
	cudaMalloc(&Td, SizeofTd);

	size_t SizeofKParameter = NumofTd*(NumofLight + OtherLight)*sizeof(float);
	cudaMalloc(&KParameter, SizeofKParameter);

	//�N��Ʋ���GPU
	cudaMemcpy(I, Host_I + NumofInstruction, SizeofI, cudaMemcpyHostToDevice);
	cudaMemcpy(Td, Host_Td, SizeofTd, cudaMemcpyHostToDevice);
	cudaMemcpy(KParameter, Host_KParameter + NumofTd, SizeofKParameter, cudaMemcpyHostToDevice);
	//���ʧ���

	int* d_GBest;
	Particle* d_Group;
	float* d_GBestFitness;

	size_t SizeofGroup = NumofParticle*sizeof(Particle);
	size_t SizeofGBestFitness = sizeof(float);
	size_t SizeofGBest = NumofLight*sizeof(int);

	size_t SizeofPBest = sizeof(int[NumofLight]);


	cudaMalloc(&d_Group, SizeofGroup);
	
	for (int i = 0; i < NumofParticle; i++)
	{
		cudaMalloc(&d_Group[i].PBest, SizeofPBest);    //CUDA�S��PBest�Ŷ��s����
	}

	cudaMalloc(&d_GBestFitness, SizeofGBestFitness);
	cudaMalloc(&d_GBest, SizeofGBest);

	///////////////////////////////////////////////////////
	StartTime = clock(); //�}�l�ɶ�


	/********************************�H�����ͪ�l��*****************************************/
	int* RandomVelocity=new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ProduceParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)
		{
			int* j=new int[NumofLight];
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
	//��l��m���ͧ���

	/*****��z�ɤl�b�O���餤����m(�[�t�N��ƶǤJGPU���t��)*****
	for(int GroupNo=0;GroupNo<NumofGroup;GroupNo++)
	{
	for(int i=1;i<NumofParticle;i++)
	{
	GroupInfo[GroupNo][i]=*Group[GroupNo][i];  //�h�����
	delete Group[GroupNo][i];                  //����Ŷ�
	}
	}
	/*****��z����*****/

	//��Ȳ��ͧ���
	//�P�ɡA�Ĥ@��PSO�����A�ѤUCircleofPSO-1��
	//�s�W����B�z�A�b�̥~�h
	/*******************************************************************************************************/
	for (int GroupClicle = 0; GroupClicle<NumofGroupCircle; GroupClicle++)
	{
		cout << "�e:";
		for (int i = 0; i < NumofLight; i++)
		{
			cout << GroupInfo[0][0].PBest[i] << ",";
		}cout << GroupInfo[0][0].PBestFitness << endl;
		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			//////////////CUDA
			for (int i = 0; i < NumofParticle;i++)
			{
				cudaMemcpy(d_Group[i].GBest, GroupInfo[GroupNo][i].GBest, sizeof(int[NumofLight]), cudaMemcpyHostToDevice);
				cudaMemcpy(&d_Group[i].GBestFitness, &GroupInfo[GroupNo][i].GBestFitness, sizeof(float), cudaMemcpyHostToDevice);
				//cudaMemcpy(d_Group[i].PBest, GroupInfo[GroupNo][i].PBest, sizeof(int[NumofLight]), cudaMemcpyHostToDevice);  
				
				cudaMemcpy(d_Group[i].PBest, GroupInfo[GroupNo][i].PBest, sizeof(int[NumofLight]), cudaMemcpyHostToDevice);

				cudaMemcpy(&d_Group[i].PBestFitness, &GroupInfo[GroupNo][i].PBestFitness, sizeof(float), cudaMemcpyHostToDevice);
				cudaMemcpy(d_Group[i].Position, GroupInfo[GroupNo][i].Position, sizeof(int[NumofLight]), cudaMemcpyHostToDevice);
				cudaMemcpy(d_Group[i].Velocity, GroupInfo[GroupNo][i].Velocity, sizeof(int[NumofLight]), cudaMemcpyHostToDevice);
			}
			//cudaMemcpy(d_Group, GroupInfo[GroupNo], SizeofGroup, cudaMemcpyHostToDevice);
			cudaMemcpy(d_GBestFitness, &GBestFitness[GroupNo], SizeofGBestFitness, cudaMemcpyHostToDevice);
			cudaMemcpy(d_GBest, GBest[GroupNo], SizeofGBest, cudaMemcpyHostToDevice);


			PSOKernel <<< NumofBlocks, ThreadPerBlocks >>>(d_Group, d_GBest, d_GBestFitness, GroupClicle, I, Td, KParameter);


			for (int i = 0; i < NumofParticle; i++)
			{
				cudaMemcpy(d_Group[i].GBest, GroupInfo[GroupNo][i].GBest, sizeof(int[NumofLight]), cudaMemcpyDeviceToHost);
				cudaMemcpy(&d_Group[i].GBestFitness, &GroupInfo[GroupNo][i].GBestFitness, sizeof(float), cudaMemcpyDeviceToHost);

				cudaMemcpy(d_Group[i].PBest, GroupInfo[GroupNo][i].PBest, sizeof(int[NumofLight]), cudaMemcpyDeviceToHost);

				cudaMemcpy(&d_Group[i].PBestFitness, &GroupInfo[GroupNo][i].PBestFitness, sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(d_Group[i].Position, GroupInfo[GroupNo][i].Position, sizeof(int[NumofLight]), cudaMemcpyDeviceToHost);
				cudaMemcpy(d_Group[i].Velocity, GroupInfo[GroupNo][i].Velocity, sizeof(int[NumofLight]), cudaMemcpyDeviceToHost);
			}
			//cudaMemcpy(GroupInfo[GroupNo], d_Group, SizeofGroup, cudaMemcpyDeviceToHost);
			cudaMemcpy(&GBestFitness[GroupNo], d_GBestFitness, SizeofGBestFitness, cudaMemcpyDeviceToHost);
			cudaMemcpy(GBest[GroupNo], d_GBest, SizeofGBest, cudaMemcpyDeviceToHost);
			///////CODA over///////
		}
		cout << "��:";
		for (int i = 0; i < NumofLight; i++)
		{
			cout << GroupInfo[0][0].PBest[i] << ",";
		}cout << GroupInfo[0][0].PBestFitness << endl;


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
		/*
		//�����覡
		//��ҭ���
		int RedistributeList[NumofGroup];
		int RedistributeParticleList[NumofGroup][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;
		Particle TempMemory0[NumofRedistribute];
		Particle TempMemory1[NumofRedistribute];

		//�إ߸s��list
		for(int i=0;i<NumofGroup;i++)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		for(int j=0;j<i;j++)
		{
		if(RedistributeList[j]==RandomGroup)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		j=0;
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
		j=0;
		}
		}
		RedistributeParticleList[k][i]=RandomParticle;
		}
		}

		//�}�l�洫
		for(int i=0,k=0;k<NumofGroup;k++) //i�O�Ҧb���s,k�O�ĴX���@�s�������洫
		{
		if(k==0)
		{
		for(int j=0;j<NumofRedistribute;j++)
		{
		TempMemory0[j] = GroupInfo[i][ RedistributeParticleList[i][j] ];
		}
		}
		if((k%2)==0)
		{
		for(int j=0;j<NumofRedistribute;j++)
		{
		TempMemory1[j] = GroupInfo[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ];
		GroupInfo[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempMemory0[j];
		}
		}
		else
		{
		for(int j=0;j<NumofRedistribute;j++)
		{
		TempMemory0[j] = GroupInfo[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ];
		GroupInfo[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempMemory1[j];
		}
		}
		i=RedistributeList[i];
		}

		/********************************/
		/************��������************
		int RandomGroup;
		for (int i = 0; i<NumofGroup; i++)
		{
		for (int j = 0; j<NumofParticle; j++)
		{
		RandomGroup = RandomNumber(0, NumofGroup - 1);
		while (TempGroupNo[RandomGroup] == NumofParticle)//�Ȧs�s�w��
		{
		RandomGroup = RandomNumber(0, NumofGroup - 1);
		}
		TempGroup[RandomGroup][TempGroupNo[RandomGroup]++] = GroupInfo[i][j];
		}
		}
		for (int i = 0; i<NumofGroup; i++)
		{
		for (int j = 0; j<NumofParticle; j++)
		{
		GroupInfo[i][j] = TempGroup[i][j];
		}
		TempGroupNo[i] = 0;
		}
		/*******************************/

		/********************************
		//////////////�����������s�A���s�զX
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
		///////////*///��������

		/*�M�䭫���᪺�U�s��gbest*/
		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			for (int i = 0; i<NumofParticle; i++)
			{
				GBestFitness[GroupNo] = ImpossibleResume;
			}
			for (int i = 0; i<NumofParticle; i++)
			{
				if (GBestFitness[GroupNo] > GroupInfo[GroupNo][i].GetPBestFitness() || GBestFitness[GroupNo] == ImpossibleResume)
				{
					GBest[GroupNo] = GroupInfo[GroupNo][i].GetPBest();
					GBestFitness[GroupNo] = GroupInfo[GroupNo][i].GetPBestFitness();
				}
			}
		}

		/****�M�䧹��****/

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			/*�U�s���s��X��gbest�g�^*/
			for(int i=0;i<NumofParticle;i++)
			{
			    GroupInfo[GroupNo][i].SetGBest(GBest[GroupNo],GBestFitness[GroupNo]);
			}
			/******/
			/*�N��X��Gbest�g�^���U�s*
			for (int i = 0; i<NumofLight; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			GBestFitness[GroupNo] = TotalGBestFitness;
			/*************************/
		}

		/*�{���D�X���Ѥwí�w*
		if (RepeatValue == -1 && TotalGBestFitness != ImpossibleResume)
			RepeatValue = TotalGBestFitness;
		else if (RepeatValue - TotalGBestFitness < TotalGBestFitness * DifferentRange && RepeatValue != -1)
			break;
		/****/

	}
	EndTime = clock();//�����ɶ�


	/***PSO���浲���A�[�ݵ��G***/
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
	//fw.close();
	/*****************************/

	cudaFree(d_GBest);
	cudaFree(d_Group);
	cudaFree(d_GBestFitness);

	cudaFree(I);
	cudaFree(Td);
	cudaFree(KParameter);

	/*****write file*****/
	filename = "DPSOResult.txt";
	fw.open(filename, ios::out);//�}���ɮ�
	if (!fw){//�p�G�}���ɮץ��ѡAfw��0�F���\�Afw���D0
		cout << "Fail to open file: " << filename << endl;
	}
	if (TotalGBestFitness != ImpossibleResume)
	{
		fw << "���O:";                                                                                 //fw<<"PSO�j�M"<<endl<<"���O:"; 
		for (int i = 0; i<NumofLight; i++)
		{
			fw << TotalGBest[i] << ",";                                                                            //fw<<GBest[i];  
		}
		fw << "�ӯ�:" << TotalGBestFitness << endl;                                                  //fw<<endl<<"�ӯ�:"<<GBest[NumofLight]<<endl;    
	}
	else
	{
		fw << "�L��" << endl;
	}
	fw << "�O��:" << double(EndTime - StartTime) / CLOCKS_PER_SEC << "��" << endl << endl;
	fw.close();//�����ɮ�

	//����Ŷ�             //GroupInfozo �D�ʺA�}�C�A���Φۦ�R��
	//for (int i = 0; i < NumofGroup; i++)
	//{
	//	for (int j = 0; j < NumofParticle; j++)
	//	{
	//		delete[] GroupInfo;
	//	}
	//}
	//

	return TotalGBestFitness;
	//return (float(EndTime-StartTime)/CLOCKS_PER_SEC);
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

/********************���ͰʺA2���}�C****************************/
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