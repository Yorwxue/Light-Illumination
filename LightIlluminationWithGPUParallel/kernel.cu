
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//使用者定義
#define NumofLight 12        //燈的個數
#define NumofInstruction 16 //亮度指令個數
#define NumofTd 5          //目標點個數

#define NumofGroup 10  //分群數目，增加此參數能提升效能


#define NumofParticle 100       //取NumofParticle組解作為PSO的初值(粒子數) ，每個particle大小為200
#define NumofRedistribute 50



#define OtherLight 0        //額外光源為太陽

#define ThreadX 4  ///擺放座標X的thread
#define ThreadY 4  ///擺放座標Y的thread

#define BlockX  100 ///擺放座標X的block 
#define BlockY  100 ///擺放座標Y的block

dim3 ThreadPerBlocks(ThreadX, ThreadY);//threads在block中的擺放方式，座標X與Y
dim3 NumofBlocks(BlockX / ThreadPerBlocks.x, BlockY / ThreadPerBlocks.y);//block的擺放方式，座標X與Y
//float OtherLightT[OtherLight]={0};        //額外光源的照度


//測試用
int RunTime = 1;  //跑多少次
float RealBest = -1;//246.34;


//程式自定義
#define ImpossibleResume 99999
#define InitialW 0.9             //慣性係數，由大到小
#define FinalW 0.4
#define InitialCOne 2.5          //本身的C係數，由大到小，因為搜尋初期重視個體
#define FinalCOne 1
#define InitialCTwo  1       //群體的C係數，由小到大，因為搜尋末期重視群體
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
CircleofPSO=NumofParticle*(0.02);      //PSO執行的循環次數，粒子數的2%
if(CircleofPSO<10)   //探索次數過少會造成無解
CircleofPSO=10;

VelocityLimit=NumofInstruction*(0.5);  //移動速度上限，移動空間的50%
if(VelocityLimit<3)   //移動限制太嚴，會導致個體停止不動
VelocityLimit=3;

NumofGroupCircle=NumofGroup;    //各處理器間交流所需，約為處理器個數
}
*/
//
//__shared__ float I[NumofInstruction][2];   //(耗能,亮度)    //__shared__指令，在device讀取值時會造成錯誤
//__shared__ float Td[NumofTd];
//__shared__ float KParameter[NumofTd][NumofLight+OtherLight];


float** Host_I;   //(耗能,亮度)
float* Host_Td;
float** Host_KParameter;

int CheckT(float T, int MatrixRow, int i, int j, int k);
float PSO();
int RandomNumber(int MinValue, int MaxValue);
void* new2d(int h, int w, int size);

class Particle  //每個粒子的基本設定，變動調光指令版本
{
public:
	int* PBest;  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,個體最低耗能)
	float PBestFitness;
	int GBest[NumofLight];
	float GBestFitness;
	int Position[NumofLight];  //記錄個體，其各盞燈光的調光指令 <----> Position   //其中燈光數也是solution的維度，第幾盞燈表示第幾維
	int Velocity[NumofLight];             //調光指令上升或下降 <----> Velocity

public:
	Particle(int* PPointer, int* VPointer);
	Particle(){ PBestFitness = ImpossibleResume; }
	__host__  int* GetPBest(){ return PBest; }
	__host__  float GetPBestFitness(){ return PBestFitness; }
	__host__  void CheckPosition()
	{
		float Trow = 0;  //根據目前粒子的位置所算出來的T(照度)值
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
	for(int i=0;i<NumofLight;i++)//i是維度
	{
	Position[i]=Position[i]+Velocity[i];
	if(Position[i]<0)                                //調光指令的上下限
	{
	Position[i]=0;
	}
	else if(Position[i]>(NumofInstruction-1))             //調光指令的上下限
	{
	Position[i]=(NumofInstruction-1);
	}
	}
	}
	__device__  void UpdateVelocity(float C0,float C1,float C2)
	{
	for(int i=0;i<NumofLight;i++)//i是維度
	{
	Velocity[i]=(C0*Velocity[i]+C1*(PBest[i]-Position[i])+C2*(GBest[i]-Position[i]));
	if(Velocity[i]<(-VelocityLimit))                //速度上下限
	{
	Velocity[i]=(-VelocityLimit);
	}
	if(Velocity[i]>VelocityLimit)                   //速度上下限
	{
	Velocity[i]=VelocityLimit;
	}
	}
	}
	__device__  void d_CheckPosition(float* I,float* Td,float* KParameter)  //存在問題
	{
	float Trow=0;  //根據目前粒子的位置所算出來的T(照度)值

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
		for (int i = 0; i<NumofLight; i++)//i是維度
		{
			//SetGBest
			GBest[i] = *(GBestPointer + i);

			//UpdateVelocity
			Velocity[i] = (C0*Velocity[i] + C1*(PBest[i] - Position[i]) + C2*(GBest[i] - Position[i]));
			if (Velocity[i]<(-VelocityLimit))                //速度上下限
			{
				Velocity[i] = (-VelocityLimit);
			}
			else if (Velocity[i]>VelocityLimit)                   //速度上下限
			{
				Velocity[i] = VelocityLimit;
			}

			//UpdatePosition 
			Position[i] = Position[i] + Velocity[i];
			if (Position[i]<0)                                //調光指令的上下限
			{
				Position[i] = 0;
			}
			else if (Position[i]>(NumofInstruction - 1))             //調光指令的上下限
			{
				Position[i] = (NumofInstruction - 1);
			}
		}
		//CheckPosition
		float Trow = 0;  //根據目前粒子的位置所算出來的T(照度)值     

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
}//PSO類別結束


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
		//更新GBest
		(Group+ParticleNoinGroup)->SetGBest(GBest,*(GBestFitness));
		//更新速度
		(Group+ParticleNoinGroup)->UpdateVelocity(C0,C1,C2);
		//更新位置
		(Group+ParticleNoinGroup)->UpdatePosition();
		//檢查合格
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
	/***********************讀取檔案**********************/
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofInstruction
		Host_I = (float **)new2d(NumofInstruction, 2, sizeof(float));
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_I[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofTd
		Host_Td = new float[NumofTd];
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf(FileInput, "%f", &Host_Td[i]);        //字串轉數字
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣高
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣寬
		Host_KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_KParameter[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();
	//資料讀取完畢

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
	/*********************須修改PSO函式的回傳值************************/
	cout << "Best:" << NumberofBest << "次" << endl;
	cout << "5%:" << Numberof5 << "次" << endl;
	cout << "Others:" << Others << "次" << endl;
	//cout<<"平均費時"<< ResultTime/RunTime<<"秒"<<endl;

	/***釋放空間***/
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

	/*********全部重分配時使用************
	Particle TempGroup[NumofGroup][NumofParticle];
	int TempGroupNo[NumofGroup];  //計算交換粒子用的暫存群是否已滿
	for (int i = 0; i<NumofGroup; i++)
	{
	TempGroupNo[i] = 0;
	}
	/**************************************/
	////////for CUDA//////////////////////

	//GPU準備空間
	float* I;   //(耗能,亮度)
	float* Td;
	float* KParameter;

	size_t SizeofI = NumofInstruction * 2 * sizeof(float);
	cudaMalloc(&I, SizeofI);

	size_t SizeofTd = NumofTd*sizeof(float);
	cudaMalloc(&Td, SizeofTd);

	size_t SizeofKParameter = NumofTd*(NumofLight + OtherLight)*sizeof(float);
	cudaMalloc(&KParameter, SizeofKParameter);

	//將資料移至GPU
	cudaMemcpy(I, Host_I + NumofInstruction, SizeofI, cudaMemcpyHostToDevice);
	cudaMemcpy(Td, Host_Td, SizeofTd, cudaMemcpyHostToDevice);
	cudaMemcpy(KParameter, Host_KParameter + NumofTd, SizeofKParameter, cudaMemcpyHostToDevice);
	//移動完畢

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
		cudaMalloc(&d_Group[i].PBest, SizeofPBest);    //CUDA沒給PBest空間存放資料
	}

	cudaMalloc(&d_GBestFitness, SizeofGBestFitness);
	cudaMalloc(&d_GBest, SizeofGBest);

	///////////////////////////////////////////////////////
	StartTime = clock(); //開始時間


	/********************************隨機產生初始值*****************************************/
	int* RandomVelocity=new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ProduceParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)
		{
			int* j=new int[NumofLight];
			//粒子初始位置隨機產生方式
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
			GroupInfo[GroupNo][ProduceParticleNo] = *GroupPointer;  //搬移資料
			delete GroupPointer;                  //釋放空間
			ProduceParticleNo++;
		}
	}
	//初始位置產生完畢

	/*****整理粒子在記憶體中的位置(加速將資料傳入GPU的速度)*****
	for(int GroupNo=0;GroupNo<NumofGroup;GroupNo++)
	{
	for(int i=1;i<NumofParticle;i++)
	{
	GroupInfo[GroupNo][i]=*Group[GroupNo][i];  //搬移資料
	delete Group[GroupNo][i];                  //釋放空間
	}
	}
	/*****整理完畢*****/

	//初值產生完畢
	//同時，第一輪PSO結束，剩下CircleofPSO-1次
	//新增平行處理，在最外層
	/*******************************************************************************************************/
	for (int GroupClicle = 0; GroupClicle<NumofGroupCircle; GroupClicle++)
	{
		cout << "前:";
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
		cout << "後:";
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
		//重分方式
		//比例重分
		int RedistributeList[NumofGroup];
		int RedistributeParticleList[NumofGroup][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;
		Particle TempMemory0[NumofRedistribute];
		Particle TempMemory1[NumofRedistribute];

		//建立群的list
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

		//選擇要交換的粒子，建立對照表
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

		//開始交換
		for(int i=0,k=0;k<NumofGroup;k++) //i是所在的群,k是第幾次作群之間的交換
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
		/************全部重分************
		int RandomGroup;
		for (int i = 0; i<NumofGroup; i++)
		{
		for (int j = 0; j<NumofParticle; j++)
		{
		RandomGroup = RandomNumber(0, NumofGroup - 1);
		while (TempGroupNo[RandomGroup] == NumofParticle)//暫存群已滿
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
		//////////////打散部分分群，重新組合
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
		///////////*///重分完畢

		/*尋找重分後的各群的gbest*/
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

		/****尋找完畢****/

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			/*各群重新找出的gbest寫回*/
			for(int i=0;i<NumofParticle;i++)
			{
			    GroupInfo[GroupNo][i].SetGBest(GBest[GroupNo],GBestFitness[GroupNo]);
			}
			/******/
			/*將整合的Gbest寫回給各群*
			for (int i = 0; i<NumofLight; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			GBestFitness[GroupNo] = TotalGBestFitness;
			/*************************/
		}

		/*認為求出的解已穩定*
		if (RepeatValue == -1 && TotalGBestFitness != ImpossibleResume)
			RepeatValue = TotalGBestFitness;
		else if (RepeatValue - TotalGBestFitness < TotalGBestFitness * DifferentRange && RepeatValue != -1)
			break;
		/****/

	}
	EndTime = clock();//結束時間


	/***PSO執行結束，觀看結果***/
	if (TotalGBestFitness != ImpossibleResume)
	{
		//cout<<"PSO搜尋"<<endl
		cout << "指令:";                                                                                 //fw<<"PSO搜尋"<<endl<<"指令:"; 
		for (int i = 0; i<NumofLight; i++)
		{
			cout << TotalGBest[i] << ",";                                                                            //fw<<GBest[i];  
		}
		cout << "耗能:" << TotalGBestFitness << endl;                                                  //fw<<endl<<"耗能:"<<GBest[NumofLight]<<endl;    
	}
	else
	{
		cout << "無解" << endl;
	}
	cout << "費時:" << double(EndTime - StartTime) / CLOCKS_PER_SEC << "秒" << endl << endl;
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
	fw.open(filename, ios::out);//開啟檔案
	if (!fw){//如果開啟檔案失敗，fw為0；成功，fw為非0
		cout << "Fail to open file: " << filename << endl;
	}
	if (TotalGBestFitness != ImpossibleResume)
	{
		fw << "指令:";                                                                                 //fw<<"PSO搜尋"<<endl<<"指令:"; 
		for (int i = 0; i<NumofLight; i++)
		{
			fw << TotalGBest[i] << ",";                                                                            //fw<<GBest[i];  
		}
		fw << "耗能:" << TotalGBestFitness << endl;                                                  //fw<<endl<<"耗能:"<<GBest[NumofLight]<<endl;    
	}
	else
	{
		fw << "無解" << endl;
	}
	fw << "費時:" << double(EndTime - StartTime) / CLOCKS_PER_SEC << "秒" << endl << endl;
	fw.close();//關閉檔案

	//釋放空間             //GroupInfozo 非動態陣列，不用自行刪除
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
給予希望的最大值與最小值，就能產生介於這兩個值中的隨機數，這個隨機亂數也可能等於最大值，或等於最小值。
其亂數值是依照時間所產生，所以重複呼叫時，會給予不同的亂數值(整數)。
***************************************************************************************************/
int RandomNumber(int MinValue, int MaxValue)   //Both MinValue and MaxValue are included
{
	int R = (rand() % (MaxValue - MinValue + 1)) + MinValue;
	return R;
}

/********************產生動態2為陣列****************************/
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