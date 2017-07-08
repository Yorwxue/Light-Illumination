
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

//使用者定義
int NumofLight;        //燈的個數
int NumofInstruction; //亮度指令個數
int NumofTd;          //目標點個數

#define MaxNumofLight 12        //燈的個數
#define MaxNumofInstruction 16 //亮度指令個數
#define MaxNumofTd 5          //目標點個數

int NumofGroup;  //分群數目，增加此參數能提升效能

float RepeatValue = -1;
float DifferentRange;
int NumofParticle;       //取NumofParticle組解作為PSO的初值(粒子數) ，每個particle大小為200
int NumofRedistribute;



#define OtherLight 0        //額外光源為太陽

#define ThreadX 8  ///擺放座標X的thread
#define ThreadY 8  ///擺放座標Y的thread

#define BlockX  128 ///擺放座標X的block 
#define BlockY  128 ///擺放座標Y的block

dim3 ThreadPerBlocks(ThreadX, ThreadY);//threads在block中的擺放方式，座標X與Y
dim3 NumofBlocks(BlockX / ThreadPerBlocks.x, BlockY / ThreadPerBlocks.y);//block的擺放方式，座標X與Y
//float OtherLightT[OtherLight]={0};        //額外光源的照度


//測試用
int RunTime = 1;  //執行多少次
float RealBest = -1;
#define in_cir_stable 2
#define in_cir_threshold 0.5


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



float Host_I[MaxNumofInstruction][2];   //(耗能,亮度)
float* Host_Td;
float Host_KParameter[MaxNumofTd][MaxNumofLight + OtherLight];

int CheckT(float T, int MatrixRow, int i, int j, int k);
void PSO();
int RandomNumber(int MinValue, int MaxValue);
void* new2d(int h, int w, int size);
float similarity(int* A, int* B, int Length);
__device__ float device_similarity(int* A, int* B, int Length);
__device__ int device_RandomNumber(unsigned int thread_id, int MinValue, int MaxValue);

class Particle  //每個粒子的基本設定，變動調光指令版本
{
private:
	int PBest[MaxNumofLight];  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,個體最低耗能)
	float PBestFitness;
	int GBest[MaxNumofLight];
	float GBestFitness;
	int Position[MaxNumofLight];  //記錄個體，其各盞燈光的調光指令 <----> Position   //其中燈光數也是solution的維度，第幾盞燈表示第幾維
	int Velocity[MaxNumofLight];             //調光指令上升或下降 <----> Velocity

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

	__device__  int* d_GetPBest(){ return PBest; }
	__device__  float d_GetPBestFitness(){ return PBestFitness; }
	__device__  void DeviceFunction(int* GBestPointer, float NewGBestFitness, float C0, float C1, float C2, float* I, float* Td, float* KParameter, int NumofLight, int NumofInstruction, int NumofTd, int ParticleNoinGroup)
	{
		//SetGBest
		GBestFitness = NewGBestFitness;

		//SetGBest & UpdateVelocity & Position 
		for (int i = 0; i<NumofLight; i++)//i是維度
		{
			//SetGBest
			GBest[i] = *(GBestPointer + i);

			//UpdateVelocity
			Velocity[i] = (C0*Velocity[i] + C1*(float)device_RandomNumber(ParticleNoinGroup, 0, 100) / 100 * (PBest[i] - Position[i]) + C2*(float)device_RandomNumber(ParticleNoinGroup,0, 100) / 100 * (GBest[i] - Position[i]));
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
}//PSO類別結束


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
	/***********************讀取檔案**********************/
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));

	//Number of Particles
	fr.open("PSOParameter.txt", ios::in); //C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\PSOParameter.txt
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
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
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofInstruction
		sscanf(FileInput, "%d", &NumofInstruction);
		//Host_I = (float**)new2d(MaxNumofInstruction, 2, sizeof(float));
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_I[i][j]);        //字串轉數字
				//cout << Host_I[i][j] << ",";
			}
			//cout<<endl;
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\Td.txt
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofTd
		sscanf(FileInput, "%d", &NumofTd);
		Host_Td = new float[NumofTd];
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf(FileInput, "%f", &Host_Td[i]);        //字串轉數字
			//cout<<Host_Td[i]<<endl;
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\KParameter.txt
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣高
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣寬
		sscanf(FileInput, "%d", &NumofLight);
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf(FileInput, "%f", &Host_KParameter[i][j]);        //字串轉數字
				//cout<<Host_KParameter[i][j]<<",";
			}
			//cout<<endl;
		}
	}
	fr.close();
	//資料讀取完畢

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

	//GPU準備空間
	float* I;   //(耗能,亮度)
	float* Td;
	float* KParameter;

	size_t SizeofI = MaxNumofInstruction * 2 * sizeof(float);
	cudaMalloc(&I, SizeofI);

	size_t SizeofTd = NumofTd*sizeof(float);
	cudaMalloc(&Td, SizeofTd);

	size_t SizeofKParameter = MaxNumofTd*(MaxNumofLight + OtherLight)*sizeof(float);
	cudaMalloc(&KParameter, SizeofKParameter);

	//將資料移至GPU
	cudaMemcpy(I, Host_I, SizeofI, cudaMemcpyHostToDevice);
	cudaMemcpy(Td, Host_Td, SizeofTd, cudaMemcpyHostToDevice);
	cudaMemcpy(KParameter, Host_KParameter, SizeofKParameter, cudaMemcpyHostToDevice);
	//移動完畢

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
	StartTime = clock(); //開始時間


	/********************************隨機產生初始值*****************************************/
	int* RandomVelocity = new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ProduceParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)
		{
			int* j = new int[NumofLight];
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

	//初值產生完畢
	//同時，第一輪PSO結束，剩下CircleofPSO-1次
	//新增平行處理，在最外層
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

			//cout << "Group" << GroupNo + 1 << ", 現在的平均相似度:" << sim_ave << endl;
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

		//認為求出的解已穩定
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
			RepeatValue = TotalGBest[NumofLight];
		else if (RepeatValue - TotalGBest[NumofLight] < TotalGBest[NumofLight] * DifferentRange && RepeatValue != -1)
		{
			cout << "外迴圈次數: " << GroupCircle << endl;
			break;
		}

		//打散部分分群，重新組合
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
		//重分完畢//

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			for (int i = 0; i<NumofLight; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			GBestFitness[GroupNo] = TotalGBestFitness;
		}
	}
	EndTime = clock();//結束時間


	/***PSO執行結束，觀看結果***
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
	/*****************************/

	cudaFree(d_GBest);
	cudaFree(d_Group);
	cudaFree(d_GBestFitness);

	cudaFree(I);
	cudaFree(Td);
	cudaFree(KParameter);

	/*****write file*****/
	filename = "The_Best_Solution.txt";//C:\\ITLab_cll\\LightIllumination\\LightIlluminationProgram\\The_Best_Solution.txt
	fw.open(filename, ios::out);//開啟檔案
	if (!fw){//如果開啟檔案失敗，fw為0；成功，fw為非0
		cout << "Fail to open file: " << filename << endl;
	}
	//cout << "燈具指令:" << endl;
	if (TotalGBestFitness != ImpossibleResume)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			//cout << TotalGBest[i] << ", ";
			fw << TotalGBest[i] << endl;
		}
		//cout << endl;
		//cout << "耗能:" << TotalGBestFitness << endl;
		fw << TotalGBestFitness << endl;
	}
	else
	{
		//cout << "無解" << endl;
		fw << "無解" << endl;
	}
	//cout << "耗時: " << double(EndTime - StartTime) / CLOCKS_PER_SEC << endl;
	fw << double(EndTime - StartTime) / CLOCKS_PER_SEC << endl;
	fw.close();//關閉檔案

	//cout << "整體平均跳出內圈50時的相似度: " << total_ave_sim / (NumofGroupCircle*NumofGroup) << endl;

	//釋放空間             //GroupInfo 非動態陣列，不用自行刪除  //現在是動態了
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
給予希望的最大值與最小值，就能產生介於這兩個值中的隨機數，這個隨機亂數也可能等於最大值，或等於最小值。
其亂數值是依照時間所產生，所以重複呼叫時，會給予不同的亂數值(整數)。
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
動態產生2為array
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
計算相似度
/**********************************************/
float similarity(int* A, int* B, int Length)
{
	//for case 1
	float variation = 0, innerProduct = 0;
	float Similarity;

	//Euclidean distance(越小越像)
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

	//Euclidean distance(越小越像)
	for (int i = 0; i < Length; i++)
	{
		innerProduct += A[i] * A[i];
		variation += (A[i] - B[i]) * (A[i] - B[i]);
	}
	Similarity = sqrt(variation) / sqrt(innerProduct);

	//cout <<"variation:"<< fixed << setprecision(5) << Similarity << ", ";
	return Similarity;
}