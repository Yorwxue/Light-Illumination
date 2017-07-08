#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;



//讀檔取得
int NumofLight;       //燈的個數
int NumofInstruction; //亮度指令個數
int NumofTd;          //目標點個數

int RunTime = 1;  //跑多少次
float RealBest = -1;

#define OtherLight 0        //額外光源-太陽
//float OtherLightT[OtherLight]={0};        //額外光源的亮度

int NumofParticle;        //取NumofParticle組解作為PSO的初值(粒子數)
int NumofRedistribute;

float RepeatValue = -1;
float DifferentRange;
int NumofGroup;

//演算法定義
#define ImpossibleResume 99999
#define InitialW 0.9             //慣性係數，由大到小
#define FinalW 0.4
#define InitialCOne 2.5          //本身的C係數，由大到小，因為搜尋初期重視個體
#define FinalCOne 1
#define InitialCTwo  1       //群體的C係數，由小到大，因為搜尋末期重視群體
#define FinalCTwo 2.5
int CircleofPSO;
int VelocityLimit;
int NumofGroupCircle;
double ResumeTime;

void Parameter()
{
	CircleofPSO = 50;//NumofParticle*(0.02);      //PSO執行的循環次數，粒子數的2%

	VelocityLimit = NumofInstruction/3;//NumofInstruction*(0.5);  //移動速度上限

	NumofGroupCircle = 10;//NumofGroup;    //各處理器間交流所需，約為處理器個數

	NumofRedistribute = NumofParticle / 2;  //要交換的粒子數
}
//


int CheckT(float T, int MatrixRow, int i, int j, int k);

void PSO();

void* new2d(int h, int w, int size);
int RandomNumber(int MinValue, int MaxValue);


float** I;
float* Td;
float** KParameter;

//額外不可控制光源處理(陽光)
/*
void OtherLightProcess()
{
for(int i=0;i<NumofTd;i++)
{
Td[i]-=KParameter[i][NumofLight+OtherLight-1]*OtherLightT[OtherLight];
}
}
*/
class Particle  //每個粒子的基本設定，變動調光指令版本
{
private:
	float* PBest=new float[NumofLight + 1],  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,個體最低耗能)
		*GBest = new float[NumofLight + 1];  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,群體最低耗能)
	int* Position = new int[NumofLight];  //記錄個體，其各盞燈光的調光指令 <----> Position   //其中燈光數也是solution的維度，第幾盞燈表示第幾維
	int* Velocity = new int[NumofLight];             //調光指令上升或下降 <----> Velocity
	float W;//模擬慣性用的weighting，其值越大，越適合全域搜尋，建議由0.9~0.4
	float C1, C2;//係數，調整偏向各體，或偏向群體
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
		for (int i = 0; i<NumofLight; i++)//i是維度
		{
			Position[i] = Position[i] + Velocity[i];
			if (Position[i]>NumofInstruction - 1)             //調光指令的上下限
				Position[i] = NumofInstruction - 1;
			if (Position[i]<0)
				Position[i] = 0;
		}
	}
	void UpdateVelocity()
	{
		for (int i = 0; i<NumofLight; i++)//i是維度
		{
			Velocity[i] = (W*Velocity[i] + C1*((float)RandomNumber(0, 100) / 100)*(PBest[i] - Position[i]) + C2*((float)RandomNumber(0, 100) / 100)*(GBest[i] - Position[i]));
			if (Velocity[i]>VelocityLimit)                   //速度上下限
				Velocity[i] = VelocityLimit;
			if (Velocity[i]<(-VelocityLimit))
				Velocity[i] = (-VelocityLimit);
		}
	}
	void CheckPosition()
	{
		float Tpso = 0;  //根據目前粒子的位置所算出來的T(照度)值
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
//PSO類別結束



int main()
{
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));

	//Number of Particles
	fr.open("PSOParameter.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
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
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open LightInstruction file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofInstruction
		sscanf_s(FileInput, "%d", &NumofInstruction);

		I = (float **)new2d(NumofInstruction, 2, sizeof(float));

		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf_s(FileInput, "%f", &I[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open Td file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofTd
		sscanf_s(FileInput, "%d", &NumofTd);

		Td = new float[NumofTd];

		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			sscanf_s(FileInput, "%f", &Td[i]);        //字串轉數字
		}
	}
	fr.close();

	//KParameter
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open KParameter file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣高
		fr.getline(FileInput, sizeof(FileInput), ',');//矩陣寬
		sscanf_s(FileInput, "%d", &NumofLight);

		KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));

		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<(NumofLight + OtherLight); j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				sscanf_s(FileInput, "%f", &KParameter[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();
	//資料讀取完畢

	Parameter();

	//OtherLightProcess();//額外光源處理

	/*******PSO的相關結果**********/
	PSO();

	/***釋放空間***/
	delete[] I;
	delete[] Td;
	delete[] KParameter;

	//system("pause");
	return 0;
}

//確認亮度合格與否
int CheckT(int MatrixRow, int i[])     //MatrixRow為要檢查是否合格的row，i為紀錄指令的矩陣
{
	float T = 0;
	for (int k = 0; k<NumofLight; k++)
	{
		T += KParameter[MatrixRow][k] * I[i[k]][1];//[表示燈的指令,1表示亮度]
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
	int ParticleNo = 0;                  //目標是找到NumofParticle個解作為初始值，所以設一個旗標用以計數
	Particle*** Group=(Particle***)new2d(NumofGroup,NumofParticle,sizeof(Particle*));
	float T = 0;

	/*********全部重分配時使用************/
	Particle*** TempGroup=(Particle***)new2d(NumofGroup,NumofParticle,sizeof(Particle*));
	int* TempGroupNo=new int[NumofGroup];  //計算交換粒子用的暫存群是否已滿
	for (int i = 0; i<NumofGroup; i++)
	{
		TempGroupNo[i] = 0;
	}
	/**************************************/
	PSOstart = clock();
	/********************************隨機產生初始值，不論其是否合格*****************************************/
	int* RandomVelocity=new int[NumofLight];
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
	{
		ParticleNo = 0;
		for (int i = 0; i<NumofParticle; i++)                                       //關於j的說明請搜尋"高維轉單維"，j存了所有燈的指令資訊
		{
			int* j=new int[NumofLight];
			//粒子初始位置隨機產生方式3
			for (int h = 0; h<NumofLight; h++)
			{
				j[h] = RandomNumber(0, NumofInstruction - 1);
			}

			//
			//初始位置產生完畢

			for (int MatrixRow = 0; MatrixRow<NumofTd; MatrixRow++)
			{
				T = CheckT(MatrixRow, j);
				if (T == 0)                                                     //不合格，Power設ImpossibleResume，不可能發生的情況，作為初值
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
				if (MatrixRow == NumofTd - 1)//全部Row都符合條件，開始計算耗能
				{
					for (int k = 0; k<NumofLight; k++)                                                 //計算耗能
					{
						Power += I[j[k]][0];
					}
					for (int h = 0; h<NumofLight; h++)                                                 //隨機初始速度
					{
						RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
					}
					Group[GroupNo][ParticleNo++] = new Particle(j, Power, &RandomVelocity[0]);                //產生粒子
					Power = 0;
				}
			}
		}
		for (int i = 0; i<NumofParticle; i++)
		{
			PBestPoint = Group[GroupNo][i]->GetPBest();
			if (*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume)    // GBest[NumofLight]是耗能
			{
				for (int j = 0; j<NumofLight + 1; j++)
				{
					GBest[GroupNo][j] = *(PBestPoint + j);
				}
			}
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
					if ((*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume) && *(PBestPoint + NumofLight) != ImpossibleResume)    // GBest[NumofLight]是耗能
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


		//重分方式
		/****************************************************
		//比例打包重分
		int RedistributeList[NumofGroup];
		int RedistributeParticleList[NumofGroup][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;

		//建立群的list
		for(int i=0;i<NumofGroup;i++)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		for(int j=0;j<i;j++)
		{
		if(RedistributeList[j]==RandomGroup)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		j=-1;                                                             //for loop關係，檢查前會加1，所以設-1，讓他從0開始檢查
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
		j=-1;                                                //for loop關係，檢查前會加1，所以設-1，讓他從0開始檢查
		}
		}
		RedistributeParticleList[k][i]=RandomParticle;
		}
		}

		//開始交換
		for(int i=0;i<NumofGroup;i++)  //備份
		{
		for(int j=0;j<NumofParticle;j++)
		{
		TempGroup[i][j]=Group[i][j];
		}
		}
		for(int i=0;i<NumofGroup;i++) //i是所在的群,k是第幾次作群之間的交換
		{
		for(int j=0;j<NumofRedistribute;j++)
		{
		Group[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempGroup[i][ RedistributeParticleList[i][j] ];
		}
		}
		/************全部重分************
		int RandomGroup;
		for(int i=0;i<NumofGroup;i++)
		{
		for(int j=0;j<NumofParticle;j++)
		{
		RandomGroup=RandomNumber(0,NumofGroup-1);
		while(TempGroupNo[RandomGroup]==NumofParticle)//暫存群已滿
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
		//////////////打散部分分群，重新組合
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
		///////////*///重分完畢

		/*尋找重分後的各群的gbest*/
		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
		    for(int i=0;i<NumofParticle;i++)
		    {
				GBest[GroupNo][NumofLight] = ImpossibleResume;
		    }
		    for(int i=0;i<NumofParticle;i++)
		    {
				PBestPoint = Group[GroupNo][i]->GetPBest();
				if ((*(PBestPoint + NumofLight)< GBest[GroupNo][NumofLight] || GBest[GroupNo][NumofLight] == ImpossibleResume) && *(PBestPoint + NumofLight) != ImpossibleResume)    // GBest[NumofLight]是耗能
		        {
		            for(int j=0;j<NumofLight+1;j++)
		            {
						GBest[GroupNo][j] = *(PBestPoint + j);
		            }
		        }
		    }
		}

		/****尋找完畢****/

		for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)
		{
			/*各群重新找出的gbest寫回*/
			for(int i=0;i<NumofParticle;i++)
			{
				Group[GroupNo][i]->SetGBest(&GBest[GroupNo][0]);
			}
			/******/
			/*將整合的Gbest寫回給各群*
			for (int i = 0; i<NumofLight + 1; i++)
			{
				GBest[GroupNo][i] = TotalGBest[i];
			}
			/*********/
		}

		//認為求出的解已穩定
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
			RepeatValue = TotalGBest[NumofLight];
		else if (RepeatValue - TotalGBest[NumofLight]<TotalGBest[NumofLight] * DifferentRange && RepeatValue != -1)
			break;
	}

	PSOstop = clock(); //結束時間

	/***PSO執行結束，觀看結果***
	if (TotalGBest[NumofLight] != ImpossibleResume)
	{
		//cout<<"PSO搜尋"<<endl
		cout << "指令:";
		for (int i = 0; i<NumofLight; i++)
		{
			cout << TotalGBest[i] << ",";
		}
		cout << "耗能:" << TotalGBest[NumofLight] << endl;
	}
	else
	{
		cout << "無解" << endl;
	}

	/****寫檔****/
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
			fw << "無解" << endl;
		}
		fw << double(PSOstop - PSOstart) / CLOCKS_PER_SEC << endl;
		fw.close();
	}

	/*各點照度*/
	fstream fw;
	string filename = "RealIllumination.txt";
	fw.open(filename, ios::out);
	float RealT = 0;
	for (int MatrixRow = 0; MatrixRow < NumofTd; MatrixRow++)
	{
		for (int k = 0; k < NumofLight; k++)
		{
			RealT += KParameter[MatrixRow][k] * I[(int)TotalGBest[k]][1];//[表示燈的指令,1表示亮度]
		}
		fw << RealT << endl;
		RealT = 0;
	}
	fw.close();
	/****寫檔完畢****/

	cout << "費時:" << double(PSOstop - PSOstart) / CLOCKS_PER_SEC << "秒" << endl << endl;

	/*****************************/
	for (int GroupNo = 0; GroupNo<NumofGroup; GroupNo++)  //釋放空間
	{
		for (int i = 0; i<NumofParticle; i++)
		{
			delete Group[GroupNo][i];
		}
	} /**/

	//return TotalGBest[NumofLight];
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