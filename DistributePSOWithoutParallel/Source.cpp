#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//使用者定義
#define NumofLight 12       //燈的個數
#define NumofInstruction 16 //亮度指令個數
#define NumofTd 5          //目標點個數
#define NumofCPU 10

int RunTime=10;  //跑多少次
float RealBest=246.34;
float totaltime=0;
//float LongestTime=0;
//float ShortestTime=999999;

#define OtherLight 0        //額外光源-太陽
//float OtherLightT[OtherLight]={0};        //額外光源的亮度


#define NumofParticle 1000        //取NumofParticle組解作為PSO的初值(粒子數)
#define NumofRedistribute 100     //10%

//#define EndFlagRepeatMAXTimes 2
//int EndFlagRepeatTimes=0;
float RepeatValue=-1;
float DifferentRange=0.0005;

float ave_power = 0;

//程式自定義
#define ImpossibleResume 99999
#define InitialW 0.9             //慣性係數，由大到小
#define FinalW 0.4
#define InitialCOne 2.5          //本身的C係數，由大到小，因為搜尋初期重視個體
#define FinalCOne 1
#define InitialCTwo  1       //群體的C係數，由小到大，因為搜尋末期重視群體
#define FinalCTwo 2.5
int CircleofPSO;
int VelocityLimit;
int NumofCPUCircle;
double ResumeTime;

void Parameter()
{
    CircleofPSO = 50;//NumofParticle*(0.02);      //PSO執行的循環次數，粒子數的2%
	//if(CircleofPSO<10)   //探索次數過少會造成無解
	//	CircleofPSO=10;

    VelocityLimit=5;//NumofInstruction*(0.5);  //移動速度上限，移動空間的50%
	//if(VelocityLimit<3)   //移動限制太嚴，會導致個體停止不動
	//	VelocityLimit=3;  

    NumofCPUCircle=10;//NumofCPU;    //各處理器間交流所需，約為處理器個數
}
//



float I[NumofInstruction][2];   //(耗能,亮度)
float Td[NumofTd];
float KParameter[NumofTd][NumofLight+OtherLight];

float start_point_range[NumofLight][2];    //(低照度指令,高照度指令)

int CheckT(float T,int MatrixRow,int i,int j,int k);
void TestAllSollution();
float PSO();
int RandomNumber(int MinValue,int MaxValue);
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
	float PBest[NumofLight+1],  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,個體最低耗能)
		  GBest[NumofLight+1];  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,群體最低耗能)
	int Position[NumofLight];  //記錄個體，其各盞燈光的調光指令 <----> Position   //其中燈光數也是solution的維度，第幾盞燈表示第幾維
	int Velocity[NumofLight];             //調光指令上升或下降 <----> Velocity
	float W;//模擬慣性用的weighting，其值越大，越適合全域搜尋，建議由0.9~0.4
	float C1,C2;//係數，調整偏向各體，或偏向群體
public:
	Particle(float* PPointer,float Power,int* VPointer);
	Particle(int RangeofLightwithInstruction[],float Power,int* VPointer);
	void SetW(float NewW,float NewC1,float NewC2)
	{
		W=NewW;
		C1=NewC1;
		C2=NewC2;
	}
	float* GetPBest(){return &PBest[0];}
	void SetGBest(float* GBestPointer)
	{
		for(int i=0;i<(NumofLight+1);i++)
		{
			GBest[i]=*(GBestPointer+i);                                                                 //cout<<"Gbest["<<i<<"]:"<<GBest[i]<<",";
		}                                                                                               //system("pause");              
	}
	void UpdatePosition()
	{
		for(int i=0;i<NumofLight;i++)//i是維度
		{                                                                                               //cout<<"原P"<<i<<":"<<Position[i]<<",";
		Position[i]=Position[i]+Velocity[i];                                                            //cout<<"V:"<<Velocity[i]<<endl;
			if(Position[i]>NumofInstruction-1)             //調光指令的上下限
				Position[i]=NumofInstruction-1;
			if(Position[i]<0)
				Position[i]=0;                                                                           //cout<<"後P"<<i<<":"<<Position[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void UpdateVelocity()
	{
		for(int i=0;i<NumofLight;i++)//i是維度
		{                                                                                               //cout<<"原V"<<i<<":"<<Velocity[i]<<",";
			Velocity[i] = (W*Velocity[i] + C1*((float)RandomNumber(0, 100) / 100)*(PBest[i] - Position[i]) + C2*((float)RandomNumber(0, 100) / 100)*(GBest[i] - Position[i]));         //cout<<"後V"<<i<<":"<<Velocity[i]<<endl;
			if(Velocity[i]>VelocityLimit)                   //速度上下限
				Velocity[i]=VelocityLimit;
			if(Velocity[i]<(-VelocityLimit))
				Velocity[i]=(-VelocityLimit);                                                            //cout<<"後V"<<i<<":"<<Velocity[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void CheckPosition()
	{
		float Tpso=0;  //根據目前粒子的位置所算出來的T(照度)值
		for(int i=0;i<NumofTd;i++)
		{
			for(int j=0;j<NumofLight;j++)
			{
				    Tpso+=KParameter[i][j]*I[Position[j]][1];
			}
			if(Tpso<Td[i])
			{
				break;
			}
			Tpso=0;
			if(i==(NumofTd-1))
			{
				ComputingFitness();
			}
		}
	}
	void ComputingFitness()
	{
		float Fitness=0;
		for(int i=0;i<NumofLight;i++)
		{                                                                                            //cout<<i<<":"<<Position[i]<<",";
		    Fitness+=I[Position[i]][0];
		}                                                                                            //cout<<endl;system("pause");
		if(Fitness<PBest[NumofLight]||PBest[NumofLight]<=0)
		{                                                                                            //cout<<Fitness<<"<"<<PBest[NumofLight]<<",";
			for(int j=0;j<NumofLight;j++)
			{
				PBest[j]=Position[j];                                                               //cout<<PBest[j];
			}
			PBest[NumofLight]=Fitness;                                                              //cout<<"New PBest:"<<PBest[NumofLight]<<endl;
		}
	}
	
	
};

Particle::Particle(float* PPointer,float Power,int* VPointer)
{
	for(int i=0;i<NumofLight;i++)
	{
		Position[i]=*(PPointer+i);                                                          
		Velocity[i]=*(VPointer+i);                                                          
		PBest[i]=Position[i];
	}
	PBest[NumofLight]=Power;
	W=InitialW;
	C1=InitialCOne;
	C2=InitialCTwo;
}
Particle::Particle(int RangeofLightwithInstruction[],float Power,int* VPointer)
{
	for(int i=0;i<NumofLight;i++)
	{
		Velocity[i]=*(VPointer+i);                                                                                    //cout<<"V"<<i<<":"<<Velocity[i]<<" ,";
	}                                                                                                                 //cout<<endl;
	for(int k=0;k<NumofLight;k++)
	{
		Position[k]=RangeofLightwithInstruction[k];                                                                                                //cout<<"P"<<k<<":"<<Position[k]<<" ,";
		PBest[k]=RangeofLightwithInstruction[k];
	}                                                                                                                 //cout<<endl;system("pause");
	PBest[NumofLight]=Power;
	W=InitialW;
	C1=InitialCOne;
	C2=InitialCTwo;
}
//PSO類別結束



int main()
{
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt",ios::in);
	if(!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//這是NumofInstruction
		for(int i=0;i<NumofInstruction;i++)
		{
			for(int j=0;j<2;j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&I[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt",ios::in);
	if(!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//這是NumofTd
		for(int i=0;i<NumofTd;i++)
		{
			fr.getline(FileInput,sizeof(FileInput),',');
			sscanf(FileInput,"%f",&Td[i]);        //字串轉數字
		}
	}
	fr.close();
	
	//KParameter
	fr.open("KParameter.txt",ios::in);
	if(!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//矩陣高
		fr.getline(FileInput,sizeof(FileInput),',');//矩陣寬
		for(int i=0;i<NumofTd;i++)
		{
			for(int j=0;j<(NumofLight+OtherLight);j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&KParameter[i][j]);        //字串轉數字
			}
		}
	}
	fr.close();

	//start range
	fr.open("start_point.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		for (int i = 0; i<NumofLight; i++)
		{
				fr.getline(FileInput, sizeof(FileInput), ',');
				float prob_solu;
				sscanf(FileInput, "%f", &prob_solu); //字串轉數字
				for (int j = 0; j <NumofInstruction; j++)
				{
					if (I[j][1] == prob_solu)
					{
						start_point_range[i][0] = j;
						start_point_range[i][1] = j;
						break;
					}
					else if (I[j][1] < prob_solu)
					{
						start_point_range[i][0] = j;
						start_point_range[i][1] = j - 1;
						break;
					}
				}
		}
	}
	fr.close();

	//資料讀取完畢

	Parameter();

	//OtherLightProcess();//額外光源處理
	/*******PSO的相關結果**********/
	clock_t start, stop;
    start = clock(); //開始時間
	float AverageResult=0;
	float TestAnswer;//是否無解或非最佳解
	int BestTimes=0;
	int Range01Times=0;
	int Range02Times=0;
	int Range03Times=0;
	int Range04Times=0;
	int Range05Times=0;
	int Range1Times=0;
	int Range2Times=0;
	int Range3Times=0;
	int Range4Times=0;
	int Range5Times=0;
	int RangeOther=0;
	int NoAnswer=0;//無解總次數
	
	for(int i=0;i<RunTime;i++)
	{
		TestAnswer=PSO();
		ave_power += TestAnswer;
		if(TestAnswer<RealBest || RealBest==-1)
		{
			RealBest=TestAnswer;
		}
		if(TestAnswer==0)
		{
			//cout<<TestAnswer<<endl;
			NoAnswer++;
		}
		else if(TestAnswer<=(RealBest+0.001))
		{
		    //cout<<TestAnswer<<endl;
			BestTimes++;
		}/*
		else if( (TestAnswer-RealBest)<=(RealBest*0.001) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range01Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.002) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range02Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.003) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range03Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.004) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range04Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.005) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range05Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.01) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range1Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.02) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range2Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.03) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range3Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.04) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range4Times++;
		}
		else if( (TestAnswer-RealBest)<=(RealBest*0.05) && TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			Range5Times++;
		}
		else if(TestAnswer!=ImpossibleResume)
		{
		    //cout<<TestAnswer<<endl;
			RangeOther++;
		}
		else
		{
			NoAnswer++;
		}/**/
	    //AverageResult+=TestAnswer;
	}
	stop = clock(); //結束時間
	
	cout<<endl;
	cout<<"執行"<<RunTime<<"次"<<endl;
	cout<<"最佳解為："<<RealBest<<endl;
	//cout<<"平均耗能:"<<(AverageResult/(RunTime-NoAnswer))<<endl;
	cout << "平均耗能" << ave_power / RunTime << endl;
	cout << "平均費時:" << totaltime / RunTime << "秒" << endl;/**
		<<"最長費時:"<<LongestTime<<"秒"<<endl<<endl
		<<"最短費時:"<<ShortestTime<<"秒"<<endl<<endl;/**/
	//cout<<"分散費時:"<<ResumeTime<<"秒"<<endl<<endl;
	//cout<<"無解次數:"<<NoAnswer<<"次"<<endl;
	cout << "最佳解個數:" << BestTimes << endl/*
	    <<"0.1%內個數:"<<Range01Times<<endl
		<<"0.2%內個數:"<<Range02Times<<endl
		<<"0.3%內個數:"<<Range03Times<<endl
		<<"0.4%內個數:"<<Range04Times<<endl
		<<"0.5%內個數:"<<Range05Times<<endl
		<<"1%內個數:"<<Range1Times<<endl
		<<"2%內個數:"<<Range2Times<<endl
		<<"3%內個數:"<<Range3Times<<endl
		<<"4%內個數:"<<Range4Times<<endl
		<<"5%內個數:"<<Range5Times<<endl
		<<"Other:"<<RangeOther<<endl/**/
		<<"無解個數："<<NoAnswer<<endl;
	//cout<<"10%內個數"<<Range10Times<<endl;
	/***********************************************/
	system("pause");
	return 0;
}

//確認亮度合格與否
int CheckT(int MatrixRow,int i[])     //MatrixRow為要檢查是否合格的row，i為紀錄指令的矩陣
{
	float T=0;
	for(int k=0;k<NumofLight;k++)
	{
	        T+=KParameter[MatrixRow][k]*I[i[k]][1];//[表示燈的指令,1表示亮度]
	}
	if(T<Td[MatrixRow])
		return 0;
	else 
		return 1;
}
/**********************************************************************************************/
float PSO()
{
	clock_t PSOstart, PSOstop;
	clock_t Distribute_start, Distribute_stop;
	float TotalGBest[NumofLight+1];
	TotalGBest[NumofLight]=ImpossibleResume;
																		
	float GBest[NumofCPU][NumofLight+1],* PBestPoint,Power=0;
	for(int i=0;i<NumofCPU;i++)
	    GBest[i][NumofLight]=ImpossibleResume;
	int ParticleNo=0;                  //目標是找到NumofParticle個解作為初始值，所以設一個旗標用以計數
	Particle* Group[NumofCPU][NumofParticle];
	float T=0;

	/*********全部重分配時使用************/
	Particle* TempGroup[NumofCPU][NumofParticle];
	int TempGroupNo[NumofCPU];  //計算交換粒子用的暫存群是否已滿
	for(int i=0;i<NumofCPU;i++)
	{
		TempGroupNo[i]=0;
	}
	/**************************************/

    //Distribute_start = clock(); //開始時間

    float StartPointRatio[NumofLight];
	int PSOStartCenter[NumofLight];

	PSOstart = clock();
/********************************隨機產生初始值，不論其是否合格*****************************************/
	int RandomVelocity[NumofLight];
	for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
	{
	    ParticleNo=0;
		for(int i=0;i<NumofParticle;i++)                                       //關於j的說明請搜尋"高維轉單維"，j存了所有燈的指令資訊
	    {
			int j[NumofLight];
			//粒子初始位置隨機產生方式
			/**/
		    for(int h=0;h<NumofLight;h++)                                       
	        {                                                                                                                                                    //cout<<StartPointCheck<<",";
		    	j[h]=RandomNumber(0,NumofInstruction-1);   
	    	}/**/                                                                                          //cout<<endl;
			//特定位置週邊產生
			/**
			for (int h = 0; h<NumofLight; h++)
			{                                                                                                                                                    //cout<<StartPointCheck<<",";
				j[h] = start_point_range[h][RandomNumber(0, 1)] +(RandomNumber(0, 1) ? 1 : -1) * RandomNumber(0, (int)(NumofInstruction / 5));
				if (j[h] >= NumofInstruction)
				{
					j[h] = NumofInstruction - 1;
				}
				else if (j[h] < 0)
				{
					j[h] = 0;
				}
			}/**/
			//
	    	//初始位置產生完畢

	    	for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	{
		        T=CheckT(MatrixRow,j);
	    	    if(T==0)                                                     //不合格，Power設ImpossibleResume，不可能發生的情況，作為初值
		        {
		    		for(int h=0;h<NumofLight;h++)
	                {
		                RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                }
			        Group[CPUNo][ParticleNo++]=new Particle(j,ImpossibleResume,&RandomVelocity[0]);
			        break;
		        }
		    	else
		    	{
			    	T=0;
		    	}
		        if(MatrixRow==NumofTd-1)//全部Row都符合條件，開始計算耗能
		        {                                                                             //cout<<"前"<<Power<<endl;system("pause");
		    	    for(int k=0;k<NumofLight;k++)                                                 //計算耗能
		    	    { 
			    	    Power+=I[j[k]][0];                                                       //cout<<Power<<endl;system("pause");
			        }                                                                         //cout<<"後"<<Power<<endl;system("pause");
		    		for(int h=0;h<NumofLight;h++)                                                 //隨機初始速度
	                {
		                RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                }                                                                                     //cout<<Power<<endl;system("pause");
                    Group[CPUNo][ParticleNo++]=new Particle(j,Power,&RandomVelocity[0]);                //產生粒子
		    		Power=0;
		        }
		    }
	    }
	    for(int i=0;i<NumofParticle;i++)
	    {
	        PBestPoint=Group[CPUNo][i]->GetPBest();                                //cout<<*(PBestPoint+NumofLight)<<endl;system("pause");
	    	if(*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)    // GBest[NumofLight]是耗能
	    	{
		    	for(int j=0;j<NumofLight+1;j++)
	            {
		    		GBest[CPUNo][j]=*(PBestPoint+j);                               //cout<<GBest[CPUNo][j]<<",";system("pause");
	            }                                                           //cout<<GBest[CPUNo][NumofLight]<<",";system("pause");//cout<<endl;system("pause");
		    }
	    }
	}
	                                                                                             //cout<<"第1輪不合格粒子數:"
	//Distribute_stop = clock();//結束時間
	//ResumeTime=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
	//初值產生完畢
	//同時，第一輪PSO結束，剩下CircleofPSO-1次
	//新增平行處理，在最外層
/*******************************************************************************************************/
	for(int CPUCircle=0;CPUCircle<NumofCPUCircle;CPUCircle++)
	{
	    for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
	    {
			//if(CPUNo==0)
			//{
			//	 Distribute_start = clock(); //開始時間
			//}
	        for(int i=1;i<CircleofPSO;i++)
	        {
		        float NewW=InitialW-(InitialW-FinalW)*((i+CPUCircle)/(NumofCPUCircle+CircleofPSO-1));
		        float NewC1=InitialCOne-(InitialCOne-FinalCOne)*((i+CPUCircle)/(NumofCPUCircle+CircleofPSO-1));
		        float NewC2=InitialCTwo+(FinalCTwo-InitialCTwo)*((i+CPUCircle)/(NumofCPUCircle+CircleofPSO-1));
		        for(int i=0;i<NumofParticle;i++)
	            {
			        Group[CPUNo][i]->SetW(NewW,NewC1,NewC2);
		            Group[CPUNo][i]->SetGBest(&GBest[CPUNo][0]);
			        Group[CPUNo][i]->UpdateVelocity();
			        Group[CPUNo][i]->UpdatePosition();
			        Group[CPUNo][i]->CheckPosition();  
				}
				for(int i=0;i<NumofParticle;i++)
				{
			        PBestPoint=Group[CPUNo][i]->GetPBest();
			        if((*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)&&*(PBestPoint+NumofLight)!=ImpossibleResume)    // GBest[NumofLight]是耗能
			        {
			            for(int j=0;j<NumofLight+1;j++)
	                    {
	        	            GBest[CPUNo][j]=*(PBestPoint+j);
				       }
			        }
	            }                                                                                          
				if(GBest[CPUNo][NumofLight]<TotalGBest[NumofLight] || TotalGBest[NumofLight]==ImpossibleResume)
				{
				    for(int p=0;p<NumofLight+1;p++)
					    TotalGBest[p]=GBest[CPUNo][p];
			    }
	        }
			//if(CPUNo==0)
			//{
			//	Distribute_stop = clock();//結束時間
		    //    ResumeTime+=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
			//}
	    }

		//Distribute_start = clock(); //開始時間
		
		

		//重分方式
		/****************************************************
		//比例打包重分
		int RedistributeList[NumofCPU];
		int RedistributeParticleList[NumofCPU][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;

		//建立群的list
		for(int i=0;i<NumofCPU;i++)  
		{
			RandomGroup=RandomNumber(0,NumofCPU-1);
			for(int j=0;j<i;j++)
			{
			    if(RedistributeList[j]==RandomGroup)
			    {
			        RandomGroup=RandomNumber(0,NumofCPU-1);
				    j=-1;                                                             //for loop關係，檢查前會加1，所以設-1，讓他從0開始檢查
				}
			}
			RedistributeList[i]=RandomGroup;
		}

		//選擇要交換的粒子，建立對照表
		for(int k=0;k<NumofCPU;k++)
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
		for(int i=0;i<NumofCPU;i++)  //備份
		{
			for(int j=0;j<NumofParticle;j++)
			{
				TempGroup[i][j]=Group[i][j];
			}
		}
		for(int i=0;i<NumofCPU;i++) //i是所在的群,k是第幾次作群之間的交換
		{
			for(int j=0;j<NumofRedistribute;j++)
		    {
				Group[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempGroup[i][ RedistributeParticleList[i][j] ];
			}
		}
		/************全部重分************
		int RandomGroup;
		for(int i=0;i<NumofCPU;i++)
		{
			for(int j=0;j<NumofParticle;j++)
			{
				RandomGroup=RandomNumber(0,NumofCPU-1);
			    while(TempGroupNo[RandomGroup]==NumofParticle)//暫存群已滿
			    {
				    RandomGroup=RandomNumber(0,NumofCPU-1);
			    }
				TempGroup[RandomGroup][TempGroupNo[RandomGroup]++] = Group[i][j];
			}
		}
		for(int i=0;i<NumofCPU;i++)
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
		int RandomCPU;
		int RandomParticle1,RandomParticle2;
		Particle *TempMemory;
		for(int j=0;j<NumofCPU;j++)
		{
		    for(int i=0;i<NumofRedistribute;i++)
		    {
			    RandomCPU=RandomNumber(0,NumofCPU-1);
		    	while(RandomCPU==j)
			    {
		    		RandomCPU=RandomNumber(0,NumofCPU-1);
		    	}
		    	RandomParticle1=RandomNumber(0,NumofParticle-1);
		    	RandomParticle2=RandomNumber(0,NumofParticle-1);
			    //
			    TempMemory=Group[j][RandomParticle1];
			    Group[j][RandomParticle1]=Group[RandomCPU][RandomParticle2];
			    Group[RandomCPU][RandomParticle2]=TempMemory;
		    }
		}
		///////////*///重分完畢

		//產生新粒子，而不重分
		/****************************************************
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)  
		{
			//建表
			int RedistributeParticleList[NumofRedistribute];
			for(int i=0;i<NumofRedistribute;i++)
		    {
			    int RandomReplaced=RandomNumber(0,NumofParticle-1);
			    for(int j=0;j<i;j++)
			    {
			        if(RedistributeParticleList[j]==RandomReplaced)
			        {
			            RandomReplaced=RandomNumber(0,NumofParticle-1);
				        j=-1;                                                //for loop關係，檢查前會加1，所以設-1，讓他從0開始檢查
				    }
			    }
			    RedistributeParticleList[i]=RandomReplaced;
			}
			//簽單完成，開始重新產生
		    for(int i=0;i<NumofRedistribute;i++)                                       //關於j的說明請搜尋"高維轉單維"，j存了所有燈的指令資訊
	        {
		    	int j[NumofLight];
		    	//粒子初始位置隨機產生方式
		        for(int h=0;h<NumofLight;h++)                                       
	            {                                                                                                                                                    //cout<<StartPointCheck<<",";
		        	j[h]=RandomNumber(0,NumofInstruction-1);   
	        	}                                                                                     
			
			    //
	        	//初始位置產生完畢

	        	for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	    {
		            T=CheckT(MatrixRow,j);
	    	        if(T==0)                                                     //不合格，Power設ImpossibleResume，不可能發生的情況，作為初值
		            {
		        		for(int h=0;h<NumofLight;h++)
	                    {
		                    RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                    }
						delete Group[CPUNo][RedistributeParticleList[i]];
			            Group[CPUNo][RedistributeParticleList[i]]=new Particle(j,ImpossibleResume,&RandomVelocity[0]);
			            break;
		            }
		        	else
		        	{
			        	T=0;
		        	}
		            if(MatrixRow==NumofTd-1)//全部Row都符合條件，開始計算耗能
		            {                                                                        
		    	        for(int k=0;k<NumofLight;k++)                                                 //計算耗能
		    	        { 
			    	        Power+=I[j[k]][0];                                                  
			            }  

						if(Power<TotalGBest[NumofLight])
						{
							TotalGBest[NumofLight]=Power;
						}

		    		    for(int h=0;h<NumofLight;h++)                                                 //隨機初始速度
	                    {
		                    RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                    }                                                                               
						delete Group[CPUNo][RedistributeParticleList[i]];
                        Group[CPUNo][RedistributeParticleList[i]]=new Particle(j,Power,&RandomVelocity[0]);                //產生粒子
		    	    	Power=0;
		            }
				}
			}
	    }
		/*************************產生完畢***************************/

		/*尋找重分或產生後的各群的gbest*
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
		{
			for(int i=0;i<NumofParticle;i++)
			{
				GBest[CPUNo][NumofLight]=ImpossibleResume;       
			}
			for(int i=0;i<NumofParticle;i++)
			{
			        PBestPoint=Group[CPUNo][i]->GetPBest();
			        if((*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)&&*(PBestPoint+NumofLight)!=ImpossibleResume)    // GBest[NumofLight]是耗能
			        {
			            for(int j=0;j<NumofLight+1;j++)
	                    {
	        	            GBest[CPUNo][j]=*(PBestPoint+j);
				       }
			        }
			} 
		}
		/****尋找完畢****/
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
		{
			/*各群重新找出的gbest寫回*/
			for(int i=0;i<NumofParticle;i++)
			{
			    Group[CPUNo][i]->SetGBest(&GBest[CPUNo][0]);
			}
			/******/
			/*將整合的Gbest寫回給各群*
			for(int i=0;i<NumofLight+1;i++)
			{
			    GBest[CPUNo][i]=TotalGBest[i];
			}
			/*********/
		}
		//Distribute_stop = clock();//結束時間
	    //ResumeTime+=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
		/********/
		/***認為求出的解已穩定***
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
		{
			RepeatValue=TotalGBest[NumofLight];
		}
		else if(RepeatValue!=-1)
		{
			 //理論上TotalGBest[NumofLight]一定小於或等於RepeatValue，因為RepeatValue就是上一輪的TotalGBest[NumofLight]，這個判斷式是保險作用
			if(RepeatValue>TotalGBest[NumofLight])
			{
				if (RepeatValue - TotalGBest[NumofLight] <= RepeatValue*DifferentRange)
				{
					break;
				}
			}
			else
			{
				if(TotalGBest[NumofLight]-RepeatValue <= TotalGBest[NumofLight]*DifferentRange)
				{
					break;
				}
			}
		}

		/*********/
	}
	
	/***PSO執行結束，觀看結果***/
	if(TotalGBest[NumofLight]!=ImpossibleResume)
	{
	    //cout<<"PSO搜尋"<<endl
		cout<<"指令:";                                                                            
	    for(int i=0;i<NumofLight;i++)
	    {
		    cout<<TotalGBest[i]<<",";  
	    }
	    cout<<"耗能:"<<TotalGBest[NumofLight];                                           
	}
	else
	{
		cout<<"無解"<<endl;
	}

	/****寫檔****/
	if(TotalGBest[NumofLight]<=RealBest || RealBest==-1)
	{
	    fstream fw;
	    string filename="The_Best_Solution.txt";
	    fw.open(filename, ios::out);
	    fw<<"PSO搜尋"<<endl<<"指令:";
	    for(int i=0;i<NumofLight;i++)
	    {
		    fw<<TotalGBest[i]<<",";  
	    }
	    fw<<endl<<"耗能:"<<TotalGBest[NumofLight]<<endl;    
		fw.close();
	}
	/****寫檔完畢****/

	PSOstop = clock(); //結束時間
	float ExecutingTime=float(PSOstop-PSOstart)/CLOCKS_PER_SEC;
    cout<<"費時:"<<ExecutingTime<<"秒"<<endl;
	totaltime+=ExecutingTime; 
	/**
	if(ExecutingTime>LongestTime)
		LongestTime=ExecutingTime;
	if(ExecutingTime<ShortestTime)
		ShortestTime=ExecutingTime;/**/
    /*****************************/
	for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)  //釋放空間
	 {
		 for(int i=0;i<NumofParticle;i++)
		 {
			 delete Group[CPUNo][i];
		 }
	 } /**/

	return TotalGBest[NumofLight];
}

/****************************************Random Number**********************************************
給予希望的最大值與最小值，就能產生介於這兩個值中的隨機數，這個隨機亂數也可能等於最大值，或等於最小值。
其亂數值是依照時間所產生，所以重複呼叫時，會給予不同的亂數值(整數)。
***************************************************************************************************/
int RandomNumber(int MinValue,int MaxValue)   //Both MinValue and MaxValue are included
{
	int R=(rand()%(MaxValue-MinValue+1))+MinValue;
	return R;
}
