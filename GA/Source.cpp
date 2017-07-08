#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
#include<iomanip>
#include<sstream>
using namespace std;

#define NumofLight 10        //燈的個數
#define NumofInstruction 16 //亮度指令個數
#define NumofTd 5          //目標點個數
#define Impossible 99999
#define MutationPossibility 10000
int RunTime = 50;
float RealBest = 224.22;

float ave_power = 0;

#define NumofGenome 100000
#define NumofCPU 1
#define NumofExchangeGenome 10000  //交流的基因組取10%


int NumofCrossoverGenome;  //排名較前面，用於交配的基因組,假設50%
int CircleofGA;  //基本次數，隨著點上升，需要稍微增加，但影響不大
int NumofCPUCircle;

float I[NumofInstruction][2];   //(耗能,亮度)
float Td[NumofTd];
float KParameter[NumofTd][NumofLight];

int CheckT(int MatrixRow, int i, int j, int k);
float GA();
int RandomNumber(int MinValue, int MaxValue);


class GeneAssemble
{
private:
	int Genome[NumofLight];
	float GenomeFitness;
public:
	GeneAssemble();
	GeneAssemble(int* GenomePointer, float Fitness);
	int* GetGenomePointer(){ return &Genome[0]; }
	float GetGenomeFitness(){ return GenomeFitness; }
	void SetGenome(int* GenomePointer)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			Genome[i] = *(GenomePointer + i);
		}
	}
	void mutation(int Which, int What)                        //Which:哪一個基因突變，What:變成什麼
	{
		Genome[Which] = What;
	}
	//////////////////////////////////
	void CheckGenome()
	{
		float Tga = 0;  //根據目前基因組所算出來的T(照度)值
		for (int i = 0; i<NumofTd; i++)
		{
			for (int j = 0; j<NumofLight; j++)
			{
				Tga += KParameter[i][j] * I[Genome[j]][1];
			}
			if (Tga<Td[i])
			{
				break;
			}
			Tga = 0;
			if (i == (NumofTd - 1))
			{
				GAComputingFitness();
			}
		}
	}
	void GAComputingFitness()
	{
		float Fitness = 0;
		for (int i = 0; i<NumofLight; i++)
		{
			Fitness += I[Genome[i]][0];
		}                                                                                                                                                                                     //cout<<Fitness<<"<"<<PBest[NumofLight]<<",";
		GenomeFitness = Fitness;
	}
	//////////////////////////////////
	friend void Crossover(int* GenomePointer1, int* GenomePointer2, int Cut);   //one point交配
	friend void Crossover(int* GenomePointer1, int* GenomePointer2, int OneCut, int TwoCut);   //two point交配
};
void Crossover(int* GenomePointer1, int* GenomePointer2, int* GenomePointer3, int* GenomePointer4, int Cut)//one point交配
{
	if (GenomePointer1 != GenomePointer2)
	{
		for (int i = 0; i<NumofLight; i++)
		{
			*(GenomePointer3 + i) = *(GenomePointer1 + i);
			*(GenomePointer4 + i) = *(GenomePointer2 + i);
			if (i >= Cut)
			{
				*(GenomePointer4 + i) = *(GenomePointer1 + i);
				*(GenomePointer3 + i) = *(GenomePointer2 + i);
			}
		}
	}
}
void Crossover(int* GenomePointer1, int* GenomePointer2, int* GenomePointer3, int* GenomePointer4, int OneCut, int TwoCut)//two point交配，從OneCut換到TwoCut，包含前後
{
	if (GenomePointer1 != GenomePointer2)
	{
		if (OneCut>TwoCut)
		{
			int CutChange = OneCut;
			OneCut = TwoCut;
			TwoCut = CutChange;
		}
		for (int i = 0; i<NumofLight; i++)
		{
			*(GenomePointer3 + i) = *(GenomePointer1 + i);
			*(GenomePointer4 + i) = *(GenomePointer2 + i);
			if (i >= OneCut && i <= TwoCut)
			{
				*(GenomePointer4 + i) = *(GenomePointer1 + i);
				*(GenomePointer3 + i) = *(GenomePointer2 + i);
			}
		}
	}
}
GeneAssemble::GeneAssemble(){}
GeneAssemble::GeneAssemble(int* GenomePointer, float Fitness)
{
	for (int i = 0; i<NumofLight; i++)
	{
		Genome[i] = *(GenomePointer + i);
	}
	GenomeFitness = Fitness;
}
//GA類別結束
void quickSort(GeneAssemble* GenomeGroup[], int left, int right);

void Parameter()
{
	CircleofGA = NumofGenome*(0.1);      //GA執行的循環次數，基因組數數的10%
	if (CircleofGA<10)   //世代次數過少會造成無解
		CircleofGA = 10;

	NumofCrossoverGenome = NumofGenome / 2;//排名較前面，用於交配的基因組,假設50%

	NumofCPUCircle = NumofCPU;    //各處理器間交流所需，約為處理器個數
}

int main()
{
	char FileInput[500];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofInstruction
		for (int i = 0; i<NumofInstruction; i++)
		{
			for (int j = 0; j<2; j++)
			{
				fr.getline(FileInput, sizeof(FileInput), ',');
				stringstream(FileInput) >> I[i][j];
				//sscanf(FileInput, "%f", &I[i][j]);        //字串轉數字
				//cout << I[i][j] << ", ";
			}
			//cout << endl;
		}
		//system("pause");
	}
	fr.close();

	//Td
	fr.open("Td.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput), ',');//這是NumofTd
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput), ',');
			stringstream(FileInput) >> Td[i];
			//sscanf(FileInput, "%f", &Td[i]);        //字串轉數字
		}
	}
	fr.close();

	//KParameter
	string InputProcess;
	fr.open("KParameter.txt", ios::in);
	if (!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
		cout << "Fail to open file" << endl;
	else
	{
		fr.getline(FileInput, sizeof(FileInput));//矩陣高、寬
		//fr.getline(FileInput, sizeof(FileInput), ',');//矩陣寬
		for (int i = 0; i<NumofTd; i++)
		{
			fr.getline(FileInput, sizeof(FileInput));
			int j = 0;
			/**
			for (int j = 0; j<NumofLight; j++)
			{
			fr.getline(FileInput, sizeof(FileInput), ',');
			stringstream(FileInput) >> KParameter[i][j];
			//sscanf(FileInput, "%f", &KParameter[i][j]);        //字串轉數字
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



	Parameter();  //係數設定

	clock_t start, stop;
	start = clock(); //開始時間
	float EnergyResume;
	int BestTimes = 0;
	//int SecTimes = 0;
	//int ThirdTimes = 0;
	for (int i = 0; i<RunTime; i++)
	{
		EnergyResume = GA();
		if (EnergyResume <= RealBest+0.001)
		{
			BestTimes++;
		}
	}
	cout << "總執行次數: " << RunTime << "其中:"<< endl;
	cout << "最佳解次數: " << BestTimes << endl;
	cout << "平均耗能:" << ave_power / RunTime << endl;
	stop = clock(); //結束時間

	cout << "總計費時:" << double(stop - start) / CLOCKS_PER_SEC << "秒" << endl;
	cout << "平均費時:" << (double(stop - start) / CLOCKS_PER_SEC) / RunTime << "秒" << endl;



	system("pause");
	return 0;
}

//確認亮度合格與否
int CheckT(int MatrixRow, int i[])     //MatrixRow為要檢查是否合格的row，i為紀錄指令的矩陣
{
	float T = 0;
	for (int k = 0; k<(NumofLight); k++)
	{
		if (k >= NumofLight)
		{
			//T+=KParameter[MatrixRow][k]*OtherLightT[k-NumofLight];
		}
		else
		{
			T += KParameter[MatrixRow][k] * I[i[k]][1];//[表示燈的指令,1表示亮度]nbv
		}
	}
	if (T<Td[MatrixRow])
		return 0;
	else
		return 1;
}


float GA()
{
	float Power = 0;
	int GenomeNo = 0;                  //目標是找到NumofGene個解作為初始值，所以設一個旗標用以計數
	int Genome[NumofLight];                         //基因組[第一盞燈的指令][第二盞燈的指令]..[最後一盞燈的指令]
	float T = 0;
	GeneAssemble* GenomeGroup[NumofCPU][NumofGenome];        //眾多基因組的集合
	GeneAssemble* ExchangeTemp[NumofExchangeGenome];//各群交流時使用的暫存
	GeneAssemble* ExchangeTemp2[NumofExchangeGenome];//各群交流時使用的暫存
	GeneAssemble* Best;
	///////
	clock_t start, stop;
	start = clock(); //開始時間

	for (int CPUNo = 0; CPUNo<NumofCPU; CPUNo++)
	{
		GenomeNo = 0;
		///////
		for (int i = 0; i<NumofGenome; i++)
		{
			int j[NumofLight];

			for (int h = 0; h < NumofLight; h++)                             //初始基因組隨機產生方式
			{                                                                                                                                                    //cout<<StartPointCheck<<",";
				j[h] = RandomNumber(0, NumofInstruction - 1);
			}
			/**
			if (i == 0)
			{
			j[0] = 3;
			j[1] = 15;
			j[2] = 6;
			j[3] = 15;
			j[4] = 15;
			j[5] = 15;
			j[6] = 15;
			j[7] = 15;
			j[8] = 15;
			j[9] = 15;
			j[10] = 15;
			j[11] = 8;
			}/**/
			//////////  
			for (int MatrixRow = 0; MatrixRow<NumofTd; MatrixRow++)
			{
				T = CheckT(MatrixRow, j);
				if (T == 0)                                                     //不合格，Power設Impossible，不可能發生的情況，作為初值
				{
					for (int p = 0; p<NumofLight; p++)                            //基因組產生
					{
						Genome[p] = j[p];
					}
					GenomeGroup[CPUNo][GenomeNo++] = new GeneAssemble(&Genome[0], Impossible);                         //基因組加入集合
					break;
				}
				else
				{
					T = 0;
				}
				if (MatrixRow == NumofTd - 1)//全部Row都符合條件，開始計算耗能
				{                                                                        //cout<<"前"<<Power<<endl;system("pause");
					for (int k = 0; k<NumofLight; k++)                                                 //計算耗能
					{
						Genome[k] = j[k];   //cout<<"m:"<<m<<",";  
						Power += I[Genome[k]][0];                                                       //cout<<Power<<endl;system("pause");
					}                                                                         //cout<<"後"<<Power<<endl;system("pause");
					GenomeGroup[CPUNo][GenomeNo++] = new GeneAssemble(&Genome[0], Power);                //產生基因組
					Power = 0;
				}
			}
		}
		quickSort(GenomeGroup[CPUNo], 0, NumofGenome - 1);      //給GenomeGroup[]作由小到大的排序
	}
	//初始值產生，第一輪結束
	//還要再執行CircleofGA-1次
	//新增平行處理，在最外層
	////
	for (int CPUCircle = 0; CPUCircle < NumofCPUCircle; CPUCircle++)
	{
		for (int CPUNo = 0; CPUNo < NumofCPU; CPUNo++)
		{
			for (int i = 1; i < CircleofGA; i++)
			{
				for (int j = NumofCrossoverGenome + 1; j < NumofGenome - 1; j += 2)
				{
					//隨機選擇one point或two point的crossover，再隨機決定哪兩個基因組作交配，最後隨機選擇要交換的位置(one point一個，two point兩個)
					if (RandomNumber(0, 1))
					{
						Crossover(GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][j]->GetGenomePointer(), GenomeGroup[CPUNo][j + 1]->GetGenomePointer(), RandomNumber(1, NumofLight - 1));
					}
					else
					{
						Crossover(GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][j]->GetGenomePointer(), GenomeGroup[CPUNo][j + 1]->GetGenomePointer(), RandomNumber(0, NumofLight - 1), RandomNumber(0, NumofLight - 1));
					}
					//Crossover(第一基因組,第二基因組,要被蓋掉的第三基因組,要被蓋掉的第四基因組,Cut);
					//Crossover(第一基因組,第二基因組,要被蓋掉的第三基因組,要被蓋掉的第四基因組,OneCut,TwoCut);
				}
				for (int j = 0; j < NumofGenome; j++)
				{
					if (RandomNumber(1, MutationPossibility) == 1)
						GenomeGroup[CPUNo][j]->mutation(RandomNumber(0, NumofLight - 1), RandomNumber(0, NumofInstruction - 1));
				}
				for (int j = 0; j < NumofGenome; j++)
				{
					GenomeGroup[CPUNo][j]->CheckGenome();
				}
				quickSort(GenomeGroup[CPUNo], 0, NumofGenome - 1);      //給GenomeGroup[]作由小到大的排序
			}
		}

		if (NumofCPU != 1)//超過一個群才要交流
		{
			//各群間的交流(對象決定)
			int RandomChange[NumofCPU][2];   //交流對象表格(是否已有接收基因組,基因組交付對像)
			for (int i = 0; i < NumofCPU; i++)
			{
				RandomChange[i][0] = -1;    //-1表示尚未接收
			}
			int RandomNo;        //隨機選擇交流對象
			for (int i = 0; i < NumofCPU; i++)     //第i群的基因組交給第幾群
			{
				RandomNo = RandomNumber(0, 9);      //得到第i群基因組的是第幾群
				while (RandomChange[RandomNo][0] != -1)
				{
					RandomNo = RandomNumber(0, 9);
				}
				RandomChange[RandomNo][0] = i;
				RandomChange[i][1] = RandomNo;
			}

			//開始交流
			int Flag = 0;//交流的進度旗桿
			for (int i = 0; i < NumofExchangeGenome; i++)
			{
				ExchangeTemp[i] = GenomeGroup[RandomChange[Flag][1]][i];
				GenomeGroup[RandomChange[Flag][1]][i] = GenomeGroup[Flag][i];
			}
			Flag = RandomChange[Flag][1];
			while (Flag != 0)
			{
				for (int i = 0; i < NumofExchangeGenome; i++)
				{
					ExchangeTemp2[i] = GenomeGroup[RandomChange[Flag][1]][i];
					GenomeGroup[RandomChange[Flag][1]][i] = ExchangeTemp[i];
				}
				Flag = RandomChange[Flag][1];
			}
		}
		//最佳基因組尋找
		for (int i = 0; i<NumofCPU; i++)
		{
			if (i == 0)
				Best = GenomeGroup[i][0];//經過排序，第0個為最佳
			if (Best->GetGenomeFitness() > GenomeGroup[i][0]->GetGenomeFitness())
			{
				Best = GenomeGroup[i][0];//經過排序，第0個為最佳
			}
		}

	}
	ave_power += Best->GetGenomeFitness();
	stop = clock(); //結束時間
	
	////
	//GA執行結束，觀看結果
	fstream fw;
	string filename = "Specific Solution_best.txt";
	fw.open(filename, ios::out);
	//cout<<"GA搜尋"<<endl;
	/**/
	cout << "指令:";
	int* BestGenome = Best->GetGenomePointer();               //經過排序，所以最佳結果在GenomeGroup[0]
	for (int j = 0; j<NumofLight; j++)
	{
		fw << *(BestGenome + j) << ",";
		cout << *(BestGenome + j) << ",";
	}
	fw << "耗能:" << Best->GetGenomeFitness() << endl;
	cout << "耗能:" << Best->GetGenomeFitness() << ", ";
	cout << "費時:" << double(stop - start) / CLOCKS_PER_SEC << "秒" << endl;
	fw.close();
	/**/
	return Best->GetGenomeFitness();

}


void quickSort(GeneAssemble* GenomeGroup[], int left, int right)
{
	int i = left, j = right;
	GeneAssemble* SortingTmp;
	float Pivot = GenomeGroup[(left + right) / 2]->GetGenomeFitness();

	// partition 
	while (i <= j)
	{
		float TempFitness;
		TempFitness = GenomeGroup[i]->GetGenomeFitness();
		while (TempFitness < Pivot)  //這個排序中，將0視為最大
		{
			i++;
			TempFitness = GenomeGroup[i]->GetGenomeFitness();
		}
		TempFitness = GenomeGroup[j]->GetGenomeFitness();
		while (TempFitness > Pivot)
		{
			j--;
			TempFitness = GenomeGroup[j]->GetGenomeFitness();
		}
		if (i <= j)
		{
			SortingTmp = GenomeGroup[i];
			GenomeGroup[i] = GenomeGroup[j];
			GenomeGroup[j] = SortingTmp;
			i++;
			j--;
		}
	};
	// recursion 
	if (left < j)
		quickSort(GenomeGroup, left, j);
	if (i < right)
		quickSort(GenomeGroup, i, right);
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