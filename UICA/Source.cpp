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
#define OtherLight 0        //額外光源為太陽
//float OtherLightT[OtherLight]={0};     //額外光源的照度

#define BalanceRatio 5 //多國拉鋸平衡，國力差距5%內

float the_best = 224.22;

int RunTime = 50;  //跑多少次
#define ImpossiblePowerConsume 9999    //照度不合格，給予不可能的極大能量消耗，使其容易被照度合格者取代

#define NumofCountry 100000    //取NumofCountry組解作為初值(國家數)
#define NumofEmpire 10     //帝國數
#define VelocityLimit 3
#define InitialW 0.4             //慣性係數，由大到小
#define FinalW 0.4
#define InitialCOne 2          //本身的C係數，由大到小
#define FinalCOne 2
#define InitialCTwo  2       //群體的C係數，由小到大
#define FinalCTwo 2


float I[NumofInstruction][2];   //(耗能,亮度)
float Td[NumofTd];
float KParameter[NumofTd][NumofLight + OtherLight];

class Country  //每個粒子的基本設定，變動調光指令版本
{
private:
	float PBest[NumofLight + 1],  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,個體最低耗能)
		GBest[NumofLight + 1];  //(第一盞燈指令,第二盞燈指令,第三盞燈指令,群體最低耗能)
	int Position[NumofLight];  //記錄個體，其各盞燈光的調光指令 <----> Position   //其中燈光數也是solution的維度，第幾盞燈表示第幾維
	int Velocity[NumofLight];             //調光指令上升或下降 <----> Velocity
	float W;//模擬慣性用的weighting，其值越大，越適合全域搜尋，建議由0.9~0.4
	float C1, C2;//係數，調整偏向各體，或偏向群體
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
		for (int i = 0; i<NumofLight; i++)//i是維度
		{                                                                                               //cout<<"原P"<<i<<":"<<Position[i]<<",";
			Position[i] = Position[i] + Velocity[i];                                                            //cout<<"V:"<<Velocity[i]<<endl;
			if (Position[i]>NumofInstruction - 1)             //調光指令的上下限
				Position[i] = NumofInstruction - 1;
			if (Position[i]<0)
				Position[i] = 0;                                                                           //cout<<"後P"<<i<<":"<<Position[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void UpdateVelocity()
	{
		for (int i = 0; i<NumofLight; i++)//i是維度
		{                                                                                               //cout<<"原V"<<i<<":"<<Velocity[i]<<",";
			Velocity[i] = W*Velocity[i] + C1*(PBest[i] - Position[i]) + C2*(GBest[i] - Position[i]);         //cout<<"後V"<<i<<":"<<Velocity[i]<<endl;
			if (Velocity[i]>VelocityLimit)                   //速度上下限
				Velocity[i] = VelocityLimit;
			if (Velocity[i]<(-VelocityLimit))
				Velocity[i] = (-VelocityLimit);                                                            //cout<<"後V"<<i<<":"<<Velocity[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
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
//Country類別結束

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
	clock_t start, stop;
	start = clock(); //開始時間
	float AverageResult = 0;
	float EnergyResume;//檢察解解
	int NoAnswer = 0;//無解總次數
	int BestTimes = 0;
	int SecTimes = 0;
	int ThirdTimes = 0;

	int NotBest = 0;//非最佳解的次數
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
	cout << "最佳解個數" << BestTimes << endl;
	//cout << endl << "非最佳解的次數為" << NotBest << endl;
	stop = clock(); //結束時間
	cout << "執行" << RunTime << "次" << endl;
	cout << "平均耗能:" << (AverageResult / RunTime) << "秒" << endl;
	cout << "合計費時:" << double(stop - start) / CLOCKS_PER_SEC << "秒" << endl << endl;
	cout << "無解次數:" << NoAnswer << "次" << endl;
	system("pause");
	return 0;
}

float Empire()
{
	float ENDFlag = ImpossiblePowerConsume;//檢查剩餘兩國時，是否進入平衡狀態
	int CirCleOrNot = 0;//是否進入循環

	int remain = -1;
	int remain_count = 0;
	int flag = 0;

	Country* AllCountry[NumofCountry];
	int RandomVelocity[NumofLight];
	clock_t start, stop;
	start = clock(); //開始時間

	for (int CountryNo = 0; CountryNo<NumofCountry; CountryNo++)
	{
		int CountryConstruct[NumofLight];

		for (int h = 0; h<NumofLight; h++)
		{
			CountryConstruct[h] = RandomNumber(0, NumofInstruction - 1);                           //cout<<CountryConstruct[h]<<",";
		}

		//初始位置產生完畢

		for (int MatrixRow = 0; MatrixRow<NumofTd; MatrixRow++)
		{
			int T = CheckT(MatrixRow, CountryConstruct);
			if (T == 0)                                                     //不合格，Power設ImpossiblePowerConsume，不可能發生的情況，作為初值
			{                                                                              //cout<<"ImpossiblePower"<<","; 
				for (int h = 0; h<NumofLight; h++)
				{
					RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
				}
				AllCountry[CountryNo] = new Country(CountryConstruct, ImpossiblePowerConsume, RandomVelocity);
				break;
			}
			if (MatrixRow == NumofTd - 1)//全部Row都符合條件，開始計算耗能
			{
				float Power = 0;
				for (int k = 0; k<NumofLight; k++)                                                 //計算耗能
				{
					Power += I[CountryConstruct[k]][0];
				}                                                                             //cout<<Power<<",";
				for (int h = 0; h<NumofLight; h++)                                                 //隨機初始速度
				{
					RandomVelocity[h] = (RandomNumber(0, 1) ? 1 : -1)*RandomNumber(0, VelocityLimit);
				}
				AllCountry[CountryNo] = new Country(CountryConstruct, Power, RandomVelocity);                //國家產生
			}
		}
	}
	//依照國家成本排序，成本越小，國力越強	
	quickSort(AllCountry, 0, NumofCountry - 1);    //(data,left,right)

	Country* EmpireCountry[NumofEmpire];
	Country* ColonyCountry[NumofEmpire][NumofCountry - NumofEmpire];   //[所屬帝國聯邦編號][被殖民國編號]，一個帝國聯邦，擁有哪些殖民國
	float EmpireCost[NumofEmpire];
	float EmpirePower[NumofEmpire];
	float TotalCost = 0;
	int EmpireNumofColony[NumofEmpire];  //一個帝國聯邦所能擁有的殖民國數目
	for (int i = 0; i<NumofEmpire; i++)
	{
		EmpireCountry[i] = AllCountry[i];
		EmpireCost[i] = (EmpireCountry[i]->GetFitness()) - (AllCountry[NumofEmpire]->GetFitness());  //正規化成本:Cn=成本:cn-最大成本ci:max{ci}
		TotalCost += EmpireCost[i];
	}
	for (int i = 0; i<NumofEmpire; i++)
	{
		EmpirePower[i] = EmpireCost[i] / TotalCost;
		EmpireNumofColony[i] = EmpirePower[i] * (NumofCountry - NumofEmpire);
	}
	int Remain = NumofCountry - NumofEmpire;      //計算是否全部國家都有殖民與統治關係
	for (int i = 0; i<NumofEmpire; i++)
	{
		Remain = Remain - EmpireNumofColony[i];
	}
	for (int i = 0, j = 0; i<Remain; i++)         //未除盡的國家數，要補進
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
	int ColonyNo[NumofEmpire];                  //加入帝國聯邦的順序編號，也用於計算目前有多少殖民國,[第幾帝國聯邦]
	for (int j = 0; j<NumofEmpire; j++)
	{
		ColonyNo[j] = 0;
	}
	for (int i = NumofEmpire; i<NumofCountry; i++)
	{
		int EmpireSelection = RandomNumber(0, (NumofEmpire - 1));
		if (ColonyNo[EmpireSelection]<EmpireNumofColony[EmpireSelection])   //檢查選擇的帝國聯邦是否尚未額滿
		{
			ColonyCountry[EmpireSelection][ColonyNo[EmpireSelection]] = AllCountry[i];
			ColonyNo[EmpireSelection] = ColonyNo[EmpireSelection] + 1;
		}
		else                                                              //若是額滿的話，需要重新選擇
		{
			i--;
		}
	}
	//初始輪結束
	float BalanceCost = 0;
	int Counter = 0;
	int Collapse;
	//第二輪開始
	while (1)         //當僅存單依一帝國聯邦，透過最後的if中的break跳出
	{
		for (int j = 0; j<NumofEmpire; j++)                 //依帝國聯邦為單位處理
		{
			for (int k = 0; k<ColonyNo[j]; k++)            //帝國聯邦中的殖民國依序處理
			{
				ColonyCountry[j][k]->SetGBest(EmpireCountry[j]->GetPBest());
				ColonyCountry[j][k]->UpdateVelocity();
				ColonyCountry[j][k]->UpdatePosition();
				ColonyCountry[j][k]->CheckPosition();
			}
		}

		for (int j = 0; j<NumofEmpire; j++)                 //依帝國聯邦為單位處理
		{
			for (int k = 0; k<ColonyNo[j]; k++)            //帝國聯邦中的殖民國依序處理
			{
				if ((ColonyCountry[j][k]->GetFitness()) < EmpireCountry[j]->GetFitness())   //殖民國竄位，成為帝國核心
				{
					Country* SwitchingPointer = ColonyCountry[j][k];
					ColonyCountry[j][k] = EmpireCountry[j];
					EmpireCountry[j] = SwitchingPointer;
				}
			}
		}
		//強大帝國搶最弱小國的殖民國
		int WeakestEmpireNo = -1;
		int AwfulEmpireNo = -1;
		float TotalColonyCost = 0;
		//////////////////////////

		for (int j = 0; j<NumofEmpire; j++)
		{
			//計算各國強大程度
			if (ColonyNo[j] != 0)
			{
				TotalColonyCost = 0;
				for (int k = 0; k<ColonyNo[j]; k++)
				{
					TotalColonyCost += ColonyCountry[j][k]->GetFitness();
				}
				TotalColonyCost = (TotalColonyCost) / ColonyNo[j];
			}
			EmpireCost[j] = ((EmpireCountry[j]->GetFitness() * 7) + (TotalColonyCost * 3)) / 10;//(帝國權力*7+平均殖民國權力*3)/10
		}

		//拉鋸戰時的平衡局面破除
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
			if ((AwfulEmpireNo == -1 || EmpireCost[j] < EmpireCost[AwfulEmpireNo]))   //尋找最強大帝國(成本最小)
			{
				AwfulEmpireNo = j;
			}
			int CheckExist = 0;//剩餘殖民地個數

			{
				CheckExist += ColonyNo[j];
			}
			if ((WeakestEmpireNo == -1 || EmpireCost[j] > EmpireCost[WeakestEmpireNo]) && (CheckExist != 0))   //尋找現存最弱小帝國(成本最大)
			{
				WeakestEmpireNo = j;
			}
		}

		///////////////////////////
		////////////////////////////

		//尋找最弱小國的最弱小殖民地
		Country* WeakestColony = NULL;
		int WeakestColonyNo;
		for (int j = 0; j<ColonyNo[WeakestEmpireNo]; j++)
		{
			if (WeakestColony == NULL || (ColonyCountry[WeakestEmpireNo][j]->GetFitness()) >(WeakestColony->GetFitness()))  //尋找最弱小殖民國(成本最大)
			{
				WeakestColony = ColonyCountry[WeakestEmpireNo][j];
				WeakestColonyNo = j;
			}
		}

		//最弱小帝國的最弱小殖民國被搶
		//被誰搶
		int Marauder = -1;
		TotalCost = 0;
		int MarauderPassibility;
		float WeakestCost = EmpireCost[WeakestEmpireNo];
		for (int j = 0; j<NumofEmpire; j++)
		{
			int CheckExist = 0;//剩餘殖民地個數

			{
				CheckExist += ColonyNo[j];
			}
			if (CheckExist != 0)
			{
				EmpireCost[j] = EmpireCost[j] - WeakestCost;      //正規化成本
				TotalCost += EmpireCost[j];
			}
		}
		while (Marauder == WeakestEmpireNo || Marauder == -1)
		{
			MarauderPassibility = RandomNumber(1, 100);
			for (int j = 0; j<NumofEmpire; j++)
			{
				int CheckExist = 0;//剩餘殖民地個數

				CheckExist += ColonyNo[j];

				if (CheckExist != 0)
				{
					EmpirePower[j] = EmpireCost[j] / TotalCost;                      //權力計算
					//是不是這個帝國搶到?
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
		//被搶了

		ColonyCountry[Marauder][ColonyNo[Marauder]] = ColonyCountry[WeakestEmpireNo][WeakestColonyNo];
		//弱小國殖民地整理
		if (WeakestColonyNo != (ColonyNo[WeakestEmpireNo] - 1))//最弱小的殖民地，不是編號最後面的殖民地時
		{
			ColonyCountry[WeakestEmpireNo][WeakestColonyNo] = ColonyCountry[WeakestEmpireNo][ColonyNo[WeakestEmpireNo] - 1];
			ColonyCountry[WeakestEmpireNo][ColonyNo[WeakestEmpireNo] - 1] = NULL;
		}
		else
		{
			ColonyCountry[WeakestEmpireNo][WeakestColonyNo] = NULL;
		}
		ColonyNo[WeakestEmpireNo]--;//殖民地少一個
		//
		ColonyNo[Marauder]++;//殖民地多一個


		/////////////////////////////
		Collapse = 0;
		for (int j = 0; j<NumofEmpire; j++)
		{

			if (ColonyNo[j] == 0)
			{
				Collapse++;
			}
		}
		if (Collapse >= (NumofEmpire - 1))  //剩餘一國
		{
			break;
		}

		else if (Collapse >= (NumofEmpire - 2))  //剩餘2國，國力極接近，進入平衡狀態，強弱關係的週期為2輪
		{
			if (CirCleOrNot == 0)
			{
				//cout<<"強聯邦:"<<EmpireCost[AwfulEmpireNo]<<endl;
				//cout<<"弱聯邦:"<<EmpireCost[WeakestEmpireNo]<<endl;
				ENDFlag = EmpireCost[AwfulEmpireNo];//強國弱國不斷對調，所以特定國家成為強國的週期也是兩輪
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
		//多國拉鋸
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

		//觀看戰局
		Counter++;
		if (Counter >= 100)
		{
			flag = 1;
		}
		if (Counter == 1000)
		{
			cout << "執行中,剩餘" << NumofEmpire - Collapse << "帝國,最強:" << EmpireCost[AwfulEmpireNo] << ",最弱:" << EmpireCost[WeakestEmpireNo] << endl;
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
	cout << "指令:";
	for (int j = 0; j<NumofLight; j++)
	{
		cout << *(EmpireCountry[BetterEmpireNo]->GetPBest() + j) << ",";
	}
	cout << "耗能:" << EmpireCountry[BetterEmpireNo]->GetFitness()<<", ";

	stop = clock(); //結束時間
	cout << "費時:" << double(stop - start) / CLOCKS_PER_SEC << "秒" << endl;
	//system("pause");
	/********************/

	return EmpireCountry[BetterEmpireNo]->GetFitness();
}

//確認亮度合格與否
int CheckT(int MatrixRow, int i[])     //MatrixRow為要檢查是否合格的row，i為紀錄指令的矩陣
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
			T += KParameter[MatrixRow][k] * I[i[k]][1];//[表示燈的指令,1表示亮度]
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
給予希望的最大值與最小值，就能產生介於這兩個值中的隨機數，這個隨機亂數也可能等於最大值，或等於最小值。
其亂數值是依照時間所產生，所以重複呼叫時，會給予不同的亂數值(整數)。
***************************************************************************************************/
int RandomNumber(int MinValue, int MaxValue)   //Both MinValue and MaxValue are included
{
	int R = (rand() % (MaxValue - MinValue + 1)) + MinValue;
	return R;
}