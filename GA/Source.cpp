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
#define Impossible 99999
#define MutationPossibility 10000
int RunTime = 50;
float RealBest = 224.22;

float ave_power = 0;

#define NumofGenome 100000
#define NumofCPU 1
#define NumofExchangeGenome 10000  //��y����]�ը�10%


int NumofCrossoverGenome;  //�ƦW���e���A�Ω��t����]��,���]50%
int CircleofGA;  //�򥻦��ơA�H���I�W�ɡA�ݭn�y�L�W�[�A���v�T���j
int NumofCPUCircle;

float I[NumofInstruction][2];   //(�ӯ�,�G��)
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
	void mutation(int Which, int What)                        //Which:���@�Ӱ�]���ܡAWhat:�ܦ�����
	{
		Genome[Which] = What;
	}
	//////////////////////////////////
	void CheckGenome()
	{
		float Tga = 0;  //�ھڥثe��]�թҺ�X�Ӫ�T(�ӫ�)��
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
	friend void Crossover(int* GenomePointer1, int* GenomePointer2, int Cut);   //one point��t
	friend void Crossover(int* GenomePointer1, int* GenomePointer2, int OneCut, int TwoCut);   //two point��t
};
void Crossover(int* GenomePointer1, int* GenomePointer2, int* GenomePointer3, int* GenomePointer4, int Cut)//one point��t
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
void Crossover(int* GenomePointer1, int* GenomePointer2, int* GenomePointer3, int* GenomePointer4, int OneCut, int TwoCut)//two point��t�A�qOneCut����TwoCut�A�]�t�e��
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
//GA���O����
void quickSort(GeneAssemble* GenomeGroup[], int left, int right);

void Parameter()
{
	CircleofGA = NumofGenome*(0.1);      //GA���檺�`�����ơA��]�ռƼƪ�10%
	if (CircleofGA<10)   //�@�N���ƹL�ַ|�y���L��
		CircleofGA = 10;

	NumofCrossoverGenome = NumofGenome / 2;//�ƦW���e���A�Ω��t����]��,���]50%

	NumofCPUCircle = NumofCPU;    //�U�B�z������y�һݡA�����B�z���Ӽ�
}

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
				//cout << I[i][j] << ", ";
			}
			//cout << endl;
		}
		//system("pause");
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



	Parameter();  //�Y�Ƴ]�w

	clock_t start, stop;
	start = clock(); //�}�l�ɶ�
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
	cout << "�`���榸��: " << RunTime << "�䤤:"<< endl;
	cout << "�̨θѦ���: " << BestTimes << endl;
	cout << "�����ӯ�:" << ave_power / RunTime << endl;
	stop = clock(); //�����ɶ�

	cout << "�`�p�O��:" << double(stop - start) / CLOCKS_PER_SEC << "��" << endl;
	cout << "�����O��:" << (double(stop - start) / CLOCKS_PER_SEC) / RunTime << "��" << endl;



	system("pause");
	return 0;
}

//�T�{�G�צX��P�_
int CheckT(int MatrixRow, int i[])     //MatrixRow���n�ˬd�O�_�X�檺row�Ai���������O���x�}
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
			T += KParameter[MatrixRow][k] * I[i[k]][1];//[��ܿO�����O,1��ܫG��]nbv
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
	int GenomeNo = 0;                  //�ؼЬO���NumofGene�Ӹѧ@����l�ȡA�ҥH�]�@�ӺX�ХΥH�p��
	int Genome[NumofLight];                         //��]��[�Ĥ@���O�����O][�ĤG���O�����O]..[�̫�@���O�����O]
	float T = 0;
	GeneAssemble* GenomeGroup[NumofCPU][NumofGenome];        //���h��]�ժ����X
	GeneAssemble* ExchangeTemp[NumofExchangeGenome];//�U�s��y�ɨϥΪ��Ȧs
	GeneAssemble* ExchangeTemp2[NumofExchangeGenome];//�U�s��y�ɨϥΪ��Ȧs
	GeneAssemble* Best;
	///////
	clock_t start, stop;
	start = clock(); //�}�l�ɶ�

	for (int CPUNo = 0; CPUNo<NumofCPU; CPUNo++)
	{
		GenomeNo = 0;
		///////
		for (int i = 0; i<NumofGenome; i++)
		{
			int j[NumofLight];

			for (int h = 0; h < NumofLight; h++)                             //��l��]���H�����ͤ覡
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
				if (T == 0)                                                     //���X��APower�]Impossible�A���i��o�ͪ����p�A�@�����
				{
					for (int p = 0; p<NumofLight; p++)                            //��]�ղ���
					{
						Genome[p] = j[p];
					}
					GenomeGroup[CPUNo][GenomeNo++] = new GeneAssemble(&Genome[0], Impossible);                         //��]�ե[�J���X
					break;
				}
				else
				{
					T = 0;
				}
				if (MatrixRow == NumofTd - 1)//����Row���ŦX����A�}�l�p��ӯ�
				{                                                                        //cout<<"�e"<<Power<<endl;system("pause");
					for (int k = 0; k<NumofLight; k++)                                                 //�p��ӯ�
					{
						Genome[k] = j[k];   //cout<<"m:"<<m<<",";  
						Power += I[Genome[k]][0];                                                       //cout<<Power<<endl;system("pause");
					}                                                                         //cout<<"��"<<Power<<endl;system("pause");
					GenomeGroup[CPUNo][GenomeNo++] = new GeneAssemble(&Genome[0], Power);                //���Ͱ�]��
					Power = 0;
				}
			}
		}
		quickSort(GenomeGroup[CPUNo], 0, NumofGenome - 1);      //��GenomeGroup[]�@�Ѥp��j���Ƨ�
	}
	//��l�Ȳ��͡A�Ĥ@������
	//�٭n�A����CircleofGA-1��
	//�s�W����B�z�A�b�̥~�h
	////
	for (int CPUCircle = 0; CPUCircle < NumofCPUCircle; CPUCircle++)
	{
		for (int CPUNo = 0; CPUNo < NumofCPU; CPUNo++)
		{
			for (int i = 1; i < CircleofGA; i++)
			{
				for (int j = NumofCrossoverGenome + 1; j < NumofGenome - 1; j += 2)
				{
					//�H�����one point��two point��crossover�A�A�H���M�w����Ӱ�]�է@��t�A�̫��H����ܭn�洫����m(one point�@�ӡAtwo point���)
					if (RandomNumber(0, 1))
					{
						Crossover(GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][j]->GetGenomePointer(), GenomeGroup[CPUNo][j + 1]->GetGenomePointer(), RandomNumber(1, NumofLight - 1));
					}
					else
					{
						Crossover(GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][RandomNumber(0, NumofCrossoverGenome)]->GetGenomePointer(), GenomeGroup[CPUNo][j]->GetGenomePointer(), GenomeGroup[CPUNo][j + 1]->GetGenomePointer(), RandomNumber(0, NumofLight - 1), RandomNumber(0, NumofLight - 1));
					}
					//Crossover(�Ĥ@��]��,�ĤG��]��,�n�Q�\�����ĤT��]��,�n�Q�\�����ĥ|��]��,Cut);
					//Crossover(�Ĥ@��]��,�ĤG��]��,�n�Q�\�����ĤT��]��,�n�Q�\�����ĥ|��]��,OneCut,TwoCut);
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
				quickSort(GenomeGroup[CPUNo], 0, NumofGenome - 1);      //��GenomeGroup[]�@�Ѥp��j���Ƨ�
			}
		}

		if (NumofCPU != 1)//�W�L�@�Ӹs�~�n��y
		{
			//�U�s������y(��H�M�w)
			int RandomChange[NumofCPU][2];   //��y��H���(�O�_�w��������]��,��]�ե�I�ﹳ)
			for (int i = 0; i < NumofCPU; i++)
			{
				RandomChange[i][0] = -1;    //-1��ܩ|������
			}
			int RandomNo;        //�H����ܥ�y��H
			for (int i = 0; i < NumofCPU; i++)     //��i�s����]�ե浹�ĴX�s
			{
				RandomNo = RandomNumber(0, 9);      //�o���i�s��]�ժ��O�ĴX�s
				while (RandomChange[RandomNo][0] != -1)
				{
					RandomNo = RandomNumber(0, 9);
				}
				RandomChange[RandomNo][0] = i;
				RandomChange[i][1] = RandomNo;
			}

			//�}�l��y
			int Flag = 0;//��y���i�׺X��
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
		//�̨ΰ�]�մM��
		for (int i = 0; i<NumofCPU; i++)
		{
			if (i == 0)
				Best = GenomeGroup[i][0];//�g�L�ƧǡA��0�Ӭ��̨�
			if (Best->GetGenomeFitness() > GenomeGroup[i][0]->GetGenomeFitness())
			{
				Best = GenomeGroup[i][0];//�g�L�ƧǡA��0�Ӭ��̨�
			}
		}

	}
	ave_power += Best->GetGenomeFitness();
	stop = clock(); //�����ɶ�
	
	////
	//GA���浲���A�[�ݵ��G
	fstream fw;
	string filename = "Specific Solution_best.txt";
	fw.open(filename, ios::out);
	//cout<<"GA�j�M"<<endl;
	/**/
	cout << "���O:";
	int* BestGenome = Best->GetGenomePointer();               //�g�L�ƧǡA�ҥH�̨ε��G�bGenomeGroup[0]
	for (int j = 0; j<NumofLight; j++)
	{
		fw << *(BestGenome + j) << ",";
		cout << *(BestGenome + j) << ",";
	}
	fw << "�ӯ�:" << Best->GetGenomeFitness() << endl;
	cout << "�ӯ�:" << Best->GetGenomeFitness() << ", ";
	cout << "�O��:" << double(stop - start) / CLOCKS_PER_SEC << "��" << endl;
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
		while (TempFitness < Pivot)  //�o�ӱƧǤ��A�N0�����̤j
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
�����Ʊ檺�̤j�ȻP�̤p�ȡA�N�ಣ�ͤ���o��ӭȤ����H���ơA�o���H���üƤ]�i�൥��̤j�ȡA�ε���̤p�ȡC
��üƭȬO�̷Ӯɶ��Ҳ��͡A�ҥH���ƩI�s�ɡA�|�������P���üƭ�(���)�C
***************************************************************************************************/
int RandomNumber(int MinValue, int MaxValue)   //Both MinValue and MaxValue are included
{
	int R = (rand() % (MaxValue - MinValue + 1)) + MinValue;
	return R;
}