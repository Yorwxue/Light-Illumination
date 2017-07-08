#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//�ϥΪ̩w�q
#define NumofLight 12       //�O���Ӽ�
#define NumofInstruction 16 //�G�׫��O�Ӽ�
#define NumofTd 5          //�ؼ��I�Ӽ�
#define NumofCPU 10

int RunTime=10;  //�]�h�֦�
float RealBest=246.34;
float totaltime=0;
//float LongestTime=0;
//float ShortestTime=999999;

#define OtherLight 0        //�B�~����-�Ӷ�
//float OtherLightT[OtherLight]={0};        //�B�~�������G��


#define NumofParticle 1000        //��NumofParticle�ոѧ@��PSO�����(�ɤl��)
#define NumofRedistribute 100     //10%

//#define EndFlagRepeatMAXTimes 2
//int EndFlagRepeatTimes=0;
float RepeatValue=-1;
float DifferentRange=0.0005;

float ave_power = 0;

//�{���۩w�q
#define ImpossibleResume 99999
#define InitialW 0.9             //�D�ʫY�ơA�Ѥj��p
#define FinalW 0.4
#define InitialCOne 2.5          //������C�Y�ơA�Ѥj��p�A�]���j�M�����������
#define FinalCOne 1
#define InitialCTwo  1       //�s�骺C�Y�ơA�Ѥp��j�A�]���j�M���������s��
#define FinalCTwo 2.5
int CircleofPSO;
int VelocityLimit;
int NumofCPUCircle;
double ResumeTime;

void Parameter()
{
    CircleofPSO = 50;//NumofParticle*(0.02);      //PSO���檺�`�����ơA�ɤl�ƪ�2%
	//if(CircleofPSO<10)   //�������ƹL�ַ|�y���L��
	//	CircleofPSO=10;

    VelocityLimit=5;//NumofInstruction*(0.5);  //���ʳt�פW���A���ʪŶ���50%
	//if(VelocityLimit<3)   //���ʭ�����Y�A�|�ɭP���鰱���
	//	VelocityLimit=3;  

    NumofCPUCircle=10;//NumofCPU;    //�U�B�z������y�һݡA�����B�z���Ӽ�
}
//



float I[NumofInstruction][2];   //(�ӯ�,�G��)
float Td[NumofTd];
float KParameter[NumofTd][NumofLight+OtherLight];

float start_point_range[NumofLight][2];    //(�C�ӫ׫��O,���ӫ׫��O)

int CheckT(float T,int MatrixRow,int i,int j,int k);
void TestAllSollution();
float PSO();
int RandomNumber(int MinValue,int MaxValue);
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
	float PBest[NumofLight+1],  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,����̧C�ӯ�)
		  GBest[NumofLight+1];  //(�Ĥ@���O���O,�ĤG���O���O,�ĤT���O���O,�s��̧C�ӯ�)
	int Position[NumofLight];  //�O������A��U���O�����ե����O <----> Position   //�䤤�O���Ƥ]�Osolution�����סA�ĴX���O��ܲĴX��
	int Velocity[NumofLight];             //�ե����O�W�ɩΤU�� <----> Velocity
	float W;//�����D�ʥΪ�weighting�A��ȶV�j�A�V�A�X����j�M�A��ĳ��0.9~0.4
	float C1,C2;//�Y�ơA�վ㰾�V�U��A�ΰ��V�s��
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
		for(int i=0;i<NumofLight;i++)//i�O����
		{                                                                                               //cout<<"��P"<<i<<":"<<Position[i]<<",";
		Position[i]=Position[i]+Velocity[i];                                                            //cout<<"V:"<<Velocity[i]<<endl;
			if(Position[i]>NumofInstruction-1)             //�ե����O���W�U��
				Position[i]=NumofInstruction-1;
			if(Position[i]<0)
				Position[i]=0;                                                                           //cout<<"��P"<<i<<":"<<Position[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void UpdateVelocity()
	{
		for(int i=0;i<NumofLight;i++)//i�O����
		{                                                                                               //cout<<"��V"<<i<<":"<<Velocity[i]<<",";
			Velocity[i] = (W*Velocity[i] + C1*((float)RandomNumber(0, 100) / 100)*(PBest[i] - Position[i]) + C2*((float)RandomNumber(0, 100) / 100)*(GBest[i] - Position[i]));         //cout<<"��V"<<i<<":"<<Velocity[i]<<endl;
			if(Velocity[i]>VelocityLimit)                   //�t�פW�U��
				Velocity[i]=VelocityLimit;
			if(Velocity[i]<(-VelocityLimit))
				Velocity[i]=(-VelocityLimit);                                                            //cout<<"��V"<<i<<":"<<Velocity[i]<<endl;
		}                                                                                                //cout<<endl;system("pause");
	}
	void CheckPosition()
	{
		float Tpso=0;  //�ھڥثe�ɤl����m�Һ�X�Ӫ�T(�ӫ�)��
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
//PSO���O����



int main()
{
	char FileInput[50];
	fstream fr;
	srand((unsigned)time(NULL));
	//Light Instruction
	fr.open("LightInstruction.txt",ios::in);
	if(!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//�o�ONumofInstruction
		for(int i=0;i<NumofInstruction;i++)
		{
			for(int j=0;j<2;j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&I[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();

	//Td
	fr.open("Td.txt",ios::in);
	if(!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//�o�ONumofTd
		for(int i=0;i<NumofTd;i++)
		{
			fr.getline(FileInput,sizeof(FileInput),',');
			sscanf(FileInput,"%f",&Td[i]);        //�r����Ʀr
		}
	}
	fr.close();
	
	//KParameter
	fr.open("KParameter.txt",ios::in);
	if(!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//�x�}��
		fr.getline(FileInput,sizeof(FileInput),',');//�x�}�e
		for(int i=0;i<NumofTd;i++)
		{
			for(int j=0;j<(NumofLight+OtherLight);j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&KParameter[i][j]);        //�r����Ʀr
			}
		}
	}
	fr.close();

	//start range
	fr.open("start_point.txt", ios::in);
	if (!fr)        //�p�G�}���ɮץ��ѡAfin��0�F���\�Afin��1
		cout << "Fail to open file" << endl;
	else
	{
		for (int i = 0; i<NumofLight; i++)
		{
				fr.getline(FileInput, sizeof(FileInput), ',');
				float prob_solu;
				sscanf(FileInput, "%f", &prob_solu); //�r����Ʀr
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

	//���Ū������

	Parameter();

	//OtherLightProcess();//�B�~�����B�z
	/*******PSO���������G**********/
	clock_t start, stop;
    start = clock(); //�}�l�ɶ�
	float AverageResult=0;
	float TestAnswer;//�O�_�L�ѩΫD�̨θ�
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
	int NoAnswer=0;//�L���`����
	
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
	stop = clock(); //�����ɶ�
	
	cout<<endl;
	cout<<"����"<<RunTime<<"��"<<endl;
	cout<<"�̨θѬ��G"<<RealBest<<endl;
	//cout<<"�����ӯ�:"<<(AverageResult/(RunTime-NoAnswer))<<endl;
	cout << "�����ӯ�" << ave_power / RunTime << endl;
	cout << "�����O��:" << totaltime / RunTime << "��" << endl;/**
		<<"�̪��O��:"<<LongestTime<<"��"<<endl<<endl
		<<"�̵u�O��:"<<ShortestTime<<"��"<<endl<<endl;/**/
	//cout<<"�����O��:"<<ResumeTime<<"��"<<endl<<endl;
	//cout<<"�L�Ѧ���:"<<NoAnswer<<"��"<<endl;
	cout << "�̨θѭӼ�:" << BestTimes << endl/*
	    <<"0.1%���Ӽ�:"<<Range01Times<<endl
		<<"0.2%���Ӽ�:"<<Range02Times<<endl
		<<"0.3%���Ӽ�:"<<Range03Times<<endl
		<<"0.4%���Ӽ�:"<<Range04Times<<endl
		<<"0.5%���Ӽ�:"<<Range05Times<<endl
		<<"1%���Ӽ�:"<<Range1Times<<endl
		<<"2%���Ӽ�:"<<Range2Times<<endl
		<<"3%���Ӽ�:"<<Range3Times<<endl
		<<"4%���Ӽ�:"<<Range4Times<<endl
		<<"5%���Ӽ�:"<<Range5Times<<endl
		<<"Other:"<<RangeOther<<endl/**/
		<<"�L�ѭӼơG"<<NoAnswer<<endl;
	//cout<<"10%���Ӽ�"<<Range10Times<<endl;
	/***********************************************/
	system("pause");
	return 0;
}

//�T�{�G�צX��P�_
int CheckT(int MatrixRow,int i[])     //MatrixRow���n�ˬd�O�_�X�檺row�Ai���������O���x�}
{
	float T=0;
	for(int k=0;k<NumofLight;k++)
	{
	        T+=KParameter[MatrixRow][k]*I[i[k]][1];//[��ܿO�����O,1��ܫG��]
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
	int ParticleNo=0;                  //�ؼЬO���NumofParticle�Ӹѧ@����l�ȡA�ҥH�]�@�ӺX�ХΥH�p��
	Particle* Group[NumofCPU][NumofParticle];
	float T=0;

	/*********���������t�ɨϥ�************/
	Particle* TempGroup[NumofCPU][NumofParticle];
	int TempGroupNo[NumofCPU];  //�p��洫�ɤl�Ϊ��Ȧs�s�O�_�w��
	for(int i=0;i<NumofCPU;i++)
	{
		TempGroupNo[i]=0;
	}
	/**************************************/

    //Distribute_start = clock(); //�}�l�ɶ�

    float StartPointRatio[NumofLight];
	int PSOStartCenter[NumofLight];

	PSOstart = clock();
/********************************�H�����ͪ�l�ȡA���ר�O�_�X��*****************************************/
	int RandomVelocity[NumofLight];
	for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
	{
	    ParticleNo=0;
		for(int i=0;i<NumofParticle;i++)                                       //����j�������зj�M"��������"�Aj�s�F�Ҧ��O�����O��T
	    {
			int j[NumofLight];
			//�ɤl��l��m�H�����ͤ覡
			/**/
		    for(int h=0;h<NumofLight;h++)                                       
	        {                                                                                                                                                    //cout<<StartPointCheck<<",";
		    	j[h]=RandomNumber(0,NumofInstruction-1);   
	    	}/**/                                                                                          //cout<<endl;
			//�S�w��m�g�䲣��
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
	    	//��l��m���ͧ���

	    	for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	{
		        T=CheckT(MatrixRow,j);
	    	    if(T==0)                                                     //���X��APower�]ImpossibleResume�A���i��o�ͪ����p�A�@�����
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
		        if(MatrixRow==NumofTd-1)//����Row���ŦX����A�}�l�p��ӯ�
		        {                                                                             //cout<<"�e"<<Power<<endl;system("pause");
		    	    for(int k=0;k<NumofLight;k++)                                                 //�p��ӯ�
		    	    { 
			    	    Power+=I[j[k]][0];                                                       //cout<<Power<<endl;system("pause");
			        }                                                                         //cout<<"��"<<Power<<endl;system("pause");
		    		for(int h=0;h<NumofLight;h++)                                                 //�H����l�t��
	                {
		                RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                }                                                                                     //cout<<Power<<endl;system("pause");
                    Group[CPUNo][ParticleNo++]=new Particle(j,Power,&RandomVelocity[0]);                //���Ͳɤl
		    		Power=0;
		        }
		    }
	    }
	    for(int i=0;i<NumofParticle;i++)
	    {
	        PBestPoint=Group[CPUNo][i]->GetPBest();                                //cout<<*(PBestPoint+NumofLight)<<endl;system("pause");
	    	if(*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
	    	{
		    	for(int j=0;j<NumofLight+1;j++)
	            {
		    		GBest[CPUNo][j]=*(PBestPoint+j);                               //cout<<GBest[CPUNo][j]<<",";system("pause");
	            }                                                           //cout<<GBest[CPUNo][NumofLight]<<",";system("pause");//cout<<endl;system("pause");
		    }
	    }
	}
	                                                                                             //cout<<"��1�����X��ɤl��:"
	//Distribute_stop = clock();//�����ɶ�
	//ResumeTime=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
	//��Ȳ��ͧ���
	//�P�ɡA�Ĥ@��PSO�����A�ѤUCircleofPSO-1��
	//�s�W����B�z�A�b�̥~�h
/*******************************************************************************************************/
	for(int CPUCircle=0;CPUCircle<NumofCPUCircle;CPUCircle++)
	{
	    for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
	    {
			//if(CPUNo==0)
			//{
			//	 Distribute_start = clock(); //�}�l�ɶ�
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
			        if((*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)&&*(PBestPoint+NumofLight)!=ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
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
			//	Distribute_stop = clock();//�����ɶ�
		    //    ResumeTime+=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
			//}
	    }

		//Distribute_start = clock(); //�}�l�ɶ�
		
		

		//�����覡
		/****************************************************
		//��ҥ��]����
		int RedistributeList[NumofCPU];
		int RedistributeParticleList[NumofCPU][NumofRedistribute];
		int RandomGroup;
		int RandomParticle;

		//�إ߸s��list
		for(int i=0;i<NumofCPU;i++)  
		{
			RandomGroup=RandomNumber(0,NumofCPU-1);
			for(int j=0;j<i;j++)
			{
			    if(RedistributeList[j]==RandomGroup)
			    {
			        RandomGroup=RandomNumber(0,NumofCPU-1);
				    j=-1;                                                             //for loop���Y�A�ˬd�e�|�[1�A�ҥH�]-1�A���L�q0�}�l�ˬd
				}
			}
			RedistributeList[i]=RandomGroup;
		}

		//��ܭn�洫���ɤl�A�إ߹�Ӫ�
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
				        j=-1;                                                //for loop���Y�A�ˬd�e�|�[1�A�ҥH�]-1�A���L�q0�}�l�ˬd
				    }
			    }
			    RedistributeParticleList[k][i]=RandomParticle;
			}
		}
		
		//�}�l�洫
		for(int i=0;i<NumofCPU;i++)  //�ƥ�
		{
			for(int j=0;j<NumofParticle;j++)
			{
				TempGroup[i][j]=Group[i][j];
			}
		}
		for(int i=0;i<NumofCPU;i++) //i�O�Ҧb���s,k�O�ĴX���@�s�������洫
		{
			for(int j=0;j<NumofRedistribute;j++)
		    {
				Group[ RedistributeList[i] ][ RedistributeParticleList[ RedistributeList[i] ][j] ] = TempGroup[i][ RedistributeParticleList[i][j] ];
			}
		}
		/************��������************
		int RandomGroup;
		for(int i=0;i<NumofCPU;i++)
		{
			for(int j=0;j<NumofParticle;j++)
			{
				RandomGroup=RandomNumber(0,NumofCPU-1);
			    while(TempGroupNo[RandomGroup]==NumofParticle)//�Ȧs�s�w��
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
		//////////////�����������s�A���s�զX
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
		///////////*///��������

		//���ͷs�ɤl�A�Ӥ�����
		/****************************************************
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)  
		{
			//�ت�
			int RedistributeParticleList[NumofRedistribute];
			for(int i=0;i<NumofRedistribute;i++)
		    {
			    int RandomReplaced=RandomNumber(0,NumofParticle-1);
			    for(int j=0;j<i;j++)
			    {
			        if(RedistributeParticleList[j]==RandomReplaced)
			        {
			            RandomReplaced=RandomNumber(0,NumofParticle-1);
				        j=-1;                                                //for loop���Y�A�ˬd�e�|�[1�A�ҥH�]-1�A���L�q0�}�l�ˬd
				    }
			    }
			    RedistributeParticleList[i]=RandomReplaced;
			}
			//ñ�槹���A�}�l���s����
		    for(int i=0;i<NumofRedistribute;i++)                                       //����j�������зj�M"��������"�Aj�s�F�Ҧ��O�����O��T
	        {
		    	int j[NumofLight];
		    	//�ɤl��l��m�H�����ͤ覡
		        for(int h=0;h<NumofLight;h++)                                       
	            {                                                                                                                                                    //cout<<StartPointCheck<<",";
		        	j[h]=RandomNumber(0,NumofInstruction-1);   
	        	}                                                                                     
			
			    //
	        	//��l��m���ͧ���

	        	for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	    {
		            T=CheckT(MatrixRow,j);
	    	        if(T==0)                                                     //���X��APower�]ImpossibleResume�A���i��o�ͪ����p�A�@�����
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
		            if(MatrixRow==NumofTd-1)//����Row���ŦX����A�}�l�p��ӯ�
		            {                                                                        
		    	        for(int k=0;k<NumofLight;k++)                                                 //�p��ӯ�
		    	        { 
			    	        Power+=I[j[k]][0];                                                  
			            }  

						if(Power<TotalGBest[NumofLight])
						{
							TotalGBest[NumofLight]=Power;
						}

		    		    for(int h=0;h<NumofLight;h++)                                                 //�H����l�t��
	                    {
		                    RandomVelocity[h]=(RandomNumber(0,1)?1:-1)*RandomNumber(0,VelocityLimit);
	                    }                                                                               
						delete Group[CPUNo][RedistributeParticleList[i]];
                        Group[CPUNo][RedistributeParticleList[i]]=new Particle(j,Power,&RandomVelocity[0]);                //���Ͳɤl
		    	    	Power=0;
		            }
				}
			}
	    }
		/*************************���ͧ���***************************/

		/*�M�䭫���β��᪺ͫ�U�s��gbest*
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
		{
			for(int i=0;i<NumofParticle;i++)
			{
				GBest[CPUNo][NumofLight]=ImpossibleResume;       
			}
			for(int i=0;i<NumofParticle;i++)
			{
			        PBestPoint=Group[CPUNo][i]->GetPBest();
			        if((*(PBestPoint+NumofLight)< GBest[CPUNo][NumofLight]||GBest[CPUNo][NumofLight]==ImpossibleResume)&&*(PBestPoint+NumofLight)!=ImpossibleResume)    // GBest[NumofLight]�O�ӯ�
			        {
			            for(int j=0;j<NumofLight+1;j++)
	                    {
	        	            GBest[CPUNo][j]=*(PBestPoint+j);
				       }
			        }
			} 
		}
		/****�M�䧹��****/
		for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)
		{
			/*�U�s���s��X��gbest�g�^*/
			for(int i=0;i<NumofParticle;i++)
			{
			    Group[CPUNo][i]->SetGBest(&GBest[CPUNo][0]);
			}
			/******/
			/*�N��X��Gbest�g�^���U�s*
			for(int i=0;i<NumofLight+1;i++)
			{
			    GBest[CPUNo][i]=TotalGBest[i];
			}
			/*********/
		}
		//Distribute_stop = clock();//�����ɶ�
	    //ResumeTime+=double(Distribute_stop-Distribute_start)/CLOCKS_PER_SEC;
		/********/
		/***�{���D�X���Ѥwí�w***
		if (RepeatValue == -1 && TotalGBest[NumofLight] != ImpossibleResume)
		{
			RepeatValue=TotalGBest[NumofLight];
		}
		else if(RepeatValue!=-1)
		{
			 //�z�פWTotalGBest[NumofLight]�@�w�p��ε���RepeatValue�A�]��RepeatValue�N�O�W�@����TotalGBest[NumofLight]�A�o�ӧP�_���O�O�I�@��
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
	
	/***PSO���浲���A�[�ݵ��G***/
	if(TotalGBest[NumofLight]!=ImpossibleResume)
	{
	    //cout<<"PSO�j�M"<<endl
		cout<<"���O:";                                                                            
	    for(int i=0;i<NumofLight;i++)
	    {
		    cout<<TotalGBest[i]<<",";  
	    }
	    cout<<"�ӯ�:"<<TotalGBest[NumofLight];                                           
	}
	else
	{
		cout<<"�L��"<<endl;
	}

	/****�g��****/
	if(TotalGBest[NumofLight]<=RealBest || RealBest==-1)
	{
	    fstream fw;
	    string filename="The_Best_Solution.txt";
	    fw.open(filename, ios::out);
	    fw<<"PSO�j�M"<<endl<<"���O:";
	    for(int i=0;i<NumofLight;i++)
	    {
		    fw<<TotalGBest[i]<<",";  
	    }
	    fw<<endl<<"�ӯ�:"<<TotalGBest[NumofLight]<<endl;    
		fw.close();
	}
	/****�g�ɧ���****/

	PSOstop = clock(); //�����ɶ�
	float ExecutingTime=float(PSOstop-PSOstart)/CLOCKS_PER_SEC;
    cout<<"�O��:"<<ExecutingTime<<"��"<<endl;
	totaltime+=ExecutingTime; 
	/**
	if(ExecutingTime>LongestTime)
		LongestTime=ExecutingTime;
	if(ExecutingTime<ShortestTime)
		ShortestTime=ExecutingTime;/**/
    /*****************************/
	for(int CPUNo=0;CPUNo<NumofCPU;CPUNo++)  //����Ŷ�
	 {
		 for(int i=0;i<NumofParticle;i++)
		 {
			 delete Group[CPUNo][i];
		 }
	 } /**/

	return TotalGBest[NumofLight];
}

/****************************************Random Number**********************************************
�����Ʊ檺�̤j�ȻP�̤p�ȡA�N�ಣ�ͤ���o��ӭȤ����H���ơA�o���H���üƤ]�i�൥��̤j�ȡA�ε���̤p�ȡC
��üƭȬO�̷Ӯɶ��Ҳ��͡A�ҥH���ƩI�s�ɡA�|�������P���üƭ�(���)�C
***************************************************************************************************/
int RandomNumber(int MinValue,int MaxValue)   //Both MinValue and MaxValue are included
{
	int R=(rand()%(MaxValue-MinValue+1))+MinValue;
	return R;
}
