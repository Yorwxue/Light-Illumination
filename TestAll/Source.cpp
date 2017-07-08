#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//�ϥΪ̩w�q
//Ū�ɨ��o
int NumofLight;       //�O���Ӽ�
int NumofInstruction; //�G�׫��O�Ӽ�
int NumofTd;          //�ؼ��I�Ӽ�


#define OtherLight 0        //�B�~����-�Ӷ�
//float OtherLightT[OtherLight]={0};        //�B�~�������G��




//�{���۩w�q
#define ImpossibleResume 99999

double ResumeTime;




float** I;   //(�ӯ�,�G��)
float* Td;
float** KParameter;

int CheckT(float T,int MatrixRow,int i,int j,int k);
void TestAllSollution();
void* new2d(int h, int w, int size);

int RandomNumber(int MinValue,int MaxValue);


/*�B�~���i��������B�z(����)*
void OtherLightProcess()
{
	for(int i=0;i<NumofTd;i++)
	{
		Td[i]-=KParameter[i][NumofLight+OtherLight-1];//*OtherLightT[OtherLight];
	}
}
/*******************/



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

		sscanf_s(FileInput, "%d", &NumofInstruction);

		I = (float **)new2d(NumofInstruction, 2, sizeof(float));

		for(int i=0;i<NumofInstruction;i++)
		{
			for(int j=0;j<2;j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&I[i][j]);        //�r����Ʀr
				//cout<<I[i][j]<<",";
			}
			//cout<<endl;
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

		sscanf_s(FileInput, "%d", &NumofTd);

		Td = new float[NumofTd];

		for(int i=0;i<NumofTd;i++)
		{
			fr.getline(FileInput,sizeof(FileInput),',');
			sscanf(FileInput,"%f",&Td[i]);        //�r����Ʀr
			//cout<<Td[i]<<",";
		}
		//cout<<endl;
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

		sscanf_s(FileInput, "%d", &NumofLight);

		KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));

		for(int i=0;i<NumofTd;i++)
		{
			for(int j=0;j<(NumofLight+OtherLight);j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&KParameter[i][j]);        //�r����Ʀr
				//cout<<KParameter[i][j]<<",";
			}
			//cout<<endl;
		}
	}
	fr.close();
	//���Ū������


	//OtherLightProcess();//�B�~�����B�z


	clock_t start, stop;
    start = clock(); //�}�l�ɶ�

	TestAllSollution();
	
	stop = clock(); //�����ɶ�

	cout<<endl;
	cout<<"�O��:"<<(double(stop-start)/CLOCKS_PER_SEC)<<"��"<<endl<<endl;
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
/********************************Test all solutions******************************************/
void TestAllSollution()
{
	int NowProcess=0;
	int* BestInstruction=new int[NumofLight];//(�Ĥ@���O�A�ĤG���O�A�ĤT���O...)
	for(int i=0;i<NumofLight;i++)//�]�w���i�઺�ȧ@����ȡA��K�᭱�P�_�ثe�S����
	{
		BestInstruction[i]=-1;
	}     
	//                                                                ��������
	//[�Ĥ@���O���ե����O][�ĤG���O���ե����O][�ĤT���O���ե����O]...   <---->   RangeofLightwithInstructions
	//EX: 3���O�A10�ث��O�A�@1000�زզX�覡
	//[2][8][9]   <--->   2*10^2 + 8*10^1 + 9*10^0                             //Note:�Ĥ@���O�b�̰���!!
	//EX:6���O�A15�ث��O�A�@11390615�زզX�覡
	//[4][2][3][1][5][8]   <--->   4*15^5 + 2*15^4 + 3*15^3 + 1*15^2 + 5*15^1 + 8*15^0
	//�H6���O�|�ҡA�ե����O�W��40�ءA�A�W�L�N�|�j��int���A�ү�������Ʀr
	float BestPower=0;//�ӯ�
	int T=0;
	                                                                                                             int Count=0,Count2=0;//�p�ƾ�
																												 fstream fw,fw2,fw3;
					                                                                                             string filename="Specific Solution_best.txt";
																												 string filename2="Specific Solution_Range_10.txt";
																												 string filename3="Now Process.txt";
																										         fw.open(filename, ios::out);
																												 fw2.open(filename2, ios::out);
																												 fw3.open(filename3, ios::out);
	clock_t start, stop;
	cout<<"����j�M"<<endl;
    start = clock(); //�}�l�ɶ�
	int* TestAllInstruction=new int[NumofLight];
	for(int s=0;s<NumofLight;s++)
	{
		TestAllInstruction[s]=0;
	}
	//
	for(int Li=0;Li<NumofLight;Li++)    //Li�O�O���s��
	{
		for(int InstructionLiterally=0;InstructionLiterally<=Li;InstructionLiterally++)  //�v���W�[���O
		{
	        for(int Ii=0;Ii<NumofInstruction;Ii++)    //Ii�O�O�����O��
	        {
				//��ܥثe�i��
				NowProcess++;
				if(NowProcess==1000000)             //���O�C�@������X
				{
					NowProcess=0;
				    cout<<"�ثe�i��:";                               
				    for(int s=0;s<NumofLight;s++)
	                {
		                cout<<TestAllInstruction[s]<<",";
	                }
				    cout<<endl;
				}
	    	//�ˬd�O�_�X��
	    	    for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	    {
		    	    T=CheckT(MatrixRow,TestAllInstruction);
					                                                                                            //�ثe�i��
					                                                                                            NowProcess++;
				                                                                                                if(NowProcess==1000000)             //���O�C�@������X
				                                                                                                {
					                                                                                                NowProcess=0;
				                                                                                                    fw3<<"�ثe�i��:";                               
				                                                                                                    for(int s=0;s<NumofLight;s++)
	                                                                                                                {
		                                                                                                                fw3<<TestAllInstruction[s]<<",";
	                                                                                                                }
				                                                                                                    fw3<<endl;
				                                                                                                }
		    	    if(T==0)//���X��
		        	{                                                                                           //Count++;//���X��p��
		    	    	break;
		    	    }
		    	    else//�X��
		    	    {
		    		    T=0;
		    	    	if(MatrixRow==NumofTd-1)//�ŦX����A�}�l�p��ӯ�
		    	    	{
		    		    	float Power=0;
		    		    	for(int k=0;k<NumofLight;k++)
			    	    	{
			    			    int j=TestAllInstruction[k];
			    			    Power+=I[j][0];
			    		    }
							                                                                            /*�[�ݥثe�X��ɪ����O*/
					                                                                                    if(((Power-BestPower)<=(BestPower*0.1)) || BestPower==0)
																										{
																											Count2++;
																											fw2<<"���O:";
																											for(int k=0;k<NumofLight;k++)
					                                                                                        {
																												fw2<<TestAllInstruction[k]<<",";
					                                                                                        }
																											fw2<<"�ӯ�:"<<Power;
																											fw2<<endl;
																										}
																										/***************************/
					                                                                                    
				        	if(Power<=BestPower||BestPower==0)//�{�b���G�ײŦX�A�B�̬ٹq
				        	{
				    	    	BestPower=Power;
				    	    	for(int k=0;k<NumofLight;k++)
				    		    {
				    		    	BestInstruction[k]=TestAllInstruction[k];
				    		    }
				    		    cout<<"�ثe�̨Ϋ��O:";
				    		    for(int k=0;k<NumofLight;k++)
					            {
									cout<<TestAllInstruction[k]<<",";
					            }
				    	    	cout<<endl;
								                                                                       /*�[�ݥثe�̨ίӯ�ɪ����O*/
					                                                                                    //if(Power==30.5)
																										//{
																											Count++;
																											fw<<"���O:";
																											for(int k=0;k<NumofLight;k++)
					                                                                                        {
																												fw<<TestAllInstruction[k]<<",";
					                                                                                        }
																											fw<<"�ӯ�:"<<BestPower;
																											fw<<endl;
																										//}
																										/***************************/
																										
																										/***************************/
					        }
				        }
			        }
		        }
			    //�X���ˬd�A�P�����AEND

		    	TestAllInstruction[0]++;
				int test=0;
				for(int j=0;j<NumofLight;j++)
    			{
					if(TestAllInstruction[j]==NumofInstruction-1)
    						test++;
    			}
    			if(test==(NumofLight))
    				break;
	        }
		//
			for(int CheckLiterally=0;CheckLiterally<=InstructionLiterally;CheckLiterally++)
			{
				if( TestAllInstruction[CheckLiterally]==NumofInstruction && (CheckLiterally+1)!=(NumofLight))
		        {
			        TestAllInstruction[CheckLiterally]=0;
				    TestAllInstruction[CheckLiterally+1]++;
		        }
			}
			if(InstructionLiterally+1<NumofLight)
    			{
    				if(TestAllInstruction[InstructionLiterally+1]==0)
    				    InstructionLiterally--;
    			}
    			
			if((InstructionLiterally+1)==NumofLight || TestAllInstruction[NumofLight-1]==NumofInstruction-1)
    			{
    				int EndorNot=0;
    				for(int i=0;i<NumofLight;i++)
    				{
    					EndorNot+=TestAllInstruction[i];
    				}
					if( EndorNot <  NumofLight*(NumofInstruction-1)  )
    				{
    					InstructionLiterally--;
    				}
    				else
    				{
    					Li=NumofLight;
    					break;
    				}
    			}
    		}
    		//
    	}
	    stop = clock(); //�����ɶ�
	                                                                                                      fw<<"�ƥ�:"<<Count<<endl;
																										  fw2<<"�ƥ�:"<<Count2<<endl;
                                                                                                          fw<<"�O��:"<<double(stop-start)/CLOCKS_PER_SEC<<"��"<<endl<<endl;
	                                                                                                      fw.close();
																										  fw2.close();
																										  fw3.close();
	
	cout<<"���O:";
	for(int i=0;i<NumofLight;i++)
	{
		cout<<BestInstruction[i]<<",";
	}
	cout<<endl;
	cout<<"�ӯ�:"<<BestPower<<endl;
    cout<<"�O��:"<<double(stop-start)/CLOCKS_PER_SEC<<"��"<<endl<<endl;
}
/*************************************Test all end*********************************************/
/**********************************************************************************************/

/****************************************Random Number**********************************************
�����Ʊ檺�̤j�ȻP�̤p�ȡA�N�ಣ�ͤ���o��ӭȤ����H���ơA�o���H���üƤ]�i�൥��̤j�ȡA�ε���̤p�ȡC
��üƭȬO�̷Ӯɶ��Ҳ��͡A�ҥH���ƩI�s�ɡA�|�������P���üƭ�(���)�C
***************************************************************************************************/
int RandomNumber(int MinValue,int MaxValue)   //Both MinValue and MaxValue are included
{
	int R=(rand()%(MaxValue-MinValue+1))+MinValue;
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