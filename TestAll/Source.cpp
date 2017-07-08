#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include <ctime>
using namespace std;

//使用者定義
//讀檔取得
int NumofLight;       //燈的個數
int NumofInstruction; //亮度指令個數
int NumofTd;          //目標點個數


#define OtherLight 0        //額外光源-太陽
//float OtherLightT[OtherLight]={0};        //額外光源的亮度




//程式自定義
#define ImpossibleResume 99999

double ResumeTime;




float** I;   //(耗能,亮度)
float* Td;
float** KParameter;

int CheckT(float T,int MatrixRow,int i,int j,int k);
void TestAllSollution();
void* new2d(int h, int w, int size);

int RandomNumber(int MinValue,int MaxValue);


/*額外不可控制光源處理(陽光)*
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

	if(!fr)        //如果開啟檔案失敗，fin為0；成功，fin為1
        cout<<"Fail to open file"<<endl;
	else
	{
	    fr.getline(FileInput,sizeof(FileInput),',');//這是NumofInstruction

		sscanf_s(FileInput, "%d", &NumofInstruction);

		I = (float **)new2d(NumofInstruction, 2, sizeof(float));

		for(int i=0;i<NumofInstruction;i++)
		{
			for(int j=0;j<2;j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&I[i][j]);        //字串轉數字
				//cout<<I[i][j]<<",";
			}
			//cout<<endl;
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

		sscanf_s(FileInput, "%d", &NumofTd);

		Td = new float[NumofTd];

		for(int i=0;i<NumofTd;i++)
		{
			fr.getline(FileInput,sizeof(FileInput),',');
			sscanf(FileInput,"%f",&Td[i]);        //字串轉數字
			//cout<<Td[i]<<",";
		}
		//cout<<endl;
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

		sscanf_s(FileInput, "%d", &NumofLight);

		KParameter = (float **)new2d(NumofTd, NumofLight + OtherLight, sizeof(float));

		for(int i=0;i<NumofTd;i++)
		{
			for(int j=0;j<(NumofLight+OtherLight);j++)
			{
				fr.getline(FileInput,sizeof(FileInput),',');
				sscanf(FileInput,"%f",&KParameter[i][j]);        //字串轉數字
				//cout<<KParameter[i][j]<<",";
			}
			//cout<<endl;
		}
	}
	fr.close();
	//資料讀取完畢


	//OtherLightProcess();//額外光源處理


	clock_t start, stop;
    start = clock(); //開始時間

	TestAllSollution();
	
	stop = clock(); //結束時間

	cout<<endl;
	cout<<"費時:"<<(double(stop-start)/CLOCKS_PER_SEC)<<"秒"<<endl<<endl;
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
/********************************Test all solutions******************************************/
void TestAllSollution()
{
	int NowProcess=0;
	int* BestInstruction=new int[NumofLight];//(第一盞燈，第二盞燈，第三盞燈...)
	for(int i=0;i<NumofLight;i++)//設定不可能的值作為初值，方便後面判斷目前沒有解
	{
		BestInstruction[i]=-1;
	}     
	//                                                                高維轉單維
	//[第一盞燈的調光指令][第二盞燈的調光指令][第三盞燈的調光指令]...   <---->   RangeofLightwithInstructions
	//EX: 3盞燈，10種指令，共1000種組合方式
	//[2][8][9]   <--->   2*10^2 + 8*10^1 + 9*10^0                             //Note:第一盞燈在最高位!!
	//EX:6盞燈，15種指令，共11390615種組合方式
	//[4][2][3][1][5][8]   <--->   4*15^5 + 2*15^4 + 3*15^3 + 1*15^2 + 5*15^1 + 8*15^0
	//以6盞燈舉例，調光指令上限40種，再超過就會大於int型態所能紀錄的數字
	float BestPower=0;//耗能
	int T=0;
	                                                                                                             int Count=0,Count2=0;//計數器
																												 fstream fw,fw2,fw3;
					                                                                                             string filename="Specific Solution_best.txt";
																												 string filename2="Specific Solution_Range_10.txt";
																												 string filename3="Now Process.txt";
																										         fw.open(filename, ios::out);
																												 fw2.open(filename2, ios::out);
																												 fw3.open(filename3, ios::out);
	clock_t start, stop;
	cout<<"全域搜尋"<<endl;
    start = clock(); //開始時間
	int* TestAllInstruction=new int[NumofLight];
	for(int s=0;s<NumofLight;s++)
	{
		TestAllInstruction[s]=0;
	}
	//
	for(int Li=0;Li<NumofLight;Li++)    //Li是燈的編號
	{
		for(int InstructionLiterally=0;InstructionLiterally<=Li;InstructionLiterally++)  //逐條增加指令
		{
	        for(int Ii=0;Ii<NumofInstruction;Ii++)    //Ii是燈的指令數
	        {
				//顯示目前進度
				NowProcess++;
				if(NowProcess==1000000)             //不是每一筆都輸出
				{
					NowProcess=0;
				    cout<<"目前進度:";                               
				    for(int s=0;s<NumofLight;s++)
	                {
		                cout<<TestAllInstruction[s]<<",";
	                }
				    cout<<endl;
				}
	    	//檢查是否合格
	    	    for(int MatrixRow=0;MatrixRow<NumofTd;MatrixRow++)
	    	    {
		    	    T=CheckT(MatrixRow,TestAllInstruction);
					                                                                                            //目前進度
					                                                                                            NowProcess++;
				                                                                                                if(NowProcess==1000000)             //不是每一筆都輸出
				                                                                                                {
					                                                                                                NowProcess=0;
				                                                                                                    fw3<<"目前進度:";                               
				                                                                                                    for(int s=0;s<NumofLight;s++)
	                                                                                                                {
		                                                                                                                fw3<<TestAllInstruction[s]<<",";
	                                                                                                                }
				                                                                                                    fw3<<endl;
				                                                                                                }
		    	    if(T==0)//不合格
		        	{                                                                                           //Count++;//不合格計數
		    	    	break;
		    	    }
		    	    else//合格
		    	    {
		    		    T=0;
		    	    	if(MatrixRow==NumofTd-1)//符合條件，開始計算耗能
		    	    	{
		    		    	float Power=0;
		    		    	for(int k=0;k<NumofLight;k++)
			    	    	{
			    			    int j=TestAllInstruction[k];
			    			    Power+=I[j][0];
			    		    }
							                                                                            /*觀看目前合格時的指令*/
					                                                                                    if(((Power-BestPower)<=(BestPower*0.1)) || BestPower==0)
																										{
																											Count2++;
																											fw2<<"指令:";
																											for(int k=0;k<NumofLight;k++)
					                                                                                        {
																												fw2<<TestAllInstruction[k]<<",";
					                                                                                        }
																											fw2<<"耗能:"<<Power;
																											fw2<<endl;
																										}
																										/***************************/
					                                                                                    
				        	if(Power<=BestPower||BestPower==0)//現在的亮度符合，且最省電
				        	{
				    	    	BestPower=Power;
				    	    	for(int k=0;k<NumofLight;k++)
				    		    {
				    		    	BestInstruction[k]=TestAllInstruction[k];
				    		    }
				    		    cout<<"目前最佳指令:";
				    		    for(int k=0;k<NumofLight;k++)
					            {
									cout<<TestAllInstruction[k]<<",";
					            }
				    	    	cout<<endl;
								                                                                       /*觀看目前最佳耗能時的指令*/
					                                                                                    //if(Power==30.5)
																										//{
																											Count++;
																											fw<<"指令:";
																											for(int k=0;k<NumofLight;k++)
					                                                                                        {
																												fw<<TestAllInstruction[k]<<",";
					                                                                                        }
																											fw<<"耗能:"<<BestPower;
																											fw<<endl;
																										//}
																										/***************************/
																										
																										/***************************/
					        }
				        }
			        }
		        }
			    //合格檢查，與紀錄，END

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
	    stop = clock(); //結束時間
	                                                                                                      fw<<"數目:"<<Count<<endl;
																										  fw2<<"數目:"<<Count2<<endl;
                                                                                                          fw<<"費時:"<<double(stop-start)/CLOCKS_PER_SEC<<"秒"<<endl<<endl;
	                                                                                                      fw.close();
																										  fw2.close();
																										  fw3.close();
	
	cout<<"指令:";
	for(int i=0;i<NumofLight;i++)
	{
		cout<<BestInstruction[i]<<",";
	}
	cout<<endl;
	cout<<"耗能:"<<BestPower<<endl;
    cout<<"費時:"<<double(stop-start)/CLOCKS_PER_SEC<<"秒"<<endl<<endl;
}
/*************************************Test all end*********************************************/
/**********************************************************************************************/

/****************************************Random Number**********************************************
給予希望的最大值與最小值，就能產生介於這兩個值中的隨機數，這個隨機亂數也可能等於最大值，或等於最小值。
其亂數值是依照時間所產生，所以重複呼叫時，會給予不同的亂數值(整數)。
***************************************************************************************************/
int RandomNumber(int MinValue,int MaxValue)   //Both MinValue and MaxValue are included
{
	int R=(rand()%(MaxValue-MinValue+1))+MinValue;
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