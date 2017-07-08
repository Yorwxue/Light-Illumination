#include<stdio.h>
#include<sys/types.h>
#include<sys/socket.h>
#include<netinet/in.h>
#include<arpa/inet.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<unistd.h>
#include <pthread.h>
#include<string>

#define NUM_THREADS 2
pthread_t threads[NUM_THREADS];

void* UserThread(void* Input);
void* ProxyThread(void* Input);
void* CommunicationThread(void* Input);

char SendBuffer[1024];
char ReceiveBuffer[1024];
int HistoryMessage = 0;
int ExitMessage = 0;

int main(void)
{
	printf("Chatting Room.\n");
	printf("Enter \"History\" to see the message before.\n");
	printf("Enter \"Exit\" to leave the Chatting Room.\n");

	pthread_create(&threads[0], NULL, UserThread, NULL);
	pthread_create(&threads[1], NULL, ProxyThread, NULL);
	pthread_create(&threads[2], NULL, CommunicationThread, NULL);
}


void* UserThread(void* Input)
{
	while (1)
	{
		char message[1024];
		scanf("%s", &message);
		if (message == "Exit")
		{
			ExitMessage = 1;
			break;
		}
		else if (message == "History")
		{
			HistoryMessage = 1;
		}
		else
		{
			strncpy(SendBuffer, message, sizeof(message));
		}
	}
	pthread_exit(NULL);
}

void* ProxyThread(void* Input)
{
	while (ExitMessage != 1)
	{
		if (ReceiveBuffer != NULL)
		{
			printf("%s\n", ReceiveBuffer);
			strcpy(ReceiveBuffer, "");
		}
	}
	pthread_exit(NULL);
}

void* CommunicationThread(void* Input)
{
	int sockfd, length;
	char serverIP[15];
	struct sockaddr_in myaddr;
	struct sockaddr_in servaddr;

	/*=================================================
	htonl():

	功能說明：將32位元主機字元順序轉換成網路字元順序。

	標頭檔：#include <netinet/in.h>

	函式宣告：unsigned long int htonl(unsigned long int hostlong);

	函式說明：htonl()用來將引數指定的32位元hostlong轉換成網路字元順序。

	傳回值：傳回對應的網路字元順序。

	htons()

	功能說明：將16位元主機字元順序轉換成網路字元順序。

	標頭檔：#include <netinet/in.h>

	函式宣告：unsigned short int htons(unsigned short int hostshort);

	函式說明：htons()用來將引數指定的16位元hostshort轉換成網路字元順序。

	傳回值：傳回對應的網路字元順序。

	=================================================*/
	if ((sockfd = socket(AF_INET, SOCK_STREAM, 0))<0){
		printf("socket error\n");
		return 0;
	}
	bzero(&myaddr, sizeof(myaddr));
	myaddr.sin_family = AF_INET;
	myaddr.sin_addr.s_addr = htonl(INADDR_ANY);
	myaddr.sin_port = htons(0);
	if (bind(sockfd, (struct sockaddr*)&myaddr, sizeof(myaddr))<0){
		printf("bind error\n");
		return 0;
	}

	bzero(&servaddr, sizeof(servaddr));
	printf("Enter Server IP : ");
	scanf("%s", serverIP);
	servaddr.sin_family = AF_INET;
	servaddr.sin_port = htons(1234);
	servaddr.sin_addr.s_addr = inet_addr(serverIP);

	if (connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))<0){
		printf("connect failed!");
		return 0;
	}
	else
		printf("connect success\n");

	while (ExitMessage != 1)
	{
		if (HistoryMessage == 1)
		{
		}
		else
		{
			//傳送
			if (SendBuffer != NULL)
			{
				write(sockfd, SendBuffer, sizeof(SendBuffer));
				strcpy(SendBuffer, "");
			}
			//接收
			read(sockfd, ReceiveBuffer, sizeof(ReceiveBuffer));
		}
	}

	close(sockfd);
	pthread_exit(NULL);
}