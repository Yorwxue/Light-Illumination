#include<iostream>
#include<stdlib.h>
using namespace std;

void* new2d(int h, int w, int size);

void main()
{
	int height = 3, width = 2;
	float** I = (float **)new2d(height, width, sizeof(float));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			I[i][j] = i*width + j;
		}
	}

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			cout << I[i][j]<<",";
		}
		cout << endl;
	}
	delete[] I;
	system("pause");
}

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