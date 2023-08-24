#define  _CRT_SECURE_NO_WARNINGS  1
#pragma  warning(disable:6031)

#include<stdio.h>
#include<stdlib.h>

//Gauss-Jordan算法

void exchange(double* p1, double* p2, int n)
{
	int i;
	double tmp;
	for (i = 0; i < n; i++)
	{
		tmp = p1[i];
		p1[i] = p2[i];
		p2[i] = tmp;
	}
}

void primary_alternation(double** p1, int m, int n)
{
	int i, j, k1, k2;
	int tmp1;
	double tmp2;
	i = 0; j = 0;
	//进行初等行变换之后所得的p1矩阵，每行第一个非零项为1
	while (j < n)
	{
		tmp1 = i;
		while (!p1[i][j])
		{
			i++;
			if (i == m)
				i = tmp1, j++;
			if (j == n)
				break;
		}
		if (j == n)
			break;
		if (i != tmp1)
		{
			exchange(p1[tmp1], p1[i], n + 1);
			i = tmp1;
		}
		tmp2 = p1[i][j];
		for (k1 = j; k1 < n + 1; k1++)
			p1[i][k1] /= tmp2;
		for (k1 = i + 1; k1 < m; k1++)
		{
			if (p1[k1][j])
			{
				tmp2 = p1[k1][j];
				for (k2 = j; k2 < n + 1; k2++)
					p1[k1][k2] -= tmp2 * p1[i][k2];
			}
		}
		i++, j++;
		if (i == m)
			break;
	}
}

void Gauss_Jordan(const double** p, int m, int n)
{
	int i, j, k1, k2;
	int count = 0;
	double* x = (double*)malloc(sizeof(double) * n);
	int* random = (int*)malloc(sizeof(int) * n);
	double** p1 = (double**)malloc(sizeof(double*) * m);
	for (i = 0; i < m; i++)
		p1[i] = malloc(sizeof(double) * (n + 1));
	for (i = 0; i < m; i++)
		for (j = 0; j < n + 1; j++)
			p1[i][j] = p[i][j];//p1为导出方程组的系数矩阵
	primary_alternation(p1, m, n);//初等变换
	//确定  能够取任意值的未知数
	i = 0, j = 0;
	while (j < n)
	{
		if (i == m)
		{
			random[count] = j;
			j++, count++;
		}
		else if (!p1[i][j])
		{
			random[count] = j;
			j++, count++;
		}
		else
			i++, j++;
	}
	//判断是否有解
	 for (i = 0; i < m; i++)
	 {
		 for (j = 0; j < n; j++)
			 if (p1[i][j])
				 break;
		 if (j == n)
			 break;//这证明导出方程组有非零解
	 }
	 if (i == m)
		 printf("导出方程组没有非零解\n");
	//for (i = 0; i < m; i++)
	//{
	//	for (j = 0; j < n; j++)
	//		if (p1[i][j])
	//			break;
	//	if (j == n)
	//		if (p1[i][j])
	//		{
	//			printf("此方程无解\n");
	//			return;
	//		}
	//}
	 else
	 {
		 //求出导出方程组的解
		 printf("导出方程的基础解系如下:\n");
		 for (k1 = 0; k1 < count; k1++)
		 {
			 for (k2 = 0; k2 < count; k2++)
				 x[random[k2]] = (k2 == k1 ? 1 : 0);
			 i = m - 1, j = 0;
			 while (i >= 0)
			 {
				 while (!p1[i][j])
				 {
					 j++;
					 if (j == n)
						 i--, j = 0;
				 }
				 x[j] = 0;
				 for (k2 = j + 1; k2 < n; k2++)
					 x[j] -= x[k2] * p1[i][k2];
				 i--, j = 0;
			 }
			 //打印导出方程的基础解系
			 for (k2 = 0; k2 < n; k2++)
				 printf("%lf ", x[k2]);
			 printf("\n");
		 }
	 }
		//求出原方程的一个解
		for (k2 = 0; k2 < count; k2++)
			x[random[k2]] = 0;
		i = m - 1, j = 0;
		while (i >= 0)
		{
			while (!p1[i][j])
			{
				j++;
				if (j == n)
					i--, j = 0;
			}
			x[j] = p1[i][n];
			for (k2 = j + 1; k2 < n; k2++)
				x[j] -= x[k2] * p1[i][k2];
			i--, j = 0;
		}
		printf("原方程的一组解为:\n");
		for (k2 = 0; k2 < n; k2++)
			printf("%lf ", x[k2]);
		printf("\n");
		//if (!count)
		//{
		//	printf("方程仅有1解\n");
		//	i = m - 1, j = 0;
		//	while (i >= 0)
		//	{
		//		while (!p1[i][j])
		//		{
		//			j++;
		//			if (j == n)
		//				i--, j = 0;
		//		}
		//		x[j] = p1[i][n];
		//		for (k2 = j + 1; k2 < n; k2++)
		//			x[j] -= x[k2] * p1[i][k2];
		//		i--, j = 0;
		//	}
		//	for (k2 = 0; k2 < n; k2++)
		//		printf("%lf ", x[k2]);
		//	printf("\n");
		//}
	return;
}

int main()
{
	int m, n;
	int i, j;
	printf("请输入未知数个数:");
	scanf("%d", &n);
	printf("请输入等式个数:");
	scanf("%d", &m);
	double** p = NULL;
	p = (double**)malloc(sizeof(double*) * m);
	for (i = 0; i < m; i++)
		p[i] = (double*)malloc(sizeof(double) * (n + 1));
	printf("请输入系数矩阵:\n");
	for (i = 0; i < m; i++)
		for (j = 0; j < n + 1; j++)
			scanf("%lf", &p[i][j]);
	Gauss_Jordan(p, m, n);
	for (i = 0; i < m; i++)
		free(p[i]);
	free(p);
	return 0;
}



