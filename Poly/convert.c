/* 
	convert.c -- a program that is used to convert .002 file to .plt file 
	plot.bat -- output file
*/
#include<stdio.h>

int main ()
{
	FILE *fp_out;
	int  i, j, a[12],b[12], c[101];
// open output file
	fp_out = fopen("plot.bat","w");

for  (i = 1 ; i <12; i++)
{
	if (i >= 10)
	{
		a[i]=1;
		b[i] = 0;
	}
	else 
	{
		a[i] = 0;
		b[i] = i;
	}
	if (i==11)
	{
		a[i]=1; b[i]=1;
	}
	for (j= 1 ; j < 101; j++)
	{		if (j <10)
		{
			c[j]=0;
			fprintf(fp_out, "\n00poly.exe pol%d%d-%d%d%d.002", a[i], b[i], c[j],c[j],j);	
		}
		else
		{
			if (j<100)
			{ 
				c[j]=0;
				fprintf(fp_out, "\n00poly.exe pol%d%d-%d%d.002", a[i], b[i], c[j],j);
			}
			else
				fprintf(fp_out, "\n00poly.exe pol%d%d-100.002", a[i], b[i]);
		}

	}
}
fclose (fp_out);
return (0);
}