#include<iostream>
using namespace std;

int main()
{
	int i;
	int sum = 0;
	for (i = 0; i < 1000000000; i++)
	{
		sum += i;
	}
	cout << sum;

}