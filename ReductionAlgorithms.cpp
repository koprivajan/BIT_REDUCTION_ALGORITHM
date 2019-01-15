// ReductionAlgorithms.cpp : Defines the entry point for the console application.
//

#include <string>
#include <iostream>
#include <bitset>
#include <sstream>
#include <ctime>


typedef unsigned long long uint64;
typedef long long int64;

using namespace std;

template <class Function>
void time_call(Function&& f)
{
   clock_t start = clock();
   f();
   clock_t stop = clock();
   cout << "time : " << (stop-start) << " ms"<< endl;
}





//Knuth EGCD without negativ

unsigned int modInv(unsigned int u, unsigned int v)
{
    unsigned int inv, u1, u3, v1, v3, t1, t3, q;
    int iter;
    /* Step X1. Initialise */
    u1 = 1;
    u3 = u;
    v1 = 0;
    v3 = v;
    /* Remember odd/even iterations */
    iter = 1;
    /* Step X2. Loop while v3 != 0 */
    while (v3 != 0)
    {
        /* Step X3. Divide and "Subtract" */
        q = u3 / v3;
        t3 = u3 % v3;
        t1 = u1 + q * v1;
        /* Swap */
        u1 = v1; v1 = t1; u3 = v3; v3 = t3;
        iter = -iter;
    }
    /* Make sure u3 = gcd(u,v) == 1 */
    if (u3 != 1)
        return 0;   /* Error: No inverse exists */
    /* Ensure a positive result */
    if (iter < 0)
        inv = v - u1;
    else
        inv = u1;
    return inv;
}




string toBinary(unsigned long long n)
{
	string result;

	do 
	result.push_back( '0' + (n & 1) );
	while (n >>= 1);

	reverse( result.begin(), result.end() );
	return result;
}

/*
int MonRe(int aComma, int bComma)
{
	int t = aComma * bComma;
	int m = (t * mComma) % r;
	int u = (t + m * M ) * rExpMinusOne;
	while (u >= M)
	{
		u = u - M;
	}
	return u;
}
*/

uint64 Montgomery(uint64 A, uint64 B, uint64 M)
{
	string mNumber=toBinary(M);

	uint64 r=(uint64)pow(2,mNumber.size());
	uint64 rExpMinusOne = modInv(r,M);
	uint64 mComma= (r * rExpMinusOne - 1 ) / M;
	uint64 aComma= (A * r) % M;
	uint64 bComma= (B * r) % M;
// mon(a^, b^)
	uint64 t = aComma * bComma;
	uint64 m = (t * mComma) % r;
	uint64 u = (t + m * M ) / r;
	while (u >= M)
	{
		u = u - M;
	}
// mon(c^,1 )
	uint64 c = (u * 1 * rExpMinusOne) % M;
	
	
	return c;
}

int StandardMultiplication(int A, int B, int M)
{
	int X = A*B;
	while(X >= M)
	{
		X = X - M;  
	}
	return X;
}



uint64 Barrett(uint64 A, uint64 B, uint64 M)
{
	string mNumber=toBinary(M);
	uint64 X = A * B;

	uint64 micro = (uint64)pow(2,2*mNumber.size())/M;
/*	int x = (int)pow(2,2*mNumber.size());
	int micro = 0;
	while(x > M)
	{
		x = x - M;
		micro = micro +1;
	}
*/
	uint64 q1 = X/(uint64)pow(2,mNumber.size()-1);
	uint64 q2 = q1*micro;
	uint64 q3 = q2/(uint64)pow(2,mNumber.size()+1);
	uint64 r1 = X % (uint64)pow(2,mNumber.size()+1);
	uint64 r2 = (q3 * M) % (uint64)pow(2,mNumber.size()+1);
	int64 r = r1 - r2;
	if (r < 0)
	{
		r=r + (uint64)pow(2,mNumber.size()+1);
	}
	while(r >= M)
	{
		r = r - M;
	}
    return r;
}

uint64 Blakley(int64 A, int64 B, uint64 M)
{
	string aString=toBinary(A);
	uint64 r = 0;
	for(uint64 i = 0; i < aString.size(); i++)
	{
//conversion from string to int - ASCII
		unsigned int a = aString[i] - '0';		
		r = 2 * r + a*B;
		r = r % M;
	}
	return r;
}
//dalsi
uint64 RussianPeasant(uint64 a, uint64 b, uint64 m) {
    uint64 res = 0;
    while (a != 0) {
        if (a & 1) res = (res + b) % m;
        a >>= 1;
        b = (b << 1) % m;
    }
    return res;
}


//cisla A a B mohou byt i zaporna, cislo M je vzdy kladne
uint64 Kopriva(uint64 A, uint64 B, uint64 M)
{
//potrebujeme ziskat mensi cislo ze dvou predlozenich,
//protoze mensi z techto cisel - a - bude urcovat pocet korenu 
// vetsi cislo bude cislem korenu b
	if(A == 0 || B == 0)
	{
		return 0;
	}


//urychleni, otoceni cisel aby se nasobilo mensi cislo vetsim
	if(A > B)
	{
		uint64 v = B;
		B = A;
		A = v;
	}


	uint64 ltree=0;
	uint64 r = 0;
	uint64 targetLevel = 0;
	uint64 ftree=B;

/*
	if((A > 1) && (A & 1))
	{
//zmenseni cisla o jedna
//posledni cislice z b, ktera se nevleze do nasobku d	
		A=A-1;
		ltree = B;
	}
//zjisteni radu cisla a, ktere slouzi k vypoctu delky stromu - logaritmu
	unsigned int ftree=B;
	string aNumber=toBinary(A);
	unsigned int targetLevel = aNumber.size()-1;

//ladenie
//ziskame pocet clenu vedlejsiho stromu
//		int r=A-(pow(2,targetLevel));
	
	if(A > 1)
	{
		r=A- (2 << targetLevel);
	}
*/
	if((A > 1))
	{
//zmenseni cisla o jedna
//posledni cislice z b, ktera se nevleze do nasobku d	
		if((A & 1))
		{
			A=A-1;
			ltree = B;
		}
//zjisteni radu cisla a, ktere slouzi k vypoctu delky stromu - logaritmu

		string aNumber=toBinary(A);
		targetLevel = aNumber.size()-1;
		r=A- (2 << targetLevel);
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
	}


//jednotlivych vrstech je stejny pocet jako log2 a
// protoze prvni rad sme jiz ulozili zbyva log2 a - 1
//mozno dar do podminky A > 1
	for(uint64 i = 0; i < targetLevel; i++)
	{
		ftree=(ftree << 1);
		if(ftree >= M)
		{
			ftree -= M;
		}
//pokud je rozklad cisla A mocninou dve - pak nepocitame pravy hreben stromu 
//a pocet radu stromu mensi o jednu protoze prvni rad byl konstruovan pred for cyklem
		if (  (i != targetLevel-1))
		{
			uint64 exponent = (2 << i); 
			if(r & exponent)
			{
				ltree=ltree + ftree;
				if(ltree >=M)
				{
					ltree -= M;
				}

			}
		}
	}


//////////////////////////////////


		uint64 result = (ftree + ltree);
		if(result>=M)
		{
			result -=M;
		}
		return result; 
	

}

int main()
{

//31,61,127,509,1021,2039  b834005 854848  k734116 728911 m853162 827476

int prime = 2039;
	

	time_call([&]{
	
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j < prime; j++)
		{

			if((i*j % prime) != Kopriva(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

		}
	}

	});

	time_call([&]{
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j <prime; j++)
		{
			if((i*j % prime) != Barrett(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

//			Barrett(i,j,prime);
		}
	}

	});

	time_call([&]{
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j <prime; j++)
		{
			if((i*j % prime) != Montgomery(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

		}
	}

	});

	time_call([&]{
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j <prime; j++)
		{
			if((i*j % prime) != StandardMultiplication(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

		}
	}

	});

	time_call([&]{
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j <prime; j++)
		{
			if((i*j % prime) != Blakley(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

		}
	}

	});

	time_call([&]{
	for(int i=0; i<prime; i++)
	{
		for(int j=0; j <prime; j++)
		{
			if((i*j % prime) != RussianPeasant(i,j,prime))
			{
				cout << i << " " << j << endl;
			}

		}
	}

	});


/*
	int A = 32748;
	int B = 32748;
	int M = 32749;

	time_call([&]{
	
	int X = A*B%M;

	});

	time_call([&]{

	int resultKopriva=Kopriva(A,B,M);

	});

	time_call([&]{

	int resultMontgomery=Montgomery(A,B,M);

	});

*/	



/*
	int X = A * B;
	long long result = (A*B) % M;

	int resultBarett=Barrett(A,B,M);
	int resultStandard=StandardMultiplication(A,B,M);
	int resultBlakley=Blakley(A,B,M);
	int resultMontgomery=Montgomery(A,B,M);
	int resultKopriva=Kopriva(A,B,M);
*/
	return 0;
}

/*
//zkouska
// M^2, bin mul a add, polynom
// test
	long long stdr = (A*B)/M;
	long long zme = (a*b)/M;
	long long rozdil = stdr - zme;
*/
	//dalsi pokus

//f6997	
//	int A = 3499;
//f17	int A = 16;
//f17	int A = 12;
//f23	int A = 12;

//f6997	
//	int B = 3498;
//f17	int B = 15;
//f17	int B = 16;
//f23    int B = 20;

//f6997	
//	int M = 6997;
//	int M = 17;
//	int M = 23;




/*
	int a;
	int b;

	bool aSign;
	bool bSign;
	long long int x;
	long int res;
	if (A > (M / 2))
	{
		a = M - A;
		aSign=1;
	}
	else
	{
		a = A;
		aSign = 0;
	}


	if (B > (M / 2))
	{
		b = M - B;
		bSign = 1;
		
	}
	else
	{
		b = B;
		bSign = 0;
	}
// zajima nas situace, kdy znamenko jednoho je rozdilne od druheho a <> b
if(aSign != bSign)
{
	changeSign = 1;
}
else
{
	changeSign = 0;
}


	
// oba jsou mensi jak nula
	
	if((aPlus == 0 && bPlus == 0) || (aPlus == 1 && bPlus == 1))
	{
		 
		while(x > M)
		{
			x = x - M; 
		}
	}
// jedno je kladne druhe zaporne
	else
	{
		while(x < 0)
		{
			x = x + M;
		}
	}
	res = x;
*/

