#include "itensor/all.h"

using namespace itensor;


int main()
    {

	auto i = Index(4);
	auto j = prime(i);

	auto T = ITensor(i,j);

	PrintData(T);

	T.set(i=1,j=1, 101);
	T.set(i=1,j=4, 12);
	T.set(i=2,j=2, 2);
	T.set(i=3,j=3, 3.14);
	T.set(i=4,j=4, -10);
	T.set(i=4,j=1, 12);

	PrintData(T);

	auto [U,D] = diagHermitian(T); 

	PrintData(U);
	PrintData(D);

    return 0;
    }
