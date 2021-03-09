#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

typedef std::complex<double> complexd;
#define QubitNum 20
#define NUM_THREAD 8
void Qubit(complexd *in, complexd *out, complexd U[2][2], int nqubits, int k)
{
    int shift = nqubits - k;
    int pow2q = 1 << (shift);
    int N = 1 << nqubits;
#pragma omp parallel for num_threads(NUM_THREAD)
    for (int i = 0; i < N; i++)
    {
        int i0 = i & (~pow2q);
        int i1 = i | pow2q;
        int iq = (i & pow2q) >> shift;
        out[i] = U[iq][0] * in[i0] + U[iq][1] * in[i1];
    }
}
void init_vector(complexd *in)
{
	int size = 1 << QubitNum;
#pragma omp parallel for num_threads(128)
    for (int i = 0; i < (size - 4); i+=4)
    {
        complexd temp1(rand(), rand()),temp2(rand(),rand()),temp3(rand(),rand()),temp4(rand(),rand());
        in[i] = temp1;
	in[i+1] = temp2;
	in[i+2] = temp3;
	in[i+3] = temp4;
    }
}
int main()
{
    srand(time(NULL));
    complexd U[2][2], *in, *out;
    unsigned long long int size = 1;
    int k;
    //k = 9; //number in list
    //k = 1;
    k = QubitNum;
    size <<= QubitNum;
    U[0][0] = U[0][1] = U[1][0] = sqrt(2) / 2;
    U[1][1] = -sqrt(2) / 2;
    in = new complexd[size];
    out = new complexd[size];
    init_vector(in);
    double begin = omp_get_wtime();
    Qubit(in, out, U, QubitNum, k);
    double end = omp_get_wtime();
    cout << "Thread num: " << NUM_THREAD << endl
         << "Qubit: " << QubitNum << endl
         << "Time: " << end - begin << " s" << endl;
    delete[] in;
    delete[] out;

    return 0;
}
