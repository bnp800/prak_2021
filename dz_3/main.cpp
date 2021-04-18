#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

typedef complex<double> complexd;

double normal_dis_gen()
{
    double S = 0.;
#pragma omp parallel for
    for (int i = 0; i < 12; ++i)
    {
        S += (double)rand() / RAND_MAX;
    }
    return S - 6.;
}

double fidelity(complexd *a, complexd *b, int size)
{
    complexd sum = 0;
#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        sum += a[i] * b[i];
    }
    return norm(sum);
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Finalize();
}