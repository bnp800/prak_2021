#include <iostream>
#include <complex>
#include <fstream>
#include<cstdlib>
#include <mpi.h>

using namespace std;

typedef std::complex<double> complexd;

void Init_vector(complexd *in, unsigned long size)
{
    for (unsigned long i = 0; i < size; i++)
    {
        in[i] = complexd(rand(), rand());
    }
}

double Get_sum(complexd *in, unsigned long size)
{
    double res = 0;
    for (unsigned long i = 0; i < size; i++)
    {
        res += abs(in[i] * in[i]);
    }
    return res;
}

void Normalize(complexd *in, unsigned long size, double sum)
{
    for (unsigned long i = 0; i < size; i++)
    {
        in[i] /= sqrt(sum);
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_size, myrank, qnum = atoi(argv[1]);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    unsigned long size = 1 << qnum;
    complexd U[2][2], *in;
    int partion_size = size / world_size;
    int myleft = myrank * partion_size;
    int myright = (myrank + 1) * partion_size - 1;
    double totalsum;
    if (myright >= size)
    {
        myright = size - 1;
    }
    partion_size = myright - myleft + 1;

    in = new complexd[partion_size];

    srand(myrank * time(NULL));
    Init_vector(in, partion_size);
    double sum = Get_sum(in, partion_size);
    MPI_Reduce(&sum, &totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&totalsum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    Normalize(in, partion_size, totalsum);

    MPI_File fout;
    MPI_File_open(MPI_COMM_WORLD, "vector.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fout);
    MPI_File_write_ordered(fout, in, partion_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_close(&fout);
    MPI_Finalize();
}