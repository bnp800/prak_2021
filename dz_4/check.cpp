#include<iostream>
#include<complex>
#include<mpi.h>
#include<string>
#include<cstdlib>
#include"check.h"

typedef std::complex<double> complexd;

using namespace std;

int main(int argc, char *argv[])
{
    int qnum = atoi(argv[1]);
    string name = argv[2];
    int rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int size = 1 << qnum;
    int part = size / world_size;
    int myleft = rank * part;
    int myright = (rank + 1) * part - 1;
    //cout << "process " << rank << " has part: " << part << " of left " << myleft << " with qnum " << qnum << " and size " << size << endl;
    //check("res1.txt", "res2.txt", "res4.txt", "res8.txt", part, myleft);
    //check("res1_noise.txt", "res2_noise.txt", "res4_noise.txt", "res8_noise.txt", part, myleft);
    check(name, "vector.txt", part, myleft);
    MPI_Finalize();
}