#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

typedef std::complex<double> complexd;
void Qubit(complexd *in, complexd *out, complexd U[2][2], int nqubits, int left, int right, int k)
{
    int shift = nqubits - k;
    int pow2q = 1 << (shift);
    int N = (1 << nqubits);
    for (int i = 0; i < right - left; i++)
    {
        int i0 = i & (~pow2q);
        int i1 = i | pow2q;
        int iq = (i & pow2q) >> shift;
        out[i] = U[iq][0] * in[i0] + U[iq][1] * in[i1];
    }
}

void Init_vector(complexd *in, int left, int right)
{
    for (int i = 0; i < right - left; i++)
    {
        in[i] = complexd(rand(), rand());
    }
}

double Get_sum(complexd *in, int left, int right)
{
    double res = 0;
    for (int i = 0; i < right - left; i++)
    {
        res += abs(in[i] * in[i]);
    }
    return res;
}

void Normalize(complexd *in, int left, int right, int sum)
{
    for (int i = 0; i < right - left; i++)
    {
        in[i] /= sum;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_size, myrank, qnum = atoi(argv[1]), target = atoi(argv[2]);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    unsigned long size = 1 << qnum;
    complexd U[2][2], *in, *out;
    U[0][0] = U[0][1] = U[1][0] = sqrt(2) / 2;
    U[1][1] = -sqrt(2) / 2;
    double my_begin, my_end, my_time, total_time;
    int partion_size = size / world_size;
    int myleft = myrank * partion_size;
    int myright = (myrank + 1) * partion_size - 1;
    int totalsum;
    if (myright >= size)
    {
        myright = size - 1;
    }
    in = new complexd[myright - myleft];
    out = new complexd[myright - myleft];

    srand(myrank * time(NULL));
    Init_vector(in, myleft, myright);
    double sum = Get_sum(in, myleft, myright);
    MPI_Reduce(&sum, &totalsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    Normalize(in, myleft, myright, totalsum);
    my_begin = MPI_Wtime();
    Qubit(in, out, U, qnum, myleft, myright, target);
    my_end = MPI_Wtime();
    my_time = my_end - my_begin;
    MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!myrank)
    {
        cout << "Qnum: " << qnum << endl;
        cout << "Time: " << total_time << endl;
        cout << "Target: " << target << endl;
        cout << "Pnum: " << world_size << endl;
    }
    //Init_vector(myrank, qnum);
    MPI_Finalize();
}