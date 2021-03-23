#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

typedef std::complex<double> complexd;
void Qubit(complexd *in, complexd *out, complexd U[2][2], int nqubits, int size, int left, int k)
{
    int shift = nqubits - k;
    int pow2q = 1 << (shift);
    for (int i = 0; i < size; i++)
    {
        int i0 = (i + left) & (~pow2q);
        int i1 = (i + left) | pow2q;
        int iq = ((i + left) & pow2q) >> shift;
        out[i] = U[iq][0] * in[i0] + U[iq][1] * in[i1];
    }
}

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

void Write(complexd *in, unsigned long size, string name)
{
    ofstream myout;
    myout.open(name);
    for (unsigned long i = 0; i < size; i++)
    {
        myout << fixed << setprecision(2) << in[i] << " ";
    }
    myout.close();
}

void Read(complexd *in, unsigned long size, string name)
{
    ifstream myin(name);
    for (unsigned long i = 0; i < size; i++)
    {
        myin >> in[i];
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_size, myrank, qnum = atoi(argv[1]), target = atoi(argv[2]);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    unsigned long size = 1 << qnum;
    complexd U[2][2], *in, *out, *all;
    U[0][0] = U[0][1] = U[1][0] = sqrt(2) / 2;
    U[1][1] = -sqrt(2) / 2;
    double my_begin, my_end, my_time, total_time;
    int partion_size = size / world_size;
    int myleft = myrank * partion_size;
    int myright = (myrank + 1) * partion_size - 1;
    double totalsum;
    int mode = atoi(argv[3]);
    if (myright >= size)
    {
        myright = size - 1;
    }
    partion_size = myright - myleft + 1;
    int *recvcount, *offset;
    recvcount = new int[world_size];
    offset = new int[world_size];
    offset[myrank] = myrank * partion_size;
    recvcount[myrank] = partion_size;

    MPI_Allgather(offset + myrank, 1, MPI_INT, offset, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(recvcount + myrank, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);

    in = new complexd[partion_size];
    out = new complexd[partion_size];
    all = new complexd[size];

    if (mode == 1)
    {
        string name = "in.txt";
        srand(myrank * time(NULL));
        Init_vector(in, partion_size);
        double sum = Get_sum(in, partion_size);
        MPI_Reduce(&sum, &totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&totalsum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        Normalize(in, partion_size, totalsum);
        MPI_Allgatherv(in, partion_size, MPI_DOUBLE_COMPLEX, all, recvcount, offset, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
        Write(all, size, name);
    }

    if (mode == 2)
    {
        string name = "in.txt";
        Read(all, size, name);
        my_begin = MPI_Wtime();
        Qubit(all, out, U, qnum, partion_size, myleft, target);
        my_end = MPI_Wtime();
        my_time = my_end - my_begin;
        MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Gatherv(out, partion_size, MPI_DOUBLE_COMPLEX, all, recvcount, offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        if (!myrank)
        {
            string name = "res" + to_string(world_size) + ".txt";
            cout << "Qnum: " << qnum << endl;
            cout << "Pnum: " << world_size << endl;
            cout << "Target: " << target << endl;
            cout << "Time: " << total_time << " s" << endl;
            Write(all, size, name);
        }
    }
    delete[] recvcount;
    delete[] offset;
    delete[] all;

    MPI_Finalize();
}