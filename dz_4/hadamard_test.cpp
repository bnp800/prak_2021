#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

typedef std::complex<double> complexd;

void Qubit(complexd *op1, complexd *op2, complexd *out, complexd U[2][2], int nqubits, int size, int left, int k)
{
    int shift = nqubits - k;
    int pow2q = 1 << shift;
    //cout << i + left << endl;
    //int i0 = (i + left) & (~pow2q);
    //int i1 = (i + left) | pow2q;

    for (int i = 0; i < size; i++)
    {
        int iq = ((i + left) & pow2q) >> shift;
        out[i] = U[iq][0] * op1[i] + U[iq][1] * op2[i];
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
    if (myright >= size)
    {
        myright = size - 1;
    }
    partion_size = myright - myleft + 1;

    in = new complexd[partion_size];
    out = new complexd[partion_size];

    int *indexleft = new int[world_size], *indexright = new int[world_size];
    indexleft[myrank] = myleft;
    indexright[myrank] = myright;
    MPI_Allgather(indexleft + myrank, 1, MPI_INT, indexleft, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(indexright + myrank, 1, MPI_INT, indexright, 1, MPI_INT, MPI_COMM_WORLD);
    string name = "H_res" + to_string(world_size) + ".txt";
    MPI_File fin, fout;
    MPI_File_open(MPI_COMM_WORLD, name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
    MPI_File_open(MPI_COMM_WORLD, string(name + ".t").c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fout);
    MPI_File_read_ordered(fin, in, partion_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    complexd *op1, *op2;
    op1 = new complexd[partion_size];
    op2 = new complexd[partion_size];

    for (int cur_rank = 0; cur_rank < world_size; cur_rank++)
    {
        if (myrank == cur_rank)
        {
            //cout << "Now process " << cur_rank << " colletcting data" << endl;
            unsigned int shift = qnum - target;
            unsigned int pow2q = 1 << shift;
            for (int i = 0; i < partion_size; i++)
            {
                unsigned int i0 = (i + myleft) & (~pow2q);
                unsigned int i1 = (i + myleft) | pow2q;
                if (i0 >= myleft && i0 <= myright)
                {
                    op1[i] = in[i0 - myleft];
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i0 || indexright[rank] < i0)
                        rank--;
                    //cout << "Process " << myrank << " : need data index " << i0 << " for op1 from process " << rank << " which now it's " << op1[i] << endl;
                    MPI_Sendrecv(&i0, 1, MPI_INT, rank, 0, &op1[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i0 << " for op1 from process " << rank << " which now it's " << op1[i] << endl;
                }
                if (i1 >= myleft && i1 <= myright)
                {
                    op2[i] = in[i1 - myleft];
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i1 || indexright[rank] < i1)
                        rank--;
                    //cout << "Process " << myrank << " : need data index " << i1 << " for op2 from process " << rank << " which now it's " << op2[i] << endl;
                    MPI_Sendrecv(&i1, 1, MPI_INT, rank, 0, &op2[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i1 << " for op2 from process " << rank << " which now it's " << op2[i] << endl;
                }
            }
            for (int i = 0; i < world_size; i++)
            {
                if (i != myrank)
                {
                    int id = -1;
                    MPI_Send(&id, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            }
        }
        else
        {
            //cout << "Process " << myrank << " sending data to process " << cur_rank << /*" flag = " << flag << */ endl;
            int target_id = 0;
            while (1)
            {
                //cout << "Process " << myrank << " wait message from process " << cur_rank << endl;
                MPI_Recv(&target_id, 1, MPI_INT, cur_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (target_id == -1)
                    break;
                else
                    MPI_Send(in + (target_id - myleft), 1, MPI_DOUBLE_COMPLEX, cur_rank, 0, MPI_COMM_WORLD);
                //cout << "Process " << myrank << " : sent data index " << target_id << " to process " << cur_rank << " which now it's " << in[target_id - myleft] << endl;
            }
        }
    }

    my_begin = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    Qubit(op1, op2, out, U, qnum, partion_size, myleft, target);
    MPI_Barrier(MPI_COMM_WORLD);
    my_end = MPI_Wtime();
    my_time = my_end - my_begin;
    MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    total_time /= world_size;
    MPI_File_write_ordered(fout, out, partion_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_close(&fin);
    MPI_File_close(&fout);
    //MPI_Gatherv(out, partion_size, MPI_DOUBLE_COMPLEX, all, recvcount, offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    //name = "res" + to_string(world_size) + to_string(myrank) + ".txt";
    //Write(out, size, name);
    /*if (!myrank)
    {

        cout << "Qnum: " << qnum << endl;
        cout << "Pnum: " << world_size << endl;
        cout << "Target: " << target << endl;
        cout << "Time: " << total_time << " s" << endl;
    }*/
    delete[] op1;
    delete[] op2;
    delete[] indexleft;
    delete[] indexright;

    delete[] in;
    delete[] out;

    MPI_Finalize();
}