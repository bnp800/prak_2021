#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

typedef std::complex<double> complexd;

void Qubit(complexd *op1, complexd *op2, complexd *op3, complexd *op4, complexd *out, complexd U[4][4], int nqubits, int size, int left, int q1, int q2)
{
    int shift1 = nqubits - q1;
    int shift2 = nqubits - q2;
    int pow2q1 = 1 << (shift1);
    int pow2q2 = 1 << (shift2);

    //cout << i + left << endl;
    //int i0 = (i + left) & (~pow2q);
    //int i1 = (i + left) | pow2q;

    for (int i = 0; i < size; i++)
    {
        int iq1 = ((i + left) & pow2q1) >> shift1;
        int iq2 = ((i + left) & pow2q2) >> shift2;
        int iq = (iq1 << 1) + iq2;

        //cout << "iq: " << iq << endl;
        out[i] = U[iq][(0 << 1) + 0] * op1[i] + U[iq][(0 << 1) + 1] * op2[i] + U[iq][(1 << 1) + 0] * op3[i] + U[iq][(1 << 1) + 1] * op4[i];
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_size, myrank, qnum = atoi(argv[1]), target1 = atoi(argv[2]), target2 = atoi(argv[3]);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    unsigned long size = 1 << qnum;
    complexd U[4][4], *in, *out, *all;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            U[i][j] = 0;
    U[0][0] = U[1][1] = U[2][3] = U[3][2] = 1;

    double my_begin, my_end, my_time, total_time;
    int partion_size = size / world_size;
    int myleft = myrank * partion_size;
    int myright = (myrank + 1) * partion_size - 1;
    double totalsum;
    /*if (myright >= size)
    {
        myright = size - 1;
    }*/
    //partion_size = myright - myleft + 1;

    in = new complexd[partion_size];
    out = new complexd[partion_size];

    int *indexleft = new int[world_size], *indexright = new int[world_size];
    indexleft[myrank] = myleft;
    indexright[myrank] = myright;
    MPI_Allgather(indexleft + myrank, 1, MPI_INT, indexleft, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(indexright + myrank, 1, MPI_INT, indexright, 1, MPI_INT, MPI_COMM_WORLD);

    string name = "CNOT_res" + to_string(world_size) + ".txt";
    MPI_File fin, fout;
    MPI_File_open(MPI_COMM_WORLD, "vector.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
    MPI_File_open(MPI_COMM_WORLD, name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fout);
    MPI_File_read_ordered(fin, in, partion_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    complexd *op1, *op2, *op3, *op4;
    op1 = new complexd[partion_size];
    op2 = new complexd[partion_size];
    op3 = new complexd[partion_size];
    op4 = new complexd[partion_size];

    for (int cur_rank = 0; cur_rank < world_size; cur_rank++)
    {
        if (myrank == cur_rank)
        {
            //cout << "Now process " << cur_rank << " colletcting data" << endl;
            int shift1 = qnum - target1;
            int shift2 = qnum - target2;
            int pow2q1 = 1 << (shift1);
            int pow2q2 = 1 << (shift2);

            for (int i = 0; i < partion_size; i++)
            {
                int i00 = (i + myleft) & ~pow2q1 & ~pow2q2;
                int i01 = (i + myleft) & ~pow2q1 | pow2q2;
                int i10 = ((i + myleft) | pow2q1) & ~pow2q2;
                int i11 = (i + myleft) | pow2q1 | pow2q2;

                //cout << "process " << myrank << " has "<< "i00: " << i00 << ", i01: " << i01 << ", i10: " << i10 << ", i11: " << i11 << ", myleft" << myleft << ", myright" << myright << endl;
                //cout << "now collecting for op1" << endl;
                if (i00 >= myleft && i00 <= myright)
                {
                    op1[i] = in[i00 - myleft];
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i00 || indexright[rank] < i00)
                        rank--;
                    //cout << "op1 need data from process " << rank << endl;
                    //cout << "Process " << myrank << " : need data index " << i00 << " for op1 from process " << rank << " which now it's " << op1[i] << endl;
                    MPI_Sendrecv(&i00, 1, MPI_INT, rank, 0, &op1[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i00 << " for op1 from process " << rank << " which now it's " << op1[i] << endl;
                }

                //cout << "now collecting for op2" << endl;
                if (i01 >= myleft && i01 <= myright)
                {
                    op2[i] = in[i01 - myleft];
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i01 || indexright[rank] < i01)
                        rank--;
                    //cout << "op2 need data from process " << rank << endl;
                    //cout << "Process " << myrank << " : need data index " << i01 << " for op2 from process " << rank << " which now it's " << op2[i] << endl;
                    MPI_Sendrecv(&i01, 1, MPI_INT, rank, 0, &op2[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i01 << " for op2 from process " << rank << " which now it's " << op2[i] << endl;
                }

                //cout << "now collecting for op3" << endl;
                if (i10 >= myleft && i10 <= myright)
                {
                    op3[i] = in[i10 - myleft];
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i10 || indexright[rank] < i10)
                        rank--;
                    //cout << "op3 need data from process " << rank << endl;
                    //cout << "Process " << myrank << " : need data index " << i10 << " for op3 from process " << rank << " which now it's " << op3[i] << endl;
                    MPI_Sendrecv(&i10, 1, MPI_INT, rank, 0, &op3[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i10 << " for op3 from process " << rank << " which now it's " << op3[i] << endl;
                }

                //cout << "now collecting for op4" << endl;
                if (i11 >= myleft && i11 <= myright)
                {
                    op4[i] = in[i11 - myleft];
                    //cout << "op4 index: " << i11 - myleft << endl;
                }
                else
                {
                    int rank = world_size - 1;
                    while (indexleft[rank] > i11 || indexright[rank] < i11)
                        rank--;
                    //cout << "op4 need data from process " << rank << endl;
                    //cout << "Process " << myrank << " : need data index " << i11 << " for op4 from process " << rank << " which now it's " << op4[i] << endl;
                    MPI_Sendrecv(&i11, 1, MPI_INT, rank, 0, &op4[i], 1, MPI_DOUBLE_COMPLEX, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "Process " << myrank << " : had data index " << i11 << " for op4 from process " << rank << " which now it's " << op4[i] << endl;
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
    //cout << "collection complete" << endl;
    my_begin = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    Qubit(op1, op2, op3, op4, out, U, qnum, partion_size, myleft, target1, target2);
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
    if (!myrank)
    {

        cout << "Qnum: " << qnum << endl;
        cout << "Pnum: " << world_size << endl;
        //cout << "Target: " << target << endl;
        cout << "Time: " << total_time << " s" << endl;
    }
    delete[] op1;
    delete[] op2;
    delete[] op3;
    delete[] op4;
    delete[] indexleft;
    delete[] indexright;

    delete[] in;
    delete[] out;

    MPI_Finalize();
}