#include <iostream>
#include <mpi.h>
#include <complex>

using namespace std;
typedef std::complex<double> complexd;

void check(string name1, string name2, string name3, string name4, int size, int left)
{
    MPI_File f1, f2, f4, f8;
    MPI_File_open(MPI_COMM_WORLD, name1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f1);
    MPI_File_open(MPI_COMM_WORLD, name2.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f2);
    MPI_File_open(MPI_COMM_WORLD, name3.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f4);
    //MPI_File_open(MPI_COMM_WORLD, name4.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f8);

    bool flag = true;
    complexd tmp1[size], tmp2[size], tmp4[size], tmp8[size];
    MPI_File_read_ordered(f1, tmp1, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_read_ordered(f2, tmp2, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_read_ordered(f4, tmp4, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    //MPI_File_read_ordered(f8, tmp8, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);

    for (unsigned long i = 0; i < size; i++)
    {
        //if (!(tmp1[i] == tmp2[i] && tmp4[i] == tmp8[i] && tmp1[i] == tmp8[i]))
        if (!(tmp1[i] == tmp2[i] && tmp2[i] == tmp4[i]))
        {
            flag = false;
            //cout << "Error! " << i << " " <<  tmp1[i] << " " << tmp2[i] << " " << tmp4[i] << tmp8[i] << endl;
            cout << "Error! " << left << " " << i + left << " " << tmp1[i] << " " << tmp2[i] << " " << tmp4[i] << endl;
            break;
        }
    }
    if (flag)
    {
        cout << "Correct!" << endl;
    }

    MPI_File_close(&f1);
    MPI_File_close(&f2);
    MPI_File_close(&f4);
}

int main(int argc, char *argv[])
{
    int qnum = atoi(argv[1]);
    int rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int size = 1 << qnum;
    int part = size / world_size;
    int myleft = rank * part;
    int myright = (rank + 1) * part - 1;
    if (myright >= size)
    {
        myright = size - 1;
    }
    part = myright - myleft + 1;
    check("res1.txt", "res2.txt", "res4.txt", "res8.txt", part, myleft);
    //check("res1_noise.txt", "res2_noise.txt", "res4_noise.txt", "res8_noise.txt", part, myleft);
    MPI_Finalize();
}