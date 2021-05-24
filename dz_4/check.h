#include <iostream>
#include <complex>
#include <fstream>
#include <mpi.h>

typedef std::complex<double> complexd;

using namespace std;

void check(string name1, string name2, int size, int left)
{
    MPI_File f1, f2; //, f4, f8;
    MPI_File_open(MPI_COMM_WORLD, name1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f1);
    MPI_File_open(MPI_COMM_WORLD, name2.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f2);
    //MPI_File_open(MPI_COMM_WORLD, name3.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f4);
    //MPI_File_open(MPI_COMM_WORLD, name4.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &f8);

    bool flag = true;
    //complexd tmp1[size], tmp2[size]; //tmp4[size], tmp8[size];
    complexd *tmp1, *tmp2;
    tmp1 = new complexd[size];
    tmp2 = new complexd[size];
    MPI_File_read_ordered(f1, tmp1, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_read_ordered(f2, tmp2, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    //MPI_File_read_ordered(f4, tmp4, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    //MPI_File_read_ordered(f8, tmp8, size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);

    for (int i = 0; i < size; i++)
    {
        //if (!(tmp1[i] == tmp2[i] && tmp4[i] == tmp8[i] && tmp1[i] == tmp8[i]))
        if (abs(tmp1[i] - tmp2[i]) > 0.0000001) // && tmp2[i] == tmp4[i]))
        {
            flag = false;
            //cout << "Error! " << i << " " <<  tmp1[i] << " " << tmp2[i] << " " << tmp4[i] << tmp8[i] << endl;
            cout << "Error! " << left << " on i " << i << " : " << tmp1[i] << " " << tmp2[i] << endl; // << " " << tmp4[i] << endl;
            break;
        }
    }
    if (flag)
    {
        cout << "Correct!" << endl;
    }

    delete[] tmp1;
    delete[] tmp2;
    MPI_File_close(&f1);
    MPI_File_close(&f2);
    //MPI_File_close(&f4);
}