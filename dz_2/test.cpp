#include <bits/stdc++.h>

using namespace std;
typedef std::complex<double> complexd;

int main(int argc, char *argv[])
{
    ifstream f1("res1.txt", ios::binary | ios::in), f2("res2.txt", ios::binary | ios::in), f4("res4.txt", ios::binary | ios::in), f8("res8.txt", ios::binary | ios::in);
    int qnum = atoi(argv[1]);
    unsigned long size = 1 << qnum;
    bool flag = true;
    complexd tmp1, tmp2, tmp4, tmp8;
    for (unsigned long i = 0; i < size; i++)
    {
        f1 >> tmp1;
        f2 >> tmp2;
        f4 >> tmp4;
        f8 >> tmp8;
        if (!(tmp1 == tmp2 && tmp4 == tmp8 && tmp1 == tmp8))
        {
            flag = false;
            cout << "Error! " << i << tmp1 << tmp2 << tmp4 << tmp8 << endl;
            break;
        }
    }
    if (flag)
    {
        cout << "Correct!" << endl;
    }

    f1.close();
    f2.close();
    f4.close();
    f8.close();
}