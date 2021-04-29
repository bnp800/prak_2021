#include <bits/stdc++.h>

using namespace std;
typedef std::complex<double> complexd;

void check(string name1, string name2, string name3, string name4, int qnum)
{
    ifstream f1(name1, ios::binary | ios::in), f2(name2, ios::binary | ios::in), f4(name3, ios::binary | ios::in), f8(name4, ios::binary | ios::in);
    //int qnum = atoi(argv[1]);
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

int main(int argc, char* argv[]) {
    int qnum = atoi(argv[1]);
    check("res1.txt", "res2.txt", "res4.txt", "res8.txt", qnum);
    check("res1_noise.txt", "res2_noise.txt", "res4_noise.txt", "res8_noise.txt", qnum);
    return 0;
}