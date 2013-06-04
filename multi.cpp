#include "boost/multi_array.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>

using namespace std;

typedef boost::multi_array<int, 3> score_type;

const int MATCH = 5;
const int MISMATCH = -4;
const int INDEL = -8;

int n2n(char c)
{
    switch(c)
    {
        case 'A':
            return 0;
        case 'T':
            return 1;
        case 'G':
            return 2;
        case 'C':
            return 3;
        default:
            return 4;
    }
}

vector<int> s1;
vector<int> s2;
vector<int> s3;

int scoring[5][5][5];

struct point3d
{
    int x;
    int y;
    int z;
};
vector<point3d> midpoint;

bool ranking(const point3d &lhs, const point3d &rhs)
{
    return lhs.x < rhs.x;
}

bool maxvalue = true;
int length = 0;
int pm = 0;

void createscoring();
void alignment(int, int, int, int, int, int);
void backtrace(int, int, int, int, int, int);

int main(int argc, char ** argv)
{
    createscoring();

    ifstream fseq1(argv[1]);
    ifstream fseq2(argv[2]);
    ifstream fseq3(argv[3]);
    if(fseq1.fail()||fseq2.fail()||fseq3.fail())
    {
        cout << "Open sequence files failed!" << endl;
    }

    string line;

    while(getline(fseq1,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            s1.push_back(n2n(line[i]));
        }
    }

    while(getline(fseq2,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            s2.push_back(n2n(line[i]));
        }
    }

    while(getline(fseq3,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            s3.push_back(n2n(line[i]));
        }
    }

    cout << "sequence size : " << s1.size() << " " << s2.size() << " " << s3.size() << endl;

    point3d pp;
    pp.x = 0;
    pp.y = 0;
    pp.z = 0;
    midpoint.push_back(pp);
    pp.x = s1.size()-1;
    pp.y = s2.size()-1;
    pp.z = s3.size()-1;
    midpoint.push_back(pp);


    clock_t start = clock();
    alignment(0, 0, 0, s1.size()-1, s2.size()-1, s3.size()-1);
    clock_t end = clock();
    cout << "Execution time : " << (double)(end - start)/CLOCKS_PER_SEC << endl;

    sort(midpoint.begin(), midpoint.end(), &ranking);
    //for(int i = 0; i < midpoint.size()-1; i++)
    for(int i = 0; i < 5; i++)
    {
        point3d pp = midpoint[i];
        point3d qq = midpoint[i+1];
        printf("Find the crossing point (%d, %d, %d) to (%d, %d, %d)\n", pp.x, pp.y, pp.z, qq.x, qq.y, qq.z);
        if(pp.x == qq.x && pp.y == qq.y && pp.z == qq.z)
            continue;
        //else if(pp.x)
        backtrace(pp.x, pp.y, pp.z, qq.x, qq.y, qq.z);
    }
    cout << "Alignment length: " << length << endl;

    return 0;
}


void backtrace(int si, int sj, int sk, int ei, int ej, int ek)
{
    int li = ei - si + 1;
    int lj = ej - sj + 1;
    int lk = ek - sk + 1;

    //cout << "scoring matrix size : [" << lj << "][" << lk << "]" << endl;

    int mid = (ei+si)/2;
    int score1[2][lj][lk];
    int trace[2][lj][lk];
    for(int i = 0; i <= 1; i++)
        for(int j = 0; j < lj; j++)
            for(int k = 0; k < lk; k++)
            {
                score1[i][j][k] = -10000;
                trace[i][j][k] = -1;
            }


    score1[0][0][0] = 0;
    trace[0][0][0] = -1; // start point
    for(int k = 1; k < lk; k++)
    {
        score1[0][0][k] = scoring[4][4][s3[sk+k-1]] * k;
        trace[0][0][k] = 6;
    }
    for(int j = 1; j < lj; j++)
    {
        score1[0][j][0] = scoring[4][s2[sj+j-1]][4] * j;
        trace[0][j][0] = 5;
    }

    for(int j = 1; j < lj; j++)
        for(int k = 1; k < lk; k++)
        {
            int mm = 0;
            int ss[3] = {-10000, -10000, -10000};
            int dir = 0;

            // - j k
            ss[0] = score1[0][j-1][k-1] + scoring[4][s2[sj+j-1]][s3[sk+k-1]];
            // - j -
            ss[1] = score1[0][j-1][k] + scoring[4][s2[sj+j-1]][4];
            // - - k
            ss[2] = score1[0][j][k-1] + scoring[4][4][s3[sk+k-1]];

            mm = ss[0];
            dir = 4;
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                    if(i == 1)
                        dir = 5;
                    if(i == 2)
                        dir = 6;
                }
            }
            score1[0][j][k] = mm;
            trace[0][j][k] = dir;
        }

    int ly = 0;
    for(int i = 1; i < li; i++)  // mid
    //for(int i = 1; i <= li; i++)
    {
        ly = (ly+1)%2;
        cout << ly << endl;

        score1[ly][0][0] = score1[(ly+1)%2][0][0] + scoring[s1[si+i-1]][4][4];
        trace[ly][0][0] = 3;

        for(int t = 1; t < lj; t++)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            int dir = 0;
            // i j -
            ss[0] = score1[(ly+1)%2][t-1][0] + scoring[s1[si+i-1]][s2[sj+t-1]][4];
            // i - -
            ss[1] = score1[(ly+1)%2][t][0] + scoring[s1[si+i-1]][4][4];
            // - j -
            ss[2] = score1[ly][t-1][0] + scoring[4][s2[sj+t-1]][4];

            mm = ss[0];
            dir = 1;
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                    if(i == 1)
                        dir = 3;
                    if(i == 2)
                        dir = 5;
                }
            }
            score1[ly][t][0] = mm;
            trace[ly][t][0] = dir;
        }
        for(int t = 1; t < lk; t++)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            int dir = 0;
            // i - k
            ss[0] = score1[(ly+1)%2][0][t-1] + scoring[s1[si+i-1]][4][s2[sk+t-1]];
            // i - -
            ss[1] = score1[(ly+1)%2][0][t] + scoring[s1[si+i-1]][4][4];
            // - - k
            ss[2] = score1[ly][0][t-1] + scoring[4][4][s2[sk+t-1]];

            mm = ss[0];
            dir = 2;
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                    if(i == 1)
                        dir = 3;
                    if(i == 2)
                        dir = 6;
                }
            }
            score1[ly][0][t] = mm;
            trace[ly][0][t] = dir;
        }

        for(int j = 1; j < lj; j++)
        {
            for(int k = 1; k < lk; k++)
            {
                int ss[7] = {-10000, -10000, -10000, -10000, -10000, -10000, -10000};
                int mm = 0;
                int ply = (ly+1)%2;
                int dir = 0;

                // i j k
                ss[0] = score1[ply][j-1][k-1] + scoring[s1[si+i-1]][s2[sj+j-1]][s3[sk+k-1]];
                // i j -
                ss[1] = score1[ply][j-1][k] + scoring[s1[si+i-1]][s2[sj+j-1]][4];
                // i - k
                ss[2] = score1[ply][j][k-1] + scoring[s1[si+i-1]][4][s3[sk+k-1]];
                // - j k
                ss[4] = score1[ly][j-1][k-1] + scoring[4][s2[sj+j-1]][s3[sk+k-1]];
                // i - -
                ss[3] = score1[ply][j][k] + scoring[s1[si+i-1]][4][4];
                // - j -
                ss[5] = score1[ly][j-1][k] + scoring[4][s2[sj+j-1]][4];
                // - - k
                ss[6] = score1[ly][j][k-1] + scoring[4][4][s3[sk+k-1]];

                mm = ss[0];
                dir = 0;
                for(int i = 1; i < 7; i++)
                {
                    if(ss[i] > mm)
                    {
                        mm = ss[i];
                        dir = i;
                    }
                }
                score1[ly][j][k] = mm;
                trace[ly][j][k] = dir;
            }
        }
    }

    for(int i = 0; i < lj; i++)
    {
        for(int j = 0; j < lk; j++)
            cout << trace[0][i][j] << "\t";
        cout << endl;
    }

        for(int i = 0; i < lj; i++)
    {
    for(int j = 0; j < lk; j++)
            cout << trace[1][i][j] << "\t";
        cout << endl;
    }

    // find the path from (ei, ej, ek) to (si, sj, sk)
    int j = lj-1;
    int k = lk-1;
    int i = 1;
    int p = 0;

    while(1)
    {
        cout << "i = " << i << " j = " << j << " k = " << k << endl;
        p = trace[i][j][k];
        if(p < 0)
            return;


        switch(p)
        {
            case 0: // perfect matching
                pm++;
                i--;
                j--;
                k--;
                break;
            case 1:
                i--;
                j--;
                break;
            case 2:
                i--;
                k--;
                break;
            case 3:
                i--;
                break;
            case 4:
                j--;
                k--;
                break;
            case 5:
                j--;
                break;
            case 6:
                k--;
                break;
        }
        if(i<0 || j<0 || k<0)
            return;
        length++;
    }

    return;
}


void alignment(int si, int sj, int sk, int ei, int ej, int ek)
{
    int li = ei - si + 1;
    int lj = ej - sj + 1;
    int lk = ek - sk + 1;

    if(ei - si <= 1)
        return;

    //cout << "scoring matrix size : [" << lj << "][" << lk << "]" << endl;

    int mid = (ei+si)/2;
    int score1[2][lj+1][lk+1];
    int score2[2][lj+1][lk+1];
    for(int i = 0; i <= 1; i++)
        for(int j = 0; j <= lj; j++)
            for(int k = 0; k <= lk; k++)
            {
                score1[i][j][k] = -10000;
                score2[i][j][k] = -10000;
            }


    score1[0][0][0] = 0;
    for(int k = 1; k <= lk; k++)
    {
        score1[0][0][k] = scoring[4][4][s3[sk+k-1]] * k;
    }
    for(int j = 1; j <= lj; j++)
        score1[0][j][0] = scoring[4][s2[sj+j-1]][4] * j;

    for(int j = 1; j <= lj; j++)
        for(int k = 1; k <= lk; k++)
        {
            int mm = 0;
            int ss[3] = {-10000, -10000, -10000};
            // - j k
            ss[0] = score1[0][j-1][k-1] + scoring[4][s2[sj+j-1]][s3[sk+k-1]];
            // - j -
            ss[1] = score1[0][j-1][k] + scoring[4][s2[sj+j-1]][4];
            // - - k
            ss[2] = score1[0][j][k-1] + scoring[4][4][s3[sk+k-1]];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score1[0][j][k] = mm;
        }

    int ly = 0;
    for(int i = 1; i <= mid-si; i++)  // mid
    //for(int i = 1; i <= li; i++)
    {
        ly = (ly+1)%2;

        score1[ly][0][0] = score1[(ly+1)%2][0][0] + scoring[s1[si+i-1]][4][4];

        for(int t = 1; t <= lj; t++)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            // i j -
            ss[0] = score1[(ly+1)%2][t-1][0] + scoring[s1[si+i-1]][s2[sj+t-1]][4];
            // i - -
            ss[1] = score1[(ly+1)%2][t][0] + scoring[s1[si+i-1]][4][4];
            // - j -
            ss[2] = score1[ly][t-1][0] + scoring[4][s2[sj+t-1]][4];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score1[ly][t][0] = mm;
        }
        for(int t = 1; t <= lk; t++)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            // i - k
            ss[0] = score1[(ly+1)%2][0][t-1] + scoring[s1[si+i-1]][4][s2[sk+t-1]];
            // i - -
            ss[1] = score1[(ly+1)%2][0][t] + scoring[s1[si+i-1]][4][4];
            // - - k
            ss[2] = score1[ly][0][t-1] + scoring[4][4][s2[sk+t-1]];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score1[ly][0][t] = mm;
        }

        for(int j = 1; j <= lj; j++)
        {
            for(int k = 1; k <= lk; k++)
            {
                int ss[7] = {-10000, -10000, -10000, -10000, -10000, -10000, -10000};
                int mm = 0;
                int ply = (ly+1)%2;

                // i j k
                ss[0] = score1[ply][j-1][k-1] + scoring[s1[si+i-1]][s2[sj+j-1]][s3[sk+k-1]];
                // i j -
                ss[1] = score1[ply][j-1][k] + scoring[s1[si+i-1]][s2[sj+j-1]][4];
                // i - k
                ss[2] = score1[ply][j][k-1] + scoring[s1[si+i-1]][4][s3[sk+k-1]];
                // - j k
                ss[3] = score1[ly][j-1][k-1] + scoring[4][s2[sj+j-1]][s3[sk+k-1]];
                // i - -
                ss[4] = score1[ply][j][k] + scoring[s1[si+i-1]][4][4];
                // - j -
                ss[5] = score1[ly][j-1][k] + scoring[4][s2[sj+j-1]][4];
                // - - k
                ss[6] = score1[ly][j][k-1] + scoring[4][4][s3[sk+k-1]];
                mm = ss[0];
                for(int i = 1; i < 7; i++)
                {
                    if(ss[i] > mm)
                    {
                        mm = ss[i];
                    }
                }
                score1[ly][j][k] = mm;
            }
        }
    }
    int fwd = ly;
    //cout << "forward score: " << score1[fwd][lj][lk] << endl;
    /*
    for(int i = 0; i <= lj; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score1[fwd][i][j] << "\t";
        cout << endl;
    }
    */
    // -------------------------------------- backward --------------------------------//

    score2[0][0][0] = 0;
    for(int k = lk-1; k >= 0; k--)
        score2[0][lj][k] = scoring[4][4][s3[sk+k]]*(lk-k);
    for(int j = lj-1; j >= 0; j--)
        score2[0][j][lk] = scoring[4][s2[sj+j]][4]*(lj-j);

    for(int j = lj-1; j >= 0; j--)
        for(int k = lk-1; k >= 0; k--)
        {

            int mm = 0;
            int ss[3] = {-10000, -10000, -10000};
            // - j k
            ss[0] = score2[0][j+1][k+1] + scoring[4][s2[sj+j]][s3[sk+k]];
            // - j -
            ss[1] = score2[0][j+1][k] + scoring[4][s2[sj+j]][4];
            // - - k
            ss[2] = score2[0][j][k+1] + scoring[4][4][s3[sk+k]];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score2[0][j][k] = mm;
        }


    ly = 0;
    for(int i = li-1; i >= mid; i--)  // mid
    //for(int i = li-1; i >= 0; i--)
    {
        ly = (ly+1)%2;

        score2[ly][lj][lk] = score2[(ly+1)%2][lj][lk] + scoring[s1[si+i]][4][4];

        for(int t = lj-1; t >= 0; t--)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            // i j -
            ss[0] = score2[(ly+1)%2][t+1][lk] + scoring[s1[si+i]][s2[sj+t]][4];
            // i - -
            ss[1] = score2[(ly+1)%2][t][lk] + scoring[s1[si+i]][4][4];
            // - j -
            ss[2] = score2[ly][t+1][lk] + scoring[4][s2[sj+t]][4];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score2[ly][t][lk] = mm;
        }
        for(int t = lk-1; t >= 0; t--)
        {
            int ss[3] = {-10000, -10000, -10000};
            int mm = 0;
            // i - k
            ss[0] = score2[(ly+1)%2][lj][t+1] + scoring[s1[si+i]][4][s2[sk+t-1]];
            // i - -
            ss[1] = score2[(ly+1)%2][lj][t] + scoring[s1[si+i]][4][4];
            // - - k
            ss[2] = score2[ly][lj][t+1] + scoring[4][4][s2[sk+t]];

            mm = ss[0];
            for(int i = 1; i < 3; i++)
            {
                if(ss[i] > mm)
                {
                    mm = ss[i];
                }
            }
            score2[ly][lj][t] = mm;
        }

        for(int j = lj-1; j >= 0; j--)
        {
            for(int k = lk-1; k >= 0; k--)
            {
                int ss[7] = {-10000, -10000, -10000, -10000, -10000, -10000, -10000};
                int mm = 0;
                int ply = (ly+1)%2;

                // i j k
                ss[0] = score2[ply][j+1][k+1] + scoring[s1[si+i]][s2[sj+j]][s3[sk+k]];
                // i j -
                ss[1] = score2[ply][j+1][k] + scoring[s1[si+i]][s2[sj+j]][4];
                // i - k
                ss[2] = score2[ply][j][k+1] + scoring[s1[si+i]][4][s3[sk+k]];
                // - j k
                ss[3] = score2[ly][j+1][k+1] + scoring[4][s2[sj+j]][s3[sk+k]];
                // i - -
                ss[4] = score2[ply][j][k] + scoring[s1[si+i]][4][4];
                // - j -
                ss[5] = score2[ly][j+1][k] + scoring[4][s2[sj+j]][4];
                // - - k
                ss[6] = score2[ly][j][k+1] + scoring[4][4][s3[sk+k]];
                mm = ss[0];
                for(int i = 1; i < 7; i++)
                {
                    if(ss[i] > mm)
                    {
                        mm = ss[i];
                    }
                }
                score2[ly][j][k] = mm;
            }
        }
    }
    int bwd = ly;
    /*
    for(int i = 0; i <= lj; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score2[bwd][i][j] << "\t";
        cout << endl;
    }
    */
    //cout << "backward score: " << score2[bwd][0][0] << endl;

    // find the max sum of two score layer in mid
    int max = score1[fwd][0][0] + score2[bwd][0][0];
    int mx = mid;
    int my = sj;
    int mz = sk;
    for(int j = 0; j <= lj; j++)
        for(int k = 0; k <= lk; k++)
        {
            if(score1[fwd][j][k]+score2[bwd][j][k] > max)
            {
                max = score1[fwd][j][k]+score2[bwd][j][k];
                my = sj + j;
                mz = sk + k;
            }
        }

    point3d pp;
    pp.x = mx;
    pp.y = my;
    pp.z = mz;
    midpoint.push_back(pp);
    //cout << "Find a mid point (" << mx << ", " << my << ", " << mz << ")" << endl;
    //printf("(%d, %d, %d) -> (%d, %d, %d) = (%d, %d, %d)\n", si, sj, sk, ei, ej, ek, pp.x, pp.y, pp.z);

    if(maxvalue)
    {
        maxvalue = false;
        cout << "Final score: " << max << endl;
    }

    alignment(si, sj, sk, mx, my, mz);
    alignment(mx, my, mz, ei, ej, ek);

    return;
}

void createscoring()
{
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            for(int k = 0; k < 5; k++)
            {
                int s = 0;
                if(i == 4 && j == 4 && k == 4)
                    s = 0;
                else if((i == 4 && j == 4 && k != 4) || (i == 4 && j != 4 && k == 4) || (i != 4 && j == 4 && k == 4))
                {
                    s = INDEL*2;
                }
                else if(i == 4 && j != 4 && k != 4)
                {
                    if(j == k)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                else if(i != 4 && j == 4 && k != 4)
                {
                    if( i == k)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                else if(i != 4 && j != 4 && k == 4)
                {
                    if( i == j)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                else if(i < 4 && j < 4 && k < 4)
                {
                    if(i == j && i == k)
                        s = MATCH*3;
                    else if((i == j && i != k) || (i != j && i == k) || (j == k && j != i) || (j != k && j == i))
                        s = MATCH + 2*MISMATCH;
                    else if(i != j && j != k && k != i)
                        s = MISMATCH*3;
                }
                scoring[i][j][k] = s;
            }
        }
    }
}
