#include "boost/multi_array.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>

using namespace std;
const int DEBUG = 0;

typedef boost::multi_array<int, 3> score_type;

// test score
const int MATCH = 5;
const int MISMATCH = -4;
const int INDEL = -8;

// from ATGC- to 0 1 2 3 4
int n2n(char c)
{
    if(c == 'A')
        return 0;
    if(c == 'T')
        return 1;
    if(c == 'G')
        return 2;
    if(c == 'C')
        return 3;
    // indel
    return 4;
}


// global vector
vector<char> s1;
vector<char> s2;
vector<char> s3;

int scoring[5][5][5];

struct point3d{
    int m;
    int n;
    int k;
};
vector<point3d> midpoint;

bool ranking(const point3d &lhs, const point3d &rhs)
{
    return lhs.m < rhs.m;
}

bool maxvalue = true;

void createscoring();
void alignment(int sm, int sn, int sk, int em, int en, int ek);
int backtrace(int sm, int sn, int sk, int em, int en, int ek);

int main(int argc, char **argv)
{
    if(argc != 4)
    {
        cout << "Usage: ./program seq1 seq2 seq3" << endl;
        return -1;
    }

    // create scoring matrix
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
            s1.push_back(line[i]);
        }
    }

    while(getline(fseq2,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            s2.push_back(line[i]);
        }
    }

    while(getline(fseq3,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            s3.push_back(line[i]);
        }
    }

    point3d pp;
    pp.m = 0;
    pp.n = 0;
    pp.k = 0;
    midpoint.push_back(pp);
    // choose the third sequence to cut at middle
    int m = s1.size();
    int n = s2.size();
    int k = s3.size();

    cout << "seq 1 size : " << m << endl;
    cout << "seq 2 size : " << n << endl;
    cout << "seq 3 size : " << k << endl;


    clock_t start = clock();
    // mn space
    alignment(0, 0, 0, m, n, k);
    clock_t end = clock();

    pp.m = m-1;
    pp.n = n-1;
    pp.k = k-1;
    midpoint.push_back(pp);
    sort(midpoint.begin(), midpoint.end(),&ranking);
    if(DEBUG == 1)
    {
        for(int i = 0; i < midpoint.size(); i++)
        {
            pp = midpoint[i];
            printf("Find the max crossing point (%d, %d, %d)\n", pp.m, pp.n, pp.k);
        }
    }

    int Length = 0;
    for(int i = 0; i < midpoint.size()-1; i++)
    {
        point3d pp = midpoint[i];
        point3d qq = midpoint[i+1];
    // only one path (i, 0, 0) to (i+1, 0, 0)
    /*
        if(ln == 0  && lk == 0)
        return 1;
    if(lm == 0  && lk == 0)
        return 1;
    if(ln == 0  && lm == 0)
        return 1;
    
    if(ln == 0 && lk == 0 && lm == 0)
        return 0;
    */
        
        printf("Find the trace from point (%d, %d, %d) to (%d, %d, %d)\n", pp.m, pp.n, pp.k, qq.m, qq.n, qq.k);

        if(pp.m == qq.m && pp.n == qq.n && pp.k == qq.k)
        {
            continue;
        }
        else if((pp.m==qq.m && pp.n==qq.n) || (pp.m==qq.m && pp.k==qq.k) || (pp.n==qq.n && pp.k==qq.k))
        {
            Length++;
        }
        else
            Length += backtrace(pp.m, pp.n, pp.k, qq.m, qq.n, qq.k);
        //cout << "Length = " << backtrace(pp.m, pp.n, pp.k, qq.m, qq.n, qq.k) << endl;;
        cout << " Length = " << Length << endl;
    }

    cout << "Alignment length = : " << Length << endl;
    cout << "Execution time: " << (double)(end - start)/CLOCKS_PER_SEC << endl;

    fseq1.close();
    fseq2.close();
    fseq3.close();

    return 0;
}

// return the length of the trace, s2[j] and s3[k]
int backtrace(int sm, int sn, int sk, int em, int en, int ek)
{
    int length = 1;

    int lm = em - sm;
    int ln = en - sn;
    int lk = ek - sk;

    //cout << "lm = " << lm << endl;
    //cout << "ln = " << ln << endl;
    //cout << "lk = " << lk << endl;


    //score_type score(boost::extents[2][ln+1][lk+1]);
    //score_type trace(boost::extents[2][ln+1][lk+1]);
    int score[2][ln+1][lk+1]; // 1 j k
    int trace[2][ln+1][lk+1]; // 1 j k
    int ly = 0;

    score[0][0][0] = 0;
    // the initial layer score
    for(int i = 1; i <= lk; i++)
    {
        score[0][0][i] = scoring[n2n('S')][n2n('S')][n2n(s3[sk+i-1])]*i;
        trace[0][0][i] = 3;
    }
    for(int j = 1; j <= ln; j++)
    {
        score[0][j][0] = scoring[n2n('S')][n2n(s2[sn+j-1])][n2n('S')]*j;
        trace[0][j][0] = 2;
    }
    for(int i = 1; i <= ln; i++)
        for(int j = 1; j <=lk; j++)
        {
            int ss = 0;
            int ss1 = 0;
            int ind = 6;
            // - j k
            ss = score[0][i-1][j-1] + scoring[n2n('S')][n2n(s2[sn+i-1])][n2n(s3[sk+j-1])];
            // - j -
            ss1 = score[0][i-1][j] + scoring[n2n('S')][n2n(s2[sn+i-1])][n2n(s3[sk+j])];
            if(ss < ss1)
                ind = 2;
            ss = ss > ss1? ss : ss1;
            // - - k
            ss1 = score[0][i][j-1] + scoring[n2n('S')][n2n(s2[sn+i])][n2n(s3[sk+j-1])];
            if(ss < ss1)
                ind = 3;
            ss = ss > ss1? ss : ss1;
            score[0][i][j] = ss;
            trace[0][i][j] = ind;
        }
    /*
    cout << "first plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score[0][i][j] << '\t';
        cout << endl;
    }
    cout << "---" << endl;
    */
    for(int i = sm; i < em; i++)
    {
        cout << "i = " << i << endl;
        ly++;
        ly = ly%2;

        // corner 
        score[ly][0][0] = score[(ly+1)%2][0][0] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
        trace[1][0][0] = 1;

        for(int t = sk; t < ek; t++)
        {
            int ss = 0;
            int ss1 = 0;

            int ind = 3;
            // - - k
            ss = score[ly][0][t-sk] + scoring[n2n('S')][n2n('S')][n2n(s3[t])];
            // i - -
            ss1 = score[(ly+1)%2][0][t-sk+1] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            if(ss < ss1)
                ind = 1;
            // i - k
            ss1 = score[(ly+1)%2][0][t-sk] + scoring[n2n(s1[i])][n2n('S')][n2n(s3[t])];
            if(ss < ss1)
                ind = 5;
            ss = ss > ss1? ss : ss1;

            score[ly][0][t-sk+1] = ss; 
            trace[1][0][t-sk+1] = ind;
        }
        for(int t = sn; t < en; t++)
        {
            int ss = 0;
            int ss1 = 0;

            int ind = 2;

            // - j - 
            ss = score[ly][t-sn][0] + scoring[n2n('S')][n2n(s2[t])][n2n('S')];
            // i - -
            ss1 = score[(ly+1)%2][t-sn+1][0] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            if(ss < ss1)
                ind = 1;
            ss = ss > ss1? ss : ss1;
            // i j -
            ss1 = score[(ly+1)%2][t-sn][0] + scoring[n2n(s1[i])][n2n(s2[t])][n2n('S')];
            if(ss <ss1)
                ind = 4;
            ss = ss > ss1? ss : ss1;

            score[ly][t-sn+1][0] = ss;
            trace[1][t-sn+1][0] = ind;
        }


        for(int j = sn; j < en; j++)
        {
            //cout << " j = " << j << endl;;
            for(int k = sk; k < ek; k++)
            {
                //cout << " k = " << k << endl;;
                int ss[7] = {-100000, -100000, -100000, -100000, -100000, -100000, -100000};

                int ly2 = (ly+1)%2;
                // no indel
                ss[0] = score[ly2][j-sn][k-sk] + scoring[n2n(s1[i+1])][n2n(s2[j+1])][n2n(s3[k+1])];
                // two indel
                // i - -
                ss[1] = score[ly2][j-sn+1][k-sk+1] + scoring[n2n(s1[i+1])][n2n('S')][n2n('S')];
                // - j -
                ss[2] = score[ly][j-sn][k-sk+1] + scoring[n2n('S')][n2n(s2[j+1])][n2n('S')];
                // - - k
                ss[3] = score[ly][j-sn+1][k-sk] + scoring[n2n('S')][n2n('S')][n2n(s3[k+1])];
                // one indel
                // i j -
                ss[4] = score[ly2][j-sn][k-sk+1] + scoring[n2n(s1[i+1])][n2n(s2[j+1])][n2n('S')];
                // i - k
                ss[5] = score[ly2][j-sn+1][k-sk] + scoring[n2n(s1[i+1])][n2n('S')][n2n(s3[k+1])];
                // - j k
                ss[6] = score[ly][j-sn][k-sk] + scoring[n2n('S')][n2n(s2[j+1])][n2n(s3[k+1])];

                int ind = 0;
                for(int t = 0; t < 7; t++)
                {
                    //cout << " ss[" << t << "] = " << ss[t] << endl; 
                    if(ss[t] >= ss[0])
                    {
                        ind = t;
                        ss[0] = ss[t];
                    }
                }
                trace[1][j-sn][k-sk] = ind;
                score[ly][j-sn+1][k-sk+1] = ss[0];
    //            cout << " ind = " << trace[1][j-sn][k-sk] << " ss = " << score[ly][j-sn+1][k-sk+1] << endl;
            }
        }
        /*
        for(int w = 0; w <= ln; w++)
        {
            for(int u = 0; u <= lk; u++)
                cout << score[ly][w][u] << '\t';
            cout << endl;
        }
        cout << "---" << endl;
        */
    }

    // find the backtrace
    trace[1][0][0] = 1;
    trace[0][0][0] = -1;
   /*
    for(int k = 0; k <= lm; k++)
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << trace[k][i][j] << '\t';
        cout << endl;
    }
    cout << "---" << endl;
   */

    ly = 1;
    int i = lm;
    int j = lk;
    while(1)
    {
        if(trace[ly][i][j] < 0)
            break;

        if(i<0 || j <0 || ly < 0)
            break;

        length++;
        if(trace[ly][i][j] == 0)
        {
            ly--;
            i--;
            j--;
        }
        else if(trace[ly][i][j] == 1)
        {
            ly--;
        }
        else if(trace[ly][i][j] == 2)
        {
            j--;
        }
        else if(trace[ly][i][j] == 3)
        {
            i--;
        }
        else if(trace[ly][i][j] == 4)
        {
            ly--;
            j--;
        }
        else if(trace[ly][i][j] == 5)
        {
            ly--;
            i--;
        }
        else if(trace[ly][i][j] == 6)
        {
            i--;
            j--;
        }
    }

    return length;
}

void alignment(int sm, int sn, int sk, int em, int en, int ek)
{
    int lm = em - sm;
    int ln = en - sn;
    int lk = ek - sk;
    if(lm <= 1)
        return;
    // choose the first sequence to cut at the mid point
    int mid = sm + lm/2;

    //cout << "Search from: (" << sm << ", " << sn << ", " << sk << ") to (" << em << ", " << en << ", " << ek << ")" << endl;
    // ------------------------------------ forward --------------------------------------//
    /* 
    cout << " sm = " << sm << endl;
    cout << " sn = " << sn << endl;
    cout << " sk = " << sk << endl;
    cout << " em = " << em << endl;
    cout << " en = " << en << endl;
    cout << " ek = " << ek << endl;
    cout << " lm = " << lm << endl;
    cout << " ln = " << ln << endl;
    cout << " lk = " << lk << endl;
    cout << " mid = " << mid << endl;
    */

    //score_type score(boost::extents[2][ln+1][lk+1]);
    int score[2][ln+1][lk+1];
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < ln+1; j++)
            for(int k = 0; k < lk+1; k++)
                score[i][j][k] = 0;

    // the initial layer score
    for(int i = 1; i <= lk; i++)
        score[0][0][i] = scoring[n2n('S')][n2n('S')][n2n(s3[sk+i-1])]*i;
    for(int j = 1; j <= ln; j++)
        score[0][j][0] = scoring[n2n('S')][n2n(s2[sn+j-1])][n2n('S')]*j;
    for(int i = 1; i <= ln; i++)
        for(int j = 1; j <=lk; j++)
        {
            int ss = 0;
            int ss1 = 0;
            // - j k
            ss = score[0][i-1][j-1] + scoring[n2n('S')][n2n(s2[sn+i-1])][n2n(s3[sk+j-1])];
            // - j -
            ss1 = score[0][i-1][j] + scoring[n2n('S')][n2n(s2[sn+i-1])][n2n(s3[sk+j])];
            ss = ss > ss1? ss : ss1;
            // - - k
            ss1 = score[0][i][j-1] + scoring[n2n('S')][n2n(s2[sn+i])][n2n(s3[sk+j-1])];
            ss = ss > ss1? ss : ss1;
            score[0][i][j] = ss;
        }
    /*
    cout << "first plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score[0][i][j] << '\t';
        cout << endl;
    }
    cout << "---" << endl;
    */
     
    int ly = 0;
    for(int i = sm+1; i <= mid; i++)  // the upper half
    {
    
        ly++;  // start 1 -> 0 -> 1 ...
        ly = ly % 2;
        // corner 
        score[ly][0][0] = score[(ly+1)%2][0][0] + scoring[n2n(s1[i-1])][n2n('S')][n2n('S')];
        for(int t = sk+1; t <= ek; t++)
        {
            int ss = 0;
            int ss1 = 0;

            // - - k
            ss = score[ly][0][t-sk-1] + scoring[n2n('S')][n2n('S')][n2n(s3[t-1])];
            // i - -
            ss1 = score[(ly+1)%2][0][t-sk] + scoring[n2n(s1[i-1])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            // i - k
            ss1 = score[(ly+1)%2][0][t-sk-1] + scoring[n2n(s1[i-1])][n2n('S')][n2n(s3[t-1])];
            ss = ss > ss1? ss : ss1;

            score[ly][0][t-sk] = ss; 
        }
        for(int t = sn+1; t <= en; t++)
        {
            int ss = 0;
            int ss1 = 0;

            // - j - 
            ss = score[ly][t-sn-1][0] + scoring[n2n('S')][n2n(s2[t-1])][n2n('S')];
            // i - -
            ss1 = score[(ly+1)%2][t-sn][0] + scoring[n2n(s1[i-1])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            // i j -
            ss1 = score[(ly+1)%2][t-sn-1][0] + scoring[n2n(s1[i-1])][n2n(s2[t-1])][n2n('S')];
            ss = ss > ss1? ss : ss1;
            score[ly][t-sn][0] = ss;
        }
        for(int j = sn+1; j <= en; j++)
        {
            for(int k = sk+1; k <= ek; k++)
            {
                int ss = 0;
                int ss1 = 0;

                int ly2 = (ly+1)%2;
                // no indel
                ss = score[ly2][j-sn-1][k-sk-1] + scoring[n2n(s1[i-1])][n2n(s2[j-1])][n2n(s3[k-1])];
                // two indel
                // i - -
                ss1 = score[ly2][j-sn][k-sk] + scoring[n2n(s1[i-1])][n2n('S')][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - j -
                ss1 = score[ly][j-sn-1][k-sk] + scoring[n2n('S')][n2n(s2[j-1])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - - k
                ss1 = score[ly][j-sn][k-sk-1] + scoring[n2n('S')][n2n('S')][n2n(s3[k-1])];
                ss = ss > ss1? ss : ss1;
                // one indel
                // i j -
                ss1 = score[ly2][j-sn-1][k-sk] + scoring[n2n(s1[i-1])][n2n(s2[j-1])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // i - k
                ss1 = score[ly2][j-sn][k-sk-1] + scoring[n2n(s1[i-1])][n2n('S')][n2n(s3[k-1])];
                ss = ss > ss1? ss : ss1;
                // - j k
                ss1 = score[ly][j-sn-1][k-sk-1] + scoring[n2n('S')][n2n(s2[j-1])][n2n(s3[k-1])];
                ss = ss > ss1? ss : ss1;

                score[ly][j-sn][k-sk] = ss;
            }
        }
        /*
        for(int w = 0; w <= ln; w++)
        {
            for(int u = 0; u <= lk; u++)
                cout << score[ly][w][u] << '\t';
            cout << endl;
        }
        cout << "---" << endl;
        */
    }

    // ------------------------------------ forward --------------------------------------//

    int fwd = ly;

    /*
    cout << "fwd mid plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score[ly][i][j] << '\t';
        cout << endl;
    }
    
    */
    // -------------------------------------------------------------------------------------//
    // -------------------------------------------------------------------------------------//
    
    // --------------------------------------backward-----------------------------------------//
     
    int score2[2][ln+1][lk+1];
    //score_type score2(boost::extents[2][ln+1][lk+1]);
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < ln+1; j++)
            for(int k = 0; k < lk+1; k++)
                score2[i][j][k] = 0;

    // the initial layer score
    for(int i = ek; i > sk; i--)
    {
        score2[0][ln][i-sk-1] = scoring[n2n('S')][n2n('S')][n2n(s3[i-1])]*(ek-i+1);
    }

    for(int j = en; j > sn; j--)
        score2[0][j-sn-1][lk] = scoring[n2n('S')][n2n(s2[j-1])][n2n('S')]*(en-j+1);

    for(int i = en; i > sn; i--)
        for(int j = ek; j > sk; j--)
        {
            int ss = 0;
            int ss1 = 0;
            // - i j
            ss = score2[0][i-sn][j-sk] + scoring[n2n('S')][n2n(s2[i-1])][n2n(s3[j-1])];
            // - - k
            ss1 = score2[0][i-sn][j-sk-1] + scoring[n2n('S')][n2n(s2[i])][n2n(s3[j-1])];
            ss = ss > ss1? ss : ss1;
            // - j -
            ss1 = score2[0][i-sn-1][j-sk] + scoring[n2n('S')][n2n(s2[i-1])][n2n(s3[j])];
            ss = ss > ss1? ss : ss1;
            
            score2[0][i-sn-1][j-sk-1] = ss;
        }
    /*
    cout << "last plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score2[0][i][j] << '\t';
        cout << endl;
    }
    cout << "---" << endl;
    */
    ly = 0;
    for(int i = em-1; i >= mid; i--)  // the backward half
    {
    
        ly++;  // start 1 -> 0 -> 1 ...
        ly = ly % 2;

        score2[ly][ln][lk] = score2[(ly+1)%2][ln][lk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];

        for(int t = ek-1; t >= sk; t--)
        {
            int ss = 0;
            int ss1 = 0;

            //- - k
            ss = score2[ly][ln][t-sk+1] + scoring[n2n('S')][n2n('S')][n2n(s3[t])];
            //i - -
            ss1 = score2[(ly+1)%2][ln][t-sk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            //i - k
            ss1 = score2[(ly+1)%2][ln][t-sk+1] + scoring[n2n(s1[i])][n2n('S')][n2n(s3[t])];
            ss = ss > ss1? ss : ss1;

            score2[ly][ln][t-sk] = ss;
        }
        for(int t = en-1; t >= sn; t--)
        {
            int ss = 0;
            int ss1 = 0;

            // - j -
            ss = score2[ly][t-sn+1][lk] + scoring[n2n('S')][n2n(s2[t])][n2n('S')];
            // i - -
            ss1 = score2[(ly+1)%2][t-sn][lk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            // i j -
            ss1 = score2[(ly+1)%2][t-sn+1][lk] + scoring[n2n(s1[i])][n2n(s2[t])][n2n('S')];
            ss = ss > ss1? ss : ss1;
            score2[ly][t-sn][lk] = ss;
        }

        for(int j = en-1; j >= sn; j--)
        {
            for(int k = ek-1; k >= sk; k--)
            {
                int ss = 0;
                int ss1 = 0;

                int ly2 = (ly+1)%2;
                // no indel
                ss = score2[ly2][j-sn+1][k-sk+1] + scoring[n2n(s1[i])][n2n(s2[j])][n2n(s3[k])];
                // two indel
                // i - -
                ss1 = score2[ly2][j-sn][k-sk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - j -
                ss1 = score2[ly][j-sn+1][k-sk] + scoring[n2n('S')][n2n(s2[j])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - - k
                ss1 = score2[ly][j-sn][k-sk+1] + scoring[n2n('S')][n2n('S')][n2n(s2[k])];
                ss = ss > ss1? ss : ss1;
                // one indel
                // i j -
                ss1 = score2[ly2][j-sn+1][k-sk] + scoring[n2n(s1[i])][n2n(s2[j])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // i - k
                ss1 = score2[ly2][j-sn][k-sk+1] + scoring[n2n(s1[i])][n2n('S')][n2n(s3[k])];
                ss = ss > ss1? ss : ss1;
                // - j k
                ss1 = score2[ly][j-sn+1][k-sk+1] + scoring[n2n('S')][n2n(s2[j])][n2n(s3[k])];
                ss = ss > ss1? ss : ss1;

                score2[ly][j-sn][k-sk] = ss;
                //cout << ss << endl;
            }
        }
        /*
        for(int w = 0; w <= ln; w++)
        {
            for(int u = 0; u <= lk; u++)
                cout << score2[ly][w][u] << '\t';
            cout << endl;
        }

        cout << "---" << endl;
        */
        
    }
    
    int bwd = ly;
    // ------------------------------------ backward --------------------------------------//
    /*
    cout << "bwd mid plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score2[ly][i][j] << '\t';
        cout << endl;
    }
    */
    // find the max value in the mid layer
    int max = score[fwd][0][0] + score2[bwd][0][0];
    int mi = sn;
    int mj = sk;

    for(int i = 0; i <= ln; i++)
        for(int j = 0; j <= lk; j++)
        {
            if(score[fwd][i][j] + score2[bwd][i][j] > max)
            {
                mi = i + sn;
                mj = j + sk;
                max = score[fwd][i][j] + score2[bwd][i][j];
                //cout << max << endl;
            }
        }

    if(maxvalue)
    {
        maxvalue = false;
        cout << "Max value: " << max << endl;
    }
    //printf("Find the max crossing point (%d, %d, %d) with value : %d\n", mid, mi, mj, max);
    point3d pp;
    pp.m = mid-1;
    pp.n = mi-1;
    pp.k = mj-1;
    midpoint.push_back(pp);

    alignment(sm, sn, sk, mid, mi, mj);
    alignment(mid, mi, mj, em, en, ek);
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
                if(i == 4 && j == 4 && k == 4) // 3 indel, no use
                {
                    s = 0;
                }
                // two indel
                else if((i == 4 && j == 4 && k != 4) || (i == 4 && j != 4 && k == 4) || (i != 4 && j == 4 && k == 4))
                {
                    s = INDEL*2;
                }
                // one indel
                else if((i == 4 && j != 4 && k != 4))
                {
                    if(j == k)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                else if ((i != 4 && j == 4 && k !=4))
                {
                    if(i == k)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                else if ((i != 4 && j != 4 && k ==4))
                {
                    if( i == j)
                        s = INDEL*2 + MATCH;
                    else
                        s = INDEL*2 + MISMATCH;
                }
                // no indel
                else if(i < 4 && j < 4 && k < 4)
                {
                    // A A A
                    if(i == j && i == k)
                        s = MATCH * 3;
                    // A A T
                    else if ((i == j && i != k) || (i != j && i == k) || (j == k && j != i) || (j != k && j == i))
                        s = MATCH + 2*MISMATCH;
                    // A T G
                    else if (i != j && j != k && k != i)
                        s = MISMATCH * 3;
                }
                
                scoring[i][j][k] = s;
            }
        }
    }
}

