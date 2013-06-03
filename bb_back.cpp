#include "boost/multi_array.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;
const int DEBUG = 0;

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

typedef boost::multi_array<int, 3> score_type;
typedef score_type::index id;
score_type scoring(boost::extents[5][5][5]);

typedef boost::multi_array<int, 2> array_type;
typedef array_type::index ind;

void createscoring();
void alignment_mnk(int sm, int sn, int sk, int lm, int ln, int lk);
void alignment(int sm, int sn, int sk, int m, int n, int k);

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

    // choose the third sequence to cut at middle
    int m = s1.size();
    int n = s2.size();
    int k = s3.size();

    // mnk space
    //alignment_mnk(0, 0, 0, m, n, k);
    // mn space
    alignment(0, 0, 0, m, n, k);

    fseq1.close();
    fseq2.close();
    fseq3.close();

    return 0;
}

// index from sm with size lm, [sm, sm+lm-1]
void alignment(int sm, int sn, int sk, int lm, int ln, int lk)
{
    int em = sm + lm ;
    int en = sn + ln ;
    int ek = sk + lk ;

    // choose the first sequence to cut at the mid point
    int mid = sm + lm/2;

    //cout << "em = " << em << endl;
    //cout << "en = " << en << endl;
    //cout << "ek = " << ek << endl;
    //cout << "mid = " << mid << endl;
    // forward [0, 0, 0] to [mid, en, ek]
    // backward [mid, 0, 0] to [em, en, en]
    //
    // ------------------------------------ forward --------------------------------------//
    score_type score(boost::extents[2][ln+1][lk+1]);

    // the initial layer score
    for(int i = ek; i > sk; i--)
    {
        score[0][ln][i-sk-1] = scoring[n2n('S')][n2n('S')][n2n(s3[i-1])];
    }

    for(int j = en; j > sn; j--)
        score[0][j-sn-1][lk] = scoring[n2n('S')][n2n(s2[j-1])][n2n('S')];
    for(int i = en; i > sn; i--)
        for(int j = ek; j > sk; j--)
            score[0][i-sn-1][j-sk-1] = scoring[n2n('S')][n2n(s2[i-1])][n2n(s3[j-1])];

    // from sm to mid, mid to sm+lm-1
    // [sm, mid]
    int ly = 0;
    for(int i = em-1; i >= sm; i--)  // the upper half
    {
        ly++;  // start 1 -> 0 -> 1 ...
        ly = ly % 2;

        score[ly][ln][lk] = scoring[n2n(s1[i])][n2n('S')][n2n('S')];
//        cout << "layer = " << i << " ly = " << ly << " ln = " << ln << " lk = " << lk << " score = " <<  score[ly][ln][lk] << endl;

        for(int t = ek-1; t >= sk; t--)
        {
            int ss = 0;
            int ss1 = 0;

            //- - k
            ss = score[ly][ln][t-sk+1] + scoring[n2n('S')][n2n('S')][n2n(s3[t])];
            //i - -
            ss1 = score[(ly+1)%2][ln][t-sk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            //i - k
            ss1 = score[(ly+1)%2][ln][t-sk+1] + scoring[n2n(s1[i])][n2n('S')][n2n(s3[t])];
            ss = ss > ss1? ss : ss1;

            score[ly][ln][t-sk] = ss;

            //score[ly][ln][t-sk] = scoring[n2n(s1[i])][n2n('S')][n2n(s3[t])];
        }
        for(int t = en-1; t >= sn; t--)
        {
            int ss = 0;
            int ss1 = 0;

            // - j -
            ss = score[ly][t-sn+1][lk] + scoring[n2n('S')][n2n(s2[t])][n2n('S')];
            // i - -
            ss1 = score[(ly+1)%2][t-sn][lk] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
            ss = ss > ss1? ss : ss1;
            // i j -
            ss1 = score[(ly+1)%2][t-sn+1][lk] + scoring[n2n(s1[i])][n2n(s2[t])][n2n('S')];
            ss = ss > ss1? ss : ss1;
            score[ly][t-sn][lk] = ss;
            //score[ly][t-sn][lk] = scoring[n2n(s1[i])][n2n(s2[t])][n2n('S')];
        }

        for(int j = en-1; j >= sn; j--)
        {
            for(int k = ek-1; k >= sk; k--)
            {
                int ss = 0;
                int ss1 = 0;

                int ly2 = (ly+1)%2;
                // no indel
                ss = score[ly2][j+1][k+1] + scoring[n2n(s1[i])][n2n(s2[j])][n2n(s3[k])];
                // two indel
                // i - -
                ss1 = score[ly2][j][k] + scoring[n2n(s1[i])][n2n('S')][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - j -
                ss1 = score[ly][j+1][k] + scoring[n2n('S')][n2n(s2[j])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // - - k
                ss1 = score[ly][j][k+1] + scoring[n2n('S')][n2n('S')][n2n(s2[k])];
                ss = ss > ss1? ss : ss1;
                // one indel
                // i j -
                ss1 = score[ly2][j+1][k] + scoring[n2n(s1[i])][n2n(s2[j])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                // i - k
                ss1 = score[ly2][j][k+1] + scoring[n2n(s1[i])][n2n('S')][n2n(s3[k])];
                ss = ss > ss1? ss : ss1;
                // - j k
                ss1 = score[ly][j+1][k+1] + scoring[n2n('S')][n2n(s2[j])][n2n(s3[k])];
                ss = ss > ss1? ss : ss1;
                score[ly][j][k] = ss;
                //cout << ss << endl;
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
    /*
    cout << "last plate " << endl;
    for(int i = 0; i <= ln; i++)
    {
        for(int j = 0; j <= lk; j++)
            cout << score[ly][i][j] << '\t';
        cout << endl;
    }
    */
    cout << "final score = " << score[ly][0][0] << endl;
    //cout << "score[1] = " << score[1][0][0] << endl;
}








// O(mnk) space version
void alignment_mnk(int sm, int sn, int sk, int lm, int ln, int lk)
{
    // choose the first sequence to cut at the middle point
    int mid = lm;
    score_type score(boost::extents[lm+1][ln+1][lk+1]);

    // 0 - (mid-1) || mid - (k-1)
    // for the first plate
    
    for(id j = 1; j <= lk; j++)
        score[0][0][j] = scoring[n2n('S')][n2n('S')][n2n(s3[j-1])];
    for(id j = 1; j <= ln; j++)
        score[0][j][0] = scoring[n2n('S')][n2n(s2[j-1])][n2n('S')];
    for(id i = 1; i <= ln; i++)
    {
        for(id j = 1; j <= lk; j++)
            score[0][i][j] = scoring[n2n('S')][n2n(s2[i-1])][n2n(s3[j-1])];
    }

    for(id j = 1; j <= ln; j++)
        score[0][j][0] = scoring[n2n('S')][n2n(s2[j-1])][n2n('S')];
    for(id j = 1; j <= lm; j++)
        score[j][0][0] = scoring[n2n(s1[j-1])][n2n('S')][n2n('S')];
    for(id i = 1; i <= lm; i++)
    {
        for(id j = 1; j <= ln; j++)
        score[i][j][0] = scoring[n2n(s1[i-1])][n2n(s2[j-1])][n2n('S')];
    }

    for(id j = 1; j <= lk; j++)
        score[0][0][j] = scoring[n2n('S')][n2n('S')][n2n(s3[j-1])];
    for(id j = 1; j <= lm; j++)
        score[j][0][0] = scoring[n2n(s1[j-1])][n2n('S')][n2n('S')];
    for(id i = 1; i <= lm; i++)
    {
        for(id j = 1; j <= lk; j++)
        score[i][0][j] = scoring[n2n(s1[i-1])][n2n('S')][n2n(s3[j-1])];
    }

    if(DEBUG == 1)
    {
        for(id i = 0; i <= lm; i++)
        {
            for(id j = 0; j <= ln; j++)
            {
                for(id k = 0; k <= lk; k++)
                    printf("%6d\t", score[i][j][k]);
            //cout << score[i][j] << " ";
            }
                cout << endl;
         }
    }
    
    for(id i = 1; i <= lm; i++)
    {
        //cout << i << endl;
        for(id j = 1; j <= ln; j++)
        {
            for(id k = 1; k <= lk; k++)
            {
                int ss = 0;
                int ss1 = 0;
                int ss2 = 0;
                int ss3 = 0;

                // no indel
                ss = score[i-1][j-1][k-1] + scoring[n2n(s1[i-1])][n2n(s2[j-1])][n2n(s3[k-1])];
                // two indel
                ss1 = score[i-1][j][k] + scoring[n2n(s1[i-1])][n2n('S')][n2n('S')];
                ss = ss > ss1? ss : ss1;
                ss2 = score[i][j-1][k] + scoring[n2n('S')][n2n(s2[j-1])][n2n('S')];
                ss = ss > ss2? ss : ss2;
                ss3 = score[i][j][k-1] + scoring[n2n('S')][n2n('S')][n2n(s3[k-1])];
                ss = ss > ss3? ss : ss3;
                // one indel
                ss1 = score[i-1][j-1][k] + scoring[n2n(s1[i-1])][n2n(s2[j-1])][n2n('S')];
                ss = ss > ss1? ss : ss1;
                ss2 = score[i-1][j][k-1] + scoring[n2n(s1[i-1])][n2n('S')][n2n(s3[k-1])];
                ss = ss > ss2? ss : ss2;
                ss3 = score[i][j-1][k-1] + scoring[n2n('S')][n2n(s2[j-1])][n2n(s3[k-1])];
                ss = ss > ss3? ss : ss3;

                score[i][j][k] = ss;
                //printf("score[%ld][%ld][%ld] = %d\n", i, j, k, ss);
            }
        }
    }

    if(DEBUG == 1)
    {
        for(id i = 0; i <= lm; i++)
        {
            for(id j = 0; j <= ln; j++)
            {
                for(id k = 0; k <= lk; k++)
                    printf("%6d\t", score[i][j][k]);
            //cout << score[i][j] << " ";
                cout << endl;
            }
         }
    }

    cout << "Final score : " << score[lm][ln][lk] << endl;
    return;
}



void createscoring()
{
    for(id i = 0; i < 5; i++)
    {
        for(id j = 0; j < 5; j++)
        {
            for(id k = 0; k < 5; k++)
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
                if(DEBUG == 1)
                    printf("score[%ld][%ld][%ld] = %d\t", i, j, k, s);
            }
    //        cout << endl;
        }
    }
}

