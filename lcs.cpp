#include "boost/multi_array.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;
const int DEBUG = 0;

// test score
const int MATCH = 8;
const int MISMATCH = -5;
const int INDEL = -3;


void alignment(vector<char>& s1, vector<char>& s2, vector<char>& s3, int &m, int &n, int &k);

int main(int argc, char **argv)
{
    if(argc != 4)
    {
        cout << "Usage: ./program seq1 seq2 seq3" << endl;
        return -1;
    }

    ifstream fseq1(argv[1]);
    ifstream fseq2(argv[2]);
    ifstream fseq3(argv[3]);
    if(fseq1.fail()||fseq2.fail()||fseq3.fail())
    {
        cout << "Open sequence files failed!" << endl;
    }

    vector<char> seq1;
    vector<char> seq2;
    vector<char> seq3;
    string line;

    while(getline(fseq1,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            seq1.push_back(line[i]);
        }
    }

    while(getline(fseq2,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            seq2.push_back(line[i]);
        }
    }

    while(getline(fseq3,line))
    {
        for(int i = 0; i < line.size(); i++)
        {
            seq3.push_back(line[i]);
        }
    }
    if(DEBUG == 0)
    {
        for(int i = 0; i < seq1.size(); i++)
        {
           cout << seq1[i];
        }
        cout << endl;
        for(int i = 0; i < seq2.size(); i++)
        {
           cout << seq2[i];
        }
        cout << endl;
        for(int i = 0; i < seq3.size(); i++)
        {
           cout << seq3[i];
        }
        cout << endl;
    }


    // choose the third sequence to cut at middle
    int m = seq1.size();
    int n = seq2.size();
    int k = seq3.size();

    alignment(seq1, seq2, seq3, m, n, k);

    fseq1.close();
    fseq2.close();
    fseq3.close();

    return 0;
}

void alignment(vector<char>& s1, vector<char>& s2, vector<char>& s3, int &m, int &n, int &k)
{
    // choose the third sequence to cut at the middle point
    int mid = k;
    typedef boost::multi_array<int, 2> array_type;
    typedef array_type::index index;
    array_type score(boost::extents[m+1][n+1]);

    // 0 - (mid-1) || mid - (k-1)
    // for the first plate, all preview scores are ZERO
    
    for(index i = 0; i <= m; i++)
    {
        score[i][0] = INDEL * i;
    }
    for(index i = 0; i <= n; i++)
    {
        score[0][i] = INDEL * i;
    }

    mid = 2;
    for(index ind = 1; ind < mid; ind++)
    {
        int ss = 0;
        int ss1 = 0;
        int ss2 = 0;
        // seq3 = ATCG or indel
        char kc = s3[ind];
        //cout << kc;
        // no the first plate, already has the preview score above in ind-1
        for(int i = 1; i <= m; i++)
        {
            for(int j = 1; j <= n; j++)
            {
                ss = 0;
                ss1 = 0;
                ss2 = 0;
                if(s1[i-1] == s2[j-1])
                {
                    ss = MATCH + score[i-1][j-1];
                    //cout << "Match " << s2[i-1] << " and" << s2[j-1] << endl;
                }
                else{
                    ss = MISMATCH + score[i-1][j-1];
                }
                // mismatch or indel
                {
                    ss1 = score[i-1][j] + INDEL;
                    ss2 = score[i][j-1] + INDEL;
                    ss1 = ss1 > ss2? ss1 : ss2;
                }
                ss = ss > ss1? ss : ss1;
                score[i][j] = ss;
            }
        }
    }

    if(DEBUG == 0)
    {
        for(index i = 0; i <= m; i++)
        {
            for(index j = 0; j <= n; j++)
            {
                printf("%6d\t", score[i][j]);
            //cout << score[i][j] << " ";
            }
            cout << endl;
         }
    }
    return;
}
