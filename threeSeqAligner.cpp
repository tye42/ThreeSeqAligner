#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <string>

using namespace std;

// SP score for three sequences alignment, 0: A, 1:C, 2:G, 3:T, 4:-
int SP[5][5][5];
// convert A,C,G,T to integer
int c2i(char c);
// convert integer to A,C,G,T
char i2c[] = {'A', 'C', 'G', 'T', '-'};
struct alignment{
    string v; string w; string u;
    int score;
    alignment(){ v = ""; w = ""; u = ""; }
    alignment(string v, string w, string u, int s): v(v), w(w), u(u), score(s){}
    alignment& operator=(const alignment a){
        v = a.v; w = a.w; u = a.u;
        score = a.score;
        return *this;
    }
    // concatenate two alignment
    alignment operator+(const alignment& a) const{
        return alignment(v + a.v, w + a.w, u + a.u, score + a.score);
    }
    // prepend a char to each row of the alignment
    alignment& prepend(char i, char j, char k){
        v = i + v; w = j + w; u = k + u;
        return *this;
    }
};
// print out the alignment, score, length and identity
void printAlignment(alignment& a, int width, bool printAlign=false);
// pair-wise alignment for initialize the matrix, with traceback matrix
void twoSeqAlignTrace(int v[], int w[], int n, int m, int s[], int prev[]);
// three sequences alignment in cubic space, with tracback matrix and return an alignment
alignment threeSeqAlignTrace(int v[], int w[], int u[], int n, int m, int l);
// three sequence alignment in quadratic space via divide-and-conquer
alignment threeSeqAlignDC(int v[], int w[], int u[], int n, int m, int l);

int main(int argc, char *argv[]) {
    // create SP score
    int scoreM[5][5];
    ifstream fin("scoreMatrix.txt");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            fin >> scoreM[i][j];
        }
    }
    fin.close();
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 5; ++k) {
                SP[i][j][k] = scoreM[i][j] + scoreM[i][k] + scoreM[j][k];
            }
        }
    }
    // read in three DNA sequences
    string sv, sw, su;
    fin.open(argv[1]);
    string str;
    while(getline(fin, str)) sv += str;
    fin.close();
    fin.open(argv[2]);
    while(getline(fin, str)) sw += str;
    fin.close();
    fin.open(argv[3]);
    while(getline(fin, str)) su += str;
    fin.close();

    // convert string to int array
    int *v = new int[sv.size()];
    int *w = new int[sw.size()];
    int *u = new int[su.size()];
    int n = sv.size(), m = sw.size(), l = su.size();

    for (int i = 0; i < sv.size(); ++i) v[i] = c2i(sv[i]);
    for (int i = 0; i < sw.size(); ++i) w[i] = c2i(sw[i]);
    for (int i = 0; i < su.size(); ++i) u[i] = c2i(su[i]);

    clock_t s = clock();
    alignment align = threeSeqAlignDC(v, w, u, n, m, l);
    clock_t e = clock();
    float diff = ((float)e - (float)s) / CLOCKS_PER_SEC;
    printf("# Sequences: %s: %d  %s: %d  %s: %d\n",
           argv[1], int(sv.size()), argv[2], int(sw.size()), argv[3], int(su.size()));
    cout << "# Running time: " << diff << endl;
    printAlignment(align, 100, true);

    delete[] v; delete[] w; delete[] u;
    return 0;
}

// convert A,C,G,T to integer
int c2i(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}
// pair-wise alignment for initialize the matrix, with traceback matrix
void twoSeqAlignTrace(int v[], int w[], int n, int m, int s[], int prev[]) {
    int wd = m + 1;
    // init first row
    s[0] = 0; // i * wd + j
    prev[0] = -1; // -1 as stop
    for (int j = 1; j < m+1; ++j) {
        // the third sequence is empty, so always a space
        s[j] = s[j-1] + SP[4][w[j-1]][4];
        prev[j] = 1;
    }
    // init first column
    for (int i = 1; i < n+1; ++i) {
        s[i * wd] = s[(i-1) * wd] + SP[v[i-1]][4][4];
        prev[i * wd] = 0;
    }
    // fill the rest matrix
    for (int i = 1; i < n+1; ++i) {
        for (int j = 1; j < m+1; ++j) {
            int neighbor[] = {
                    s[(i-1) * wd + j] + SP[v[i-1]][4][4],
                    s[i * wd + j-1] + SP[4][w[j-1]][4],
                    s[(i-1) * wd + j-1] + SP[v[i-1]][w[j-1]][4]
            };
            s[i * wd + j] = *max_element(neighbor, neighbor+3);
            if (s[i * wd + j] == neighbor[0]){
                prev[i * wd + j] = 0; // v,-,-
            }else if (s[i * wd + j] == neighbor[1]){
                prev[i * wd + j] = 1; // -,w,-
            }else{
                prev[i * wd + j] = 3; // v,w,-
            }

        }
    }
}
// three sequences alignment in cubic space, with tracback matrix and return an alignment
alignment threeSeqAlignTrace(int v[], int w[], int u[], int n, int m, int l){
    int wd = m + 1;

    int** score = new int*[l+1];
    int** prev = new int*[l+1]; // traceback matrix
    for (int i = 0; i < l + 1; ++i) {
        score[i] = new int[(n+1)*(m+1)];
        prev[i] = new int[(n+1)*(m+1)];
    }
    twoSeqAlignTrace(v, w, n, m, score[0], prev[0]); // initialize the first layer

    for (int k = 1; k < l+1; ++k) {
        // init the first row of current layer with w and u
        score[k][0] = score[k-1][0] + SP[4][4][u[k-1]];
        prev[k][0] = 2;
        for (int j = 1; j < m+1; ++j) {
            int neighbor[] = {
                    score[k-1][j] + SP[4][4][u[k-1]],
                    score[k][j-1] + SP[4][w[j-1]][4],
                    score[k-1][j-1] + SP[4][w[j-1]][u[k-1]]
            };
            score[k][j] = *max_element(neighbor, neighbor+3);
            if (score[k][j] == neighbor[0]){
                prev[k][j] = 2;
            }else if (score[k][j] == neighbor[1]){
                prev[k][j] = 1;
            }else{
                prev[k][j] = 5;
            }
        }
        // init the first col of current layer with v and u
        for (int i = 1; i < n+1; ++i) {
            int neighbor[] = {
                    score[k-1][i*wd] + SP[4][4][u[k-1]],
                    score[k][(i-1)*wd] + SP[v[i-1]][4][4],
                    score[k-1][(i-1)*wd] + SP[v[i-1]][4][u[k-1]]
            };
            score[k][i*wd] = *max_element(neighbor, neighbor+3);
            if (score[k][i*wd] == neighbor[0]){
                prev[k][i*wd] = 2;
            }else if (score[k][i*wd] == neighbor[1]){
                prev[k][i*wd] = 0;
            }else{
                prev[k][i*wd] = 4;
            }
        }
        // fill the rest of current layer
        for (int i = 1; i < n+1; ++i) {
            for (int j = 1; j < m+1; ++j) {
                int neighbor[] = {
                        score[k][(i-1)*wd+j] + SP[v[i-1]][4][4],
                        score[k][i*wd+j-1] + SP[4][w[j-1]][4],
                        score[k-1][i*wd+j] + SP[4][4][u[k-1]],
                        score[k][(i-1)*wd+j-1] + SP[v[i-1]][w[j-1]][4],
                        score[k-1][(i-1)*wd+j] + SP[v[i-1]][4][u[k-1]],
                        score[k-1][i*wd+j-1] + SP[4][w[j-1]][u[k-1]],
                        score[k-1][(i-1)*wd+j-1] + SP[v[i-1]][w[j-1]][u[k-1]]
                };
                int* tmp = max_element(neighbor, neighbor+7);
                score[k][i*wd+j] = *tmp;
                prev[k][i*wd+j] = int(distance(neighbor, tmp));
            }
        }
    }
    // traceback
    alignment align;
    align.score = score[l][n*wd+m];
    int i = n, j = m, k = l;
    while (i > 0 || j > 0 || k > 0){
        if (prev[k][i*wd+j] == 0){
            align.prepend(i2c[v[i-1]], '-', '-'); i--;
        } else if (prev[k][i*wd+j] == 1){
            align.prepend('-', i2c[w[j-1]], '-'); j--;
        } else if (prev[k][i*wd+j] == 2){
            align.prepend('-', '-', i2c[u[k-1]]); k--;
        } else if (prev[k][i*wd+j] == 3){
            align.prepend(i2c[v[i-1]], i2c[w[j-1]], '-'); i--; j--;
        } else if (prev[k][i*wd+j] == 4){
            align.prepend(i2c[v[i-1]], '-', i2c[u[k-1]]); i--; k--;
        } else if (prev[k][i*wd+j] == 5){
            align.prepend('-', i2c[w[j-1]], i2c[u[k-1]]); j--; k--;
        } else{
            align.prepend(i2c[v[i-1]], i2c[w[j-1]], i2c[u[k-1]]); i--; j--; k--;
        }
    }
    for (int d = 0; d < l + 1; ++d) {
        delete[] score[d]; delete[] prev[d];
    }
    delete[] score; delete[] prev;
    return align;
}
// three sequence alignment in quadratic space via divide-and-conquer
alignment threeSeqAlignDC(int v[], int w[], int u[], int n, int m, int l){
    int midu = l / 2;
    int wd = m + 1;
    // --------------------forward part start (v, w, u_1_midu) -----------------------
    int** forward = new int*[2]; // two layers, one for previous, one for current
    forward[0] = new int[(n+1)*(m+1)]; // i * wd + j
    forward[1] = new int[(n+1)*(m+1)];
    // -------------init the first layer start-----------------------
    // init first row
    forward[0][0] = 0;
    for (int j = 1; j < m+1; ++j) {
        forward[0][j] = forward[0][j-1] + SP[4][w[j-1]][4];
    }
    // init first col
    for (int i = 1; i < n+1; ++i) {
        forward[0][i * wd] = forward[0][(i-1) * wd] + SP[v[i-1]][4][4];
    }
    // fill the rest matrix
    for (int i = 1; i < n+1; ++i) {
        for (int j = 1; j < m+1; ++j) {
            forward[0][i * wd + j] = max(
                    max(forward[0][(i-1) * wd + j] + SP[v[i-1]][4][4],
                        forward[0][i * wd + j-1] + SP[4][w[j-1]][4]
                    ),
                    forward[0][(i-1) * wd + j-1] + SP[v[i-1]][w[j-1]][4]
            );
        }
    }
    // -------------init the first layer end -----------------------
    // -------------compute the rest layers-------------------------
    for (int k = 1; k < midu+1; ++k) {
        // index of current layer, if sInd = 1, the previous ind is 1-sInd = 0; else the previous ind is 1-sInd = 1
        int sInd = k % 2;
        // init the first row of current layer with w and u
        forward[sInd][0] = forward[1-sInd][0] + SP[4][4][u[k-1]];
        for (int j = 1; j < m+1; ++j) {
            forward[sInd][j] = max(
                    max(forward[1-sInd][j] + SP[4][4][u[k-1]],
                        forward[sInd][j-1] + SP[4][w[j-1]][4]
                    ),
                    forward[1-sInd][j-1] + SP[4][w[j-1]][u[k-1]]
            );
        }
        // init the first col of current layer with v and u
        for (int i = 1; i < n+1; ++i) {
            forward[sInd][i*wd] = max(
                    max(forward[1-sInd][i*wd] + SP[4][4][u[k-1]],
                        forward[sInd][(i-1)*wd] + SP[v[i-1]][4][4]
                    ),
                    forward[1-sInd][(i-1)*wd] + SP[v[i-1]][4][u[k-1]]
            );
        }
        // fill the rest of current layer
        for (int i = 1; i < n+1; ++i) {
            for (int j = 1; j < m+1; ++j) {
                int max = forward[sInd][(i-1)*wd+j] + SP[v[i-1]][4][4];
                int tmp = forward[sInd][i*wd+j-1] + SP[4][w[j-1]][4];
                max = max > tmp? max:tmp;
                tmp =  forward[1-sInd][i*wd+j] + SP[4][4][u[k-1]];
                max = max > tmp? max:tmp;
                tmp = forward[sInd][(i-1)*wd+j-1] + SP[v[i-1]][w[j-1]][4];
                max = max > tmp? max:tmp;
                tmp =  forward[1-sInd][(i-1)*wd+j] + SP[v[i-1]][4][u[k-1]];
                max = max > tmp? max:tmp;
                tmp =  forward[1-sInd][i*wd+j-1] + SP[4][w[j-1]][u[k-1]];
                max = max > tmp? max:tmp;
                tmp =  forward[1-sInd][(i-1)*wd+j-1] + SP[v[i-1]][w[j-1]][u[k-1]];
                max = max > tmp? max:tmp;
                forward[sInd][i*wd+j] = max;
            }
        }
    }
    // --------------------forward part end (v, w, u/2) -----------------------------------
    // --------------------reverse part start (rv, rw, ru_midu+1_l) -----------------------
    int** reverse = new int*[2];
    reverse[0] = new int[(n+1)*(m+1)]; // i * wd + j
    reverse[1] = new int[(n+1)*(m+1)];
    // -------------init the first layer start-----------------------
    // init first row
    reverse[0][0] = 0; // i * wd + j
    for (int j = 1; j < m+1; ++j) {
        reverse[0][j] = reverse[0][j-1] + SP[4][w[m-j]][4];
    }
    // init first col
    for (int i = 1; i < n+1; ++i) {
        reverse[0][i * wd] = reverse[0][(i-1) * wd] + SP[v[n-i]][4][4];
    }
    // fill the rest matrix
    for (int i = 1; i < n+1; ++i) {
        for (int j = 1; j < m+1; ++j) {
            reverse[0][i * wd + j] = max(
                    max(reverse[0][(i-1) * wd + j] + SP[v[n-i]][4][4],
                        reverse[0][i * wd + j-1] + SP[4][w[m-j]][4]
                    ),
                    reverse[0][(i-1) * wd + j-1] + SP[v[n-i]][w[m-j]][4]
            );
        }
    }
    // -------------init the first layer end -----------------------
    // -------------compute the rest layers-------------------------
    for (int k = 1; k < l-midu+1; ++k) {
        int sInd = k % 2;
        reverse[sInd][0] = reverse[1-sInd][0] + SP[4][4][u[l-k]];
        for (int j = 1; j < m+1; ++j) {
            reverse[sInd][j] = max(
                    max(reverse[1-sInd][j] + SP[4][4][u[l-k]],
                        reverse[sInd][j-1] + SP[4][w[m-j]][4]
                    ),
                    reverse[1-sInd][j-1] + SP[4][w[m-j]][u[l-k]]
            );
        }
        for (int i = 1; i < n+1; ++i) {
            reverse[sInd][i*wd] = max(
                    max(reverse[1-sInd][i*wd] + SP[4][4][u[l-k]],
                        reverse[sInd][(i-1)*wd] + SP[v[n-i]][4][4]
                    ),
                    reverse[1-sInd][(i-1)*wd] + SP[v[n-i]][4][u[l-k]]
            );
        }
        // fill the rest of current layer
        for (int i = 1; i < n+1; ++i) {
            for (int j = 1; j < m+1; ++j) {
                int max = reverse[sInd][(i-1)*wd+j] + SP[v[n-i]][4][4];
                int tmp = reverse[sInd][i*wd+j-1] + SP[4][w[m-j]][4];
                max = max > tmp? max:tmp;
                tmp = reverse[1-sInd][i*wd+j] + SP[4][4][u[l-k]];
                max = max > tmp? max:tmp;
                tmp = reverse[sInd][(i-1)*wd+j-1] + SP[v[n-i]][w[m-j]][4];
                max = max > tmp? max:tmp;
                tmp = reverse[1-sInd][(i-1)*wd+j] + SP[v[n-i]][4][u[l-k]];
                max = max > tmp? max:tmp;
                tmp = reverse[1-sInd][i*wd+j-1] + SP[4][w[m-j]][u[l-k]];
                max = max > tmp? max:tmp;
                tmp = reverse[1-sInd][(i-1)*wd+j-1] + SP[v[n-i]][w[m-j]][u[l-k]];
                max = max > tmp? max:tmp;
                reverse[sInd][i*wd+j] = max;
            }
        }
    }
    // --------------------reverse part end (rv, rw, ru_midu+1_l) -------------------------
    int midv = 0, midw = 0;
    int maxSum = numeric_limits<int>::min();
    // we have two layers, find and rind compute where is the last layer locate
    int find = midu % 2, rind = (l-midu) % 2;
    for (int i = 0; i < n+1; ++i) {
        for (int j = 0; j < m+1; ++j) {
            int tmp = forward[find][i*wd+j] + reverse[rind][(n-i)*wd+m-j];
            if (tmp > maxSum){
                midv = i;
                midw = j;
                maxSum = tmp;
            }
        }
    }
    for (int d = 0; d < 2; ++d) {
        delete[] forward[d]; delete[] reverse[d];
    }
    delete[] forward; delete[] reverse;

    alignment forwardAlign, reverseAlign;
    if (midv > 1 && midw > 1 && midu > 1){
        forwardAlign = threeSeqAlignDC(v, w, u, midv, midw, midu);
    } else{
        forwardAlign = threeSeqAlignTrace(v, w, u, midv, midw, midu); }
    if (n - midv > 1 && m - midw > 1 && l - midu > 1){
        reverseAlign = threeSeqAlignDC(v+midv, w+midw, u+midu, n-midv, m-midw, l-midu);
    } else{
        reverseAlign = threeSeqAlignTrace(v+midv, w+midw, u+midu, n-midv, m-midw, l-midu); }
    return forwardAlign + reverseAlign;
}
// print out the alignment, score, length and identity
void printAlignment(alignment& a, int width, bool printAlign){
    int addWd = 15;
    int len = a.v.size();

    cout << "#" << string(width+addWd, '=') << endl << "#" << endl;
    cout << "# Length: " << len << endl;
    cout << "# Score: " << a.score << endl;

    // count the number of perfected matched columns
    int identity = 0;
    for (int i = 0; i < len; ++i) {
        if (a.v[i] == a.w[i] && a.w[i] == a.u[i]){
            identity++;
        }
    }
    printf("# Identity: %d/%d (%.2f%%)\n", identity, len, identity*100.0/len);
    cout << "#" << endl;
    cout << "#" << string(width+addWd, '=') << endl;

    // print the alignment
    int startv = 1, startw = 1, startu = 1;
    int endv, endw, endu;

    if (printAlign) {
        cout << "#" << endl;
        int lines;
        if (len % width == 0) {
            lines = len / width;
        } else {
            lines = len / width + 1;
        }
        for (int i = 0; i < lines; ++i) {
            string tmpv = a.v.substr(i * width, width);
            string tmpw = a.w.substr(i * width, width);
            string tmpu = a.u.substr(i * width, width);

            int lenv = tmpv.size() - count(tmpv.begin(), tmpv.end(), '-');
            endv = startv + lenv - 1;
            printf("# %-5d %s %5d\n", startv, tmpv.c_str(), endv);
            startv += lenv;

            int lenw = tmpw.size() - count(tmpw.begin(), tmpw.end(), '-');
            endw = startw + lenw - 1;
            printf("# %-5d %s %5d\n", startw, tmpw.c_str(), endw);
            startw += lenw;

            int lenu = tmpu.size() - count(tmpu.begin(), tmpu.end(), '-');
            endu = startu + lenu - 1;
            printf("# %-5d %s %5d\n", startu, tmpu.c_str(), endu);
            startu += lenu;

            // if column perfected matched, print *
            cout << "#       ";
            for (int j = 0; j < tmpv.size(); ++j) {
                if (tmpv[j] == tmpw[j] && tmpw[j] == tmpu[j]) {
                    cout << "*";
                } else {
                    cout << " ";
                }
            }
            cout << endl << "#" << endl;
        }
        cout << "#" << string(width + addWd, '=') << endl << endl;
    }
}