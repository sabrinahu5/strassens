#include <stdio.h>
#include <iostream>
#include <random>
#include <string.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

// standard matrix multiplication
vector<vector<int> > standard(vector<vector<int> >& M1, vector<vector<int> >& M2, int n) {
    vector<vector<int> > ans(n, vector<int>(n,0));
    
    int entry;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            entry = 0;
            for (int k = 0; k < n; k++) {
                entry += M1[i][k] * M2[k][j];
            }
            ans[i][j] = entry;

        }
    }

    return ans;
}

// implement normal strassen's?



// modified strassen's
vector<vector<int> > strassens(vector<vector<int> >& M1, vector<vector<int> >& M2, int n, int cp) {
    
    if (n <= cp) {
        return standard(M1, M2, n);
    } 

    if (n % 2 == 1) {
        // add padding of zeros to matrix
        M1.push_back(vector<int>(n,0));
        M2.push_back(vector<int>(n,0));
        for (int i = 0; i < M1.size(); i++) {
            M1[i].push_back(0);
            M2[i].push_back(0);
        }
        n++;
    }

    // strassen's!
    vector<vector<int> > P1, P2, P3, P4, P5, P6, P7;

    vector<vector<int> > sub1(n/2, vector<int>(n/2,0));
    vector<vector<int> > sub2(n/2, vector<int>(n/2,0));
    vector<vector<int> > ans(n, vector<int>(n,0));

    // P1 = A(F-H)
    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i][j];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i][j+n/2]+M2[i+n/2][j+n/2];
        }
    }

    P1 = strassens(sub1,sub2,n/2,cp);

    // P2 = (A+B)H

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i][j]+M1[i][j+n/2];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i+n/2][j+n/2];
        }
    }

    P2 = strassens(sub1,sub2,n/2,cp);

    // P3 = (C+D)E

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i+n/2][j]+M1[i+n/2][j+n/2];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i][j];
        }
    }

    P3 = strassens(sub1,sub2,n/2,cp);

    // P4 = D(G-E)

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i+n/2][j+n/2];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i+n/2][j]-M2[i][j];
        }
    }

    P4 = strassens(sub1,sub2,n/2,cp);

    // P5 = (A+D)(E+H)

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i][j]+M1[i+n/2][j+n/2];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i][j]+M2[i+n/2][j+n/2];
        }
    }

    P5 = strassens(sub1,sub2,n/2,cp);

    // P6 = (B-D)(G+H)

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i][j+n/2]-M1[i+n/2][j+n/2];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i+n/2][j]+M2[i+n/2][j+n/2];
        }
    }

    P6 = strassens(sub1,sub2,n/2,cp);

    // P7 = (C-A)(E+F)

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub1[i][j] = M1[i+n/2][j]-M1[i][j];
        }
    }

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            sub2[i][j] = M2[i][j]+M2[i][j+n/2];
        }
    }

    P7 = strassens(sub1,sub2,n/2,cp);

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            ans[i][j] = -P2[i][j]+P4[i][j]+P5[i][j]+P6[i][j];
        }
            
        for (int j = n/2; j < n; j++) {
            ans[i][j] = P1[i][j-n/2]+P2[i][j-n/2];
        }
    }
    for (int i = n/2; i < n; i++) {
        for (int j = 0; j < n/2; j++) {
            ans[i][j] = P3[i-n/2][j]+P4[i-n/2][j];
        }
            
        for (int j = n/2; j < n; j++) {
            ans[i][j] = P1[i-n/2][j-n/2]-P3[i-n/2][j-n/2]+P5[i-n/2][j-n/2]+P7[i-n/2][j-n/2];
        }
    }

    return ans;

}

void test_strassens(int n) {

    // generate random matrix
    vector<vector<int> > M1(n, vector<int>(n,0));
    vector<vector<int> > M2(n, vector<int>(n,0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int rand1 = rand() % static_cast<int>(3);
            int rand2 = rand() % static_cast<int>(3);

            M1[i][j] = rand1;
            M2[i][j] = rand2;
        }
    }

    auto start_time1 = std::chrono::high_resolution_clock::now();

    vector<vector<int> > ans1 = standard(M1,M2,n);

    auto end_time1= std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1);

    auto start_time2 = std::chrono::high_resolution_clock::now();

    vector<vector<int> > ans2 = strassens(M1,M2,n,2);

    auto end_time2= std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time2 - start_time2);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << M1[i][j] << " ";
        }
        cout << std::endl;
    }

    cout << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << M2[i][j] << " ";
        }
        cout << std::endl;
    }

    cout << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << ans1[i][j] << " ";
        }
        cout << std::endl;
    }

    std::cout << "Elapsed time: " << duration1.count() << " milliseconds\n";

    cout << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << ans2[i][j] << " ";
        }
        cout << std::endl;
    }

    std::cout << "Elapsed time: " << duration2.count() << " milliseconds\n";

}


int main(int argc, char *argv[]) {

    int n = atoi(argv[1]);

    // test stuff
    test_strassens(n);

}
