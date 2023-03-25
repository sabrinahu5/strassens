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
// change

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



// strassen's 



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

    auto start_time = std::chrono::high_resolution_clock::now();

    vector<vector<int> > ans = standard(M1,M2,n);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

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
            cout << ans[i][j] << " ";
        }
        cout << std::endl;
    }

    std::cout << "Elapsed time: " << duration.count() << " milliseconds\n";

}


int main() {

    // test stuff
    test_strassens(100);
    

}
