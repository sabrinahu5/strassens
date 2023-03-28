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
            sub2[i][j] = M2[i][j+n/2]-M2[i+n/2][j+n/2];
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

void count_triangles() {

    srand(time(0));

    // p = 0.01, 0.02, 0.03, 0.04, 0.05

    // generate graph with 1024 vertices
    vector<vector<int> > A(1024, vector<int>(1024,0));

    double p = 0.01;

    // loop through all pairs of vertices
    for (int i = 0; i < 1024; i++) {
        for (int j = i+1; j < 1024; j++) {
            double rand_num = (double) rand() / RAND_MAX; // generate a random number between 0 and 1

            if (rand_num <= p) { // if the random number is less than or equal to the probability
                A[i][j] = 1; // set the adjacency matrix entry to 1
                A[j][i] = 1; // since the graph is undirected, set the symmetric entry to 1 as well
            }
        }
    }

    int cp = 16;

    vector<vector<int> > A_squared = strassens(A, A, 1024, cp);
    vector<vector<int> > A_cubed = strassens(A, A_squared, 1024, cp);

    double num_triangles = 0;
    for (int i = 0; i < 1024; i++) {
        num_triangles += (double) A_cubed[i][i];
    }
    num_triangles = num_triangles / 6;

    cout << "Num triangles, p = " << p << ": " << num_triangles << std::endl;

}


int main(int argc, char *argv[]) {

    //count_triangles();

    if (argc != 4) {
        std::cout << "Usage: ./strassen 0 dimension inputfile\n";
        return 0;
    }

    int n = atoi(argv[2]);
    int cp = 16;

    // take in input matrix
    vector<vector<int> > M1(n, vector<int>(n,0));
    vector<vector<int> > M2(n, vector<int>(n,0));


    if (atoi(argv[1]) != 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int rand1 = rand() % static_cast<int>(3);
                int rand2 = rand() % static_cast<int>(3);

                M1[i][j] = rand1;
                M2[i][j] = rand2;
            }
        }
    } else if (atoi(argv[1]) == 0) {
        ifstream testfile(argv[3]);
        string line;

        int i = 0;
        int a, b;
        while (getline(testfile, line)) {
        
            if (i < n*n) {
                a = i / n;
                b = i % n;
                M1[a][b] = stoi(line);
                i++;
            } else {
                a = (i / n) - n;
                b = i % n;
                M2[a][b] = stoi(line);
                i++;
            }
        }

        testfile.close();
    }

    // test standard multiplication
    auto start_time1 = std::chrono::high_resolution_clock::now();

    vector<vector<int> > ans1 = standard(M1,M2,n);

    /*for (int i = 0; i < n; i++) {
        cout << ans1[i][i] << std::endl;
    }*/

    auto end_time1= std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time1);

    // std::cout << "Standard elapsed time: " << duration1.count() << " microseconds\n";

    // test modified strassen's
    auto start_time2 = std::chrono::high_resolution_clock::now();
    
    vector<vector<int> > ans2 = strassens(M1,M2,n,cp);

    for (int i = 0; i < n; i++) {
        cout << ans2[i][i] << std::endl;
    }

    auto end_time2= std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time2);
    // std::cout << "Strassen's elapsed time: " << duration2.count() << " microseconds\n"; 
    

    
}