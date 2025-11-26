#include <iostream>
#include <vector>
using namespace std;


vector<vector<int>> addMatrix(const vector<vector<int>>& A,
                              const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];

    return C;
}

vector<vector<int>> subMatrix(const vector<vector<int>>& A,
                              const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];

    return C;
}


vector<vector<int>> standardMultiply(const vector<vector<int>>& A,
                                     const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];

    return C;
}


vector<vector<int>> strassen(const vector<vector<int>>& A,
                             const vector<vector<int>>& B) {
    int n = A.size();


    if (n == 1) {
        return {{A[0][0] * B[0][0]}};
    }

    int k = n / 2;


    vector<vector<int>>
        A11(k, vector<int>(k)), A12(k, vector<int>(k)),
        A21(k, vector<int>(k)), A22(k, vector<int>(k)),
        B11(k, vector<int>(k)), B12(k, vector<int>(k)),
        B21(k, vector<int>(k)), B22(k, vector<int>(k));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j+k];
            A21[i][j] = A[i+k][j];
            A22[i][j] = A[i+k][j+k];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j+k];
            B21[i][j] = B[i+k][j];
            B22[i][j] = B[i+k][j+k];
        }
    }

    auto M1 = strassen(addMatrix(A11, A22), addMatrix(B11, B22));
    auto M2 = strassen(addMatrix(A21, A22), B11);
    auto M3 = strassen(A11, subMatrix(B12, B22));
    auto M4 = strassen(A22, subMatrix(B21, B11));
    auto M5 = strassen(addMatrix(A11, A12), B22);
    auto M6 = strassen(subMatrix(A21, A11), addMatrix(B11, B12));
    auto M7 = strassen(subMatrix(A12, A22), addMatrix(B21, B22));

    // Compute C quadrants
    auto C11 = addMatrix(subMatrix(addMatrix(M1, M4), M5), M7);
    auto C12 = addMatrix(M3, M5);
    auto C21 = addMatrix(M2, M4);
    auto C22 = addMatrix(subMatrix(addMatrix(M1, M3), M2), M6);


    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }
    }

    return C;
}



void printMatrix(const vector<vector<int>>& M) {
    for (auto& row : M) {
        for (int x : row) cout << x << " ";
        cout << endl;
    }
}



int main() {
    vector<vector<int>> A = {
        {1, 2},
        {3, 4}
    };

    vector<vector<int>> B = {
        {5, 6},
        {7, 8}
    };

    cout << "Standard Matrix Multiplication:\n";
    auto C1 = standardMultiply(A, B);
    printMatrix(C1);

    cout << "\nStrassen Matrix Multiplication:\n";
    auto C2 = strassen(A, B);
    printMatrix(C2);

    return 0;
}