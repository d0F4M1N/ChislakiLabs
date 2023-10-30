#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// ������� ��� ������ �������
void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }
}

// ������� ��� ������ �������
void printVector(const vector<double>& vector) {
    for (double value : vector) {
        cout << value << " ";
    }
    cout << endl;
}

// ��������������� ������� ��� ������ ������ � ������� �������� �������� �� �������
void swapRows(vector<vector<double>>& matrix, vector<double>& b, int row1, int row2) {
    swap(matrix[row1], matrix[row2]);
    swap(b[row1], b[row2]);
}

void gaussEliminationRecursive(vector<vector<double>>& matrix, vector<double>& b, int i) {
    int n = matrix.size();

    if (i == n) {
        return; // ������� ������: ��������� ����� �������
    }

    // ����� ������������� �������� � �������
    int maxRow = i;
    for (int k = i + 1; k < n; ++k) {
        if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
            maxRow = k;
        }
    }

    // ����� �����
    swapRows(matrix, b, i, maxRow);

    // ������ ��� ������
    for (int k = i + 1; k < n; ++k) {
        double factor = matrix[k][i] / matrix[i][i];
        for (int j = i; j < n; ++j) {
            matrix[k][j] -= factor * matrix[i][j];
        }
        b[k] -= factor * b[i];
    }

    // ����������� ����� ��� ��������� ������
    gaussEliminationRecursive(matrix, b, i + 1);
}

// �������� ������� ��� ������ ������ � ������� �������� �������� �� �������
void gaussElimination(vector<vector<double>>& matrix, vector<double>& b, vector<double>& residuals, double& norm) {
    int n = matrix.size();

    // ���������� �������� ������� � ������� b
    vector<vector<double>> originalMatrix = matrix;
    vector<double> originalB = b;

    gaussEliminationRecursive(matrix, b, 0);

    // �������� ��� ������
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / matrix[i][i];
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= matrix[k][i] * x[i];
        }
    }

    // ���������� ������� ������� F = Ax - b
    residuals.resize(n);
    for (int i = 0; i < n; ++i) {
        residuals[i] = 0;
        for (int j = 0; j < n; ++j) {
            residuals[i] += originalMatrix[i][j] * x[j];
        }
        residuals[i] -= originalB[i];
    }

    // ���������� ����� ������� ������� Delta = max(|F_i|)
    norm = 0;
    for (int i = 0; i < n; ++i) {
        norm = max(norm, abs(residuals[i]));
    }

    // ����� ����������
    cout << "Solution (x) using Gauss Elimination: ";
    printVector(x);
    cout << "Residuals (F): ";
    printVector(residuals);
    cout << "Norm of Residuals (Delta): " << norm << endl;
}

// ������� ��� LDL^T-������������ ������� A
void LDLTFactorization(const vector<vector<double>>& A, vector<vector<double>>& L, vector<double>& D) {
    int n = A.size();

    // ������������� L � D
    L.resize(n, vector<double>(n, 0.0));
    D.resize(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = A[i][j];
            for (int k = 0; k < j; ++k) {
                sum -= L[i][k] * D[k] * L[j][k];
            }

            if (i == j) {
                // ������������ �������� D
                D[i] = sum;
            }
            else {
                // �������� ���������������� ������� L
                L[i][j] = sum / D[j];
            }
        }
    }
}

// ������� ��� ������� ������� ������� LDL^T-������������
void solveLDLTFactorization(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {
    int n = A.size();

    // LDL^T-������������
    vector<vector<double>> L;
    vector<double> D;
    LDLTFactorization(A, L, D);

    // ������� ������� Ly = b
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = b[i];
        for (int j = 0; j < i; ++j) {
            sum -= L[i][j] * y[j];
        }
        y[i] = sum / L[i][i];
    }

    // ������� ������� L^Tz = y
    vector<double> z(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= L[j][i] * z[j];
        }
        z[i] = sum;
    }

    // ������� ������� Dx = z
    x.resize(n, 0.0);
    for (int i = 0; i < n; ++i) {
        x[i] = z[i] / D[i];
    }
}

int main() {
    // ���� ������� ������������� A � ������ b
    vector<vector<double>> A = { {2.30, 3.50, 1.70},
                                {5.70, -2.70, 2.30},
                                {-0.80, 5.30, -1.80} };
    vector<double> b = { -6.49, 19.20, -5.09 };

    // ����� �������� ������
    cout << "Matrix A:" << endl;
    printMatrix(A);
    cout << "Vector b:" << endl;
    printVector(b);

    // ����� ������� ��� ������� ������� ������� ������ � ����������� ������� ������� � ������������� �����������
    vector<double> residuals;
    double norm;
    gaussElimination(A, b, residuals, norm);

    // ������� ������� ������� LDL^T-������������
    vector<double> x;
    solveLDLTFactorization(A, b, x);

    // ����� ����������
    cout << "Solution (x) using LDL^T Factorization: ";
    printVector(x);

    return 0;
}
