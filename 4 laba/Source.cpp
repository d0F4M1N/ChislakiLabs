#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

void gaussElimination(vector<vector<double>> A, vector<double> B, vector<double>& x, int n) {
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[maxRow], A[i]);
        swap(B[maxRow], B[i]);

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x.push_back(B[i]);
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    for (int i = 0; i < n; i++) {
        cout << "a[" << i << "] = " << x[i] << endl;
    }
}

void drawGraph(vector<double> X, vector<double> a, int m, int n) {
    ofstream dataFile("data.txt");
    if (!dataFile.is_open()) {
        cerr << "Unable to open data file" << endl;
        return;
    }
    for (int i = 0; i < n; i++)
    {
        double Y = 0;
        for (int q = 0; q <= m; q++)
        {
            Y += a[q] * pow(X[i], q);
        }
        dataFile << X[i] << " " << Y << endl;
    }
    dataFile.close();
    ofstream scriptFile("script.plt");
    if (!scriptFile.is_open()) {
        cerr << "Unable to open script file" << endl;
        return;
    }
    scriptFile << "set term png\n";
    scriptFile << "set output 'graph.png'\n";
    scriptFile << "set multiplot\n";
    scriptFile << "plot 'data.txt' with lines\n";
    scriptFile << "plot 'oldData.txt' with lines\n";
    scriptFile << "unset multiplot";
    scriptFile.close();
    system("gnuplot script.plt");
}

int main()
{
    vector<double> X = { 0,1,2,3,4,5,6,7,8,9,10 };
    vector<double> Y = { 3, 87, 156, 210, 238, 252, 239, 211, 158, 90, -5 };
    double n = 11;
    int m = 2;
    vector<double> powerX;
    for (int i = 1; i <= 2 * m; i++)
    {
        double sum = 0;
        for (int q = 0; q < n; q++)
        {
            sum += pow(X[q], i);
        }
        powerX.push_back(sum);
    }
    vector<vector<double>> sumX = { {n,0,0},{0,0,0},{0,0,0} };
    for (int i = 0; i <= m; i++)
    {
        for (int q = 0; q <= m; q++)
        {
            if (i + q >= 1)
            {
                sumX[i][q] = powerX[i + q - 1];
            }
        }
    }
    for (int i = 0; i < m + 1; i++)
    {
        for (int q = 0; q < m + 1; q++)
        {
            cout << sumX[i][q] << " ";
        }
        cout << endl;
    }
    vector<double> praw;
    for (int i = 0; i < m + 1; i++)
    {
        double sum = 0;
        for (int q = 0; q < n; q++)
        {
            sum += Y[q] * pow(X[q], i);
        }
        praw.push_back(sum);
    }
    vector<double> a;
    gaussElimination(sumX, praw, a, m + 1);
    double deviation = 0;
    for (int i = 0; i < n; i++)
    {
        double sum = Y[i];
        for (int q = 0; q <= m; q++)
        {
            sum -= a[q] * pow(X[i], q);
        }
        deviation += sum * sum / (n - m - 1);
    }
    cout << "Deviation = " << sqrt(deviation) << endl;
    drawGraph(X, a, m, n);
}
