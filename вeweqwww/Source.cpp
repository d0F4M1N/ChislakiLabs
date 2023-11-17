#include<iostream>
#include<iomanip>
#include<math.h>
using namespace std;

double* gauss(double** matrix, double* y, const int n, int accur);
void outputMatrix(double** matrix, double* y, const int n);
void initialisationConst(double** matrix, double* y);
void residVector(double** matrix, double* resVect, double* bVect, double* xVect, int n);
void residVectorMax(double* resVect, int n);

int main() {
    setlocale(LC_ALL, "Russian");
    const int n = 3;
    int accur = 100000;
    double* bVect, * xVect, * xVect2, * residVect;
    double** matrix = new double* [n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
    }
    bVect = new double[n];
    xVect = new double[n];
    xVect2 = new double[n];
    residVect = new double[n];
    for (int i = 0; i < n; i++) {
        residVect[i] = 0;
    }
    initialisationConst(matrix, bVect);
    outputMatrix(matrix, bVect, n);
    xVect = gauss(matrix, bVect, n, accur);
    cout << "answer:";
    for (int i = 0; i < n; i++) {
        cout << "x" << (i + 1) << " = " << xVect[i] << "; ";
    }
    // residual vector (вектор невязки)
    initialisationConst(matrix, bVect);
    residVector(matrix, residVect, bVect, xVect, n);
    //---------------------------------------------------------
    xVect2 = gauss(matrix, residVect, n, accur);
    cout << "answer2:";
    for (int i = 0; i < n; i++) {
        cout << "x" << (i + 1) << " = " << xVect2[i] << "; ";
    }
    cout << endl;
    double modul1, xMod1 = 0, modul2, xMod2 = 0, pogreshnost;
    for (int i = 0; i <= 2; i++) {
        xMod1 += pow((xVect2[i] - xVect[i]), 2);
    }
    for (int i = 0; i <= 2; i++) {
        xMod2 += pow(xVect[i], 2);
    }
    modul1 = pow((xMod1), 0.5);
    modul2 = pow((xMod2), 0.5);
    pogreshnost = modul1 / modul2;
    cout << "Pogreshnost= ";
    cout << pogreshnost << endl;
    // удаление массивов
    for (int i = 0; i < n; i++) {
        delete[]matrix[i];
    }
    delete[]matrix;
    delete[]bVect;
    delete[]xVect;
    delete[]xVect2;
    delete[]residVect;
}
void outputMatrix(double** matrix, double* y, const int n) {
    cout << endl << "matrix: " << endl;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {

            cout << matrix[j][i] << setw(10);

        }
        cout << "|" << y[j] << " " << endl;
    }
}
double* gauss(double** matrix, double* y, const int n, int accur)
{
    double* xVect, max, roundH;
    int k, index, round;
    const double eps = 0.00001;  // точность
    xVect = new double[n];
    k = 0;
    while (k < n)
    {
        // Поиск максимума
        max = abs(matrix[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(matrix[i][k]) > max)
            {
                max = abs(matrix[i][k]);
                index = i;
            }
        }

        if (max < eps)
        {
            // диаг элем = 0?
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        for (int j = 0; j < n; j++)
        {
            swap(matrix[k][j], matrix[index][j]);
        }
        swap(y[k], y[index]);
        if (k > 0)
            outputMatrix(matrix, y, n);
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = matrix[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++) {
                matrix[i][j] = matrix[i][j] / temp;
                round = matrix[i][j] * accur;
                roundH = round;
                matrix[i][j] = roundH / accur;
            }
            y[i] = y[i] / temp;
            round = y[i] * accur;
            roundH = round;
            y[i] = roundH / accur;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            y[i] = y[i] - y[k];
        }
        if (k == n - 1)
            outputMatrix(matrix, y, n);

        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        xVect[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - matrix[i][k] * xVect[k];
    }
    return xVect;
}
void initialisationConst(double** matrix, double* y) {
    matrix[0][0] = 2.30; matrix[0][1] = 3.50; matrix[0][2] = 1.70;
    matrix[1][0] = 5.70; matrix[1][1] = -2.70; matrix[1][2] = 2.30;
    matrix[2][0] = -0.80; matrix[2][1] = 5.30; matrix[2][2] = -1.80;
    y[0] = -6.49; y[1] = 19.20; y[2] = -5.09;
}
void residVectorMax(double* resVect, int n) {
    double max = resVect[0];
    for (int i = 1; i < n; ++i) {
        if (max < resVect[i]) {
            max = resVect[i];
        }
    }
    cout << "maximal elem resid vect =" << max << endl;
}
void residVector(double** matrix, double* resVect, double* bVect, double* xVect, int n) {
    double* a;
    a = new double[n];
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2; j++) {
            resVect[i] += matrix[i][j] * xVect[j];
        }
    }
    cout << "vector nevyazki" << endl;
    for (int i = 0; i <= 2; i++) {
        a[i] = resVect[i] - bVect[i];
        cout << a[i] << endl;
    }
    residVectorMax(a, n);
    delete[]a;
}