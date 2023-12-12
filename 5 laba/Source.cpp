#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
double f(double x) {
    return pow(x + pow(x, 3),0.5);
}
double f2(double x, double y)
{
    return sin(x + y);
}
double integralTrapezia(double a, double b, int n,double h)
{
    double sum = h / 2 * (f(a) + f(b));
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x) * h;
    }
    return sum;
}
double integralSimpson(double a, double b, int n,double h)
{
    double sum = h / 3 * (f(a) + f(b));
    for (int i = 1; i < 2 * n; i += 2) {
        double x = a + i * h;
        sum += 4 * f(x) * h / 3;
    }
    for (int i = 2; i < 2 * n - 1; i += 2) {
        double x = a + i * h;
        sum += 2 * f(x) * h / 3;
    }
    return sum;
}
double trapezoidalIntegration(double a, double b, double eps) {
    int n = 2;
    double h = (b - a) / n;
    while (integralTrapezia(a, b, n * 2, h / 2) - integralTrapezia(a, b, n, h) > 3 * eps)
    {
        h /= 2;
        n *= 2;
    }
    return integralTrapezia(a, b, n, h);
}

double simpsonIntegration(double a, double b, double eps) {
    int n = 2;
    double h = (b - a) / (2 * n);
    while (integralSimpson(a, b, n * 2, h / 2) - integralSimpson(a, b, n, h) > 15 * eps)
    {
        h /= 2;
        n *= 2;
    }
    return integralSimpson(a, b, n, h);
}
double cubeSimpsonIntegration(double a, double b, double c, double d, int n, int m)
{
    double hx = (b - a) / (2 * n);
    double hy = (d - c) / (2 * m);
    double integral = 0;

    for (int i = 0; i < n - 1; i++)
    {
        double sum = 0;
        for (int j = 0; j < n - 1; j++)
        {
            sum += f2(a + hx * 2 * i, c + hy * 2 * j) + 4 * f2(a + hx * (2 * i + 1), c + hy * 2 * j) +
                f2(a + hx * (2 * i + 2), c + hy * 2 * j) + 4 * f2(a + hx * 2 * i, c + hy * (2 * j + 1)) +
                16 * f2(a + hx * (2 * i + 1), c + hy * (2 * j + 1)) + 4 * f2(a + hx * (2 * i + 2), c + hy * (2 * j + 1)) +
                f2(a + hx * 2 * i, c + hy * (2 * j + 2)) + 4 * f2(a + hx * (2 * i + 1), c + hy * (2 * j + 2)) +
                f2(a + hx * (2 * i + 2), c + hy * (2 * j + 2));
        }
        integral += sum;
    }
    return hx * hy / 9 * integral;
}
int main() {
    double a = 0.6;
    double b = 1.724;
    double n = 1000;
    double m = 1000;
    cout << "e = 10e-4" << endl;
    cout << "Trapezoidal: " << trapezoidalIntegration(a, b, 10e-4) << endl;
    cout << "Simpson result: " << simpsonIntegration(a, b, 10e-4) << endl;
    cout << "e = 10e-5" << endl;
    cout << "Trapezoidal: " << trapezoidalIntegration(a, b, 10e-5) << endl;
    cout << "Simpson result: " << simpsonIntegration(a, b, 10e-5) << endl;
    cout << "Cube simpson: " << cubeSimpsonIntegration(0, 3.14 / 2, 0, 3.14 / 4, n, m) << endl;
}
