##include <iostream>
#include <cmath>

// Подынтегральная функция для одной переменной
double f(double x) {
    return std::sqrt(x + std::pow(x, 3));
}

// Подынтегральная функция для двух переменных
double f2D(double x, double y) {
    // Замените на вашу функцию двух переменных
    return std::sqrt(x + std::pow(x, 3)) * y;
}

// Метод трапеций для одной переменной с критерием завершения
double trapezoidalIntegration(double a, double b, double epsilon) {
    int n = 1;
    double h = b - a;
    double integral_h, integral_h_half;

    do {
        h /= 2.0;

        // Вычисление интеграла для h
        integral_h = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i) {
            integral_h += f(a + i * h);
        }
        integral_h *= h;

        // Вычисление интеграла для h / 2
        integral_h_half = 0.5 * (f(a) + f(b));
        for (int i = 1; i < 2 * n; ++i) {
            integral_h_half += f(a + 0.5 * i * h);
        }
        integral_h_half *= 0.5 * h;

        n *= 2;

    } while (std::abs(integral_h_half - integral_h) > 3 * epsilon);

    return integral_h_half;
}

// Метод Симпсона для одной переменной с критерием завершения
double simpsonIntegration(double a, double b, double epsilon) {
    int n = 1;
    double h = b - a;
    double integral_h, integral_h_half;

    do {
        h /= 2.0;

        // Вычисление интеграла для h
        integral_h = f(a) + f(b);

        for (int i = 1; n < n; ++i) {
            integral_h += (i % 2 == 0) ? 2 * f(a + i * h) : 4 * f(a + i * h);
        }

        integral_h = (h / 3.0) * integral_h;

        // Вычисление интеграла для h / 2
        integral_h_half = f(a) + f(b);

        for (int i = 1; i < 2 * n; ++i) {
            integral_h_half += (i % 2 == 0) ? 2 * f(a + 0.5 * i * h) : 4 * f(a + 0.5 * i * h);
        }

        integral_h_half = (0.5 * h / 3.0) * integral_h_half;

        n *= 2;

    } while (std::abs(integral_h_half - integral_h) > 15 * epsilon);

    return integral_h_half;
}

// Метод Симпсона для двух переменных с критерием завершения
double simpsonIntegration2D(double a, double b, double c, double d, double epsilon) {
    int nx = 1, ny = 1;
    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    double integral_h, integral_h_half;

    do {
        hx /= 2.0;
        hy /= 2.0;

        // Вычисление интеграла для h
        integral_h = f2D(a, c) + f2D(a, d) + f2D(b, c) + f2D(b, d);

        for (int i = 1; i < nx; ++i) {
            integral_h += 2 * (f2D(a + i * hx, c) + f2D(a + i * hx, d) + f2D(b + i * hx, c) + f2D(b + i * hx, d));
        }

        for (int j = 1; j < ny; ++j) {
            integral_h += 2 * (f2D(a, c + j * hy) + f2D(a, d + j * hy) + f2D(b, c + j * hy) + f2D(b, d + j * hy));
        }

        for (int i = 1; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                integral_h += 4 * f2D(a + i * hx, c + j * hy);
            }
        }

        integral_h = (hx * hy / 9.0) * integral_h;

        // Вычисление интеграла для h / 2
        integral_h_half = f2D(a, c) + f2D(a, d) + f2D(b, c) + f2D(b, d);

        for (int i = 1; i < 2 * nx; ++i) {
            integral_h_half += 2 * (f2D(a + 0.5 * i * hx, c) + f2D(a + 0.5 * i * hx, d) + f2D(b + 0.5 * i * hx, c) + f2D(b + 0.5 * i * hx, d));
        }

        for (int j = 1; j < 2 * ny; ++j) {
            integral_h_half += 2 * (f2D(a, c + 0.5 * j * hy) + f2D(a, d + 0.5 * j * hy) + f2D(b, c + 0.5 * j * hy) + f2D(b, d + 0.5 * j * hy));
        }

        for (int i = 1; i < 2 * nx; ++i) {
            for (int j = 1; j < 2 * ny; ++j) {
                integral_h_half += 4 * f2D(a + 0.5 * i * hx, c + 0.5 * j * hy);
            }
        }

        integral_h_half = (0.5 * hx * 0.5 * hy / 9.0) * integral_h_half;

        nx *= 2;
        ny *= 2;

    } while (std::abs(integral_h_half - integral_h) > 15 * epsilon);

    return integral_h_half;
}

int main() {
    // Интервалы
    double a = 0.6, b = 1.724;
    double c = -1.0, d = 1.0; // Пример значений для интервала по второй переменной

    // Заданные значения для epsilon
    double epsilon1 = 1e-4;
    double epsilon2 = 1e-5;

    // Решение задачи методом Симпсона для двух переменных с критерием завершения
    double simpson2D_result1 = simpsonIntegration2D(a, b, c, d, epsilon1);
    double simpson2D_result2 = simpsonIntegration2D(a, b, c, d, epsilon2);

    std::cout << "Simpson 2D Integration Result (epsilon = 1e-4): " << simpson2D_result1 << std::endl;
    std::cout << "Simpson 2D Integration Result (epsilon = 1e-5): " << simpson2D_result2 << std::endl;

    return 0;
}

