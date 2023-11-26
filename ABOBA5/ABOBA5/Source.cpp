#include <iostream>
#include <cmath>

// Подынтегральная функция
double f(double x) {
    return std::sqrt(x + std::pow(x, 3));
}

// Метод трапеций с критерием завершения
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

// Метод Симпсона с критерием завершения
double simpsonIntegration(double a, double b, double epsilon) {
    int n = 1;
    double h = b - a;
    double integral_h, integral_h_half;

    do {
        h /= 2.0;

        // Вычисление интеграла для h
        integral_h = f(a) + f(b);

        for (int i = 1; i < n; ++i) {
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

int main() {
    // Интервал
    double a = 0.6, b = 1.724;

    // Заданные значения для epsilon
    double epsilon1 = 1e-4;
    double epsilon2 = 1e-5;

    // Решение задачи методом трапеций с критерием завершения
    double trap_result1 = trapezoidalIntegration(a, b, epsilon1);
    double trap_result2 = trapezoidalIntegration(a, b, epsilon2);

    std::cout << "Trapezoidal Integration Result (epsilon = 1e-4): " << trap_result1 << std::endl;
    std::cout << "Trapezoidal Integration Result (epsilon = 1e-5): " << trap_result2 << std::endl;

    // Решение задачи методом Симпсона с критерием завершения
    double simpson_result1 = simpsonIntegration(a, b, epsilon1);
    double simpson_result2 = simpsonIntegration(a, b, epsilon2);

    std::cout << "Simpson Integration Result (epsilon = 1e-4): " << simpson_result1 << std::endl;
    std::cout << "Simpson Integration Result (epsilon = 1e-5): " << simpson_result2 << std::endl;

    return 0;
}
