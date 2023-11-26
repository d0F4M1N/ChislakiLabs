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
    double integral_prev, integral_current = 0;

    do {
        integral_prev = integral_current;
        h /= 2.0;
        double sum = 0.0;

        for (int i = 0; i < n; ++i) {
            sum += f(a + (2 * i + 1) * h);
        }

        integral_current = 0.5 * integral_prev + h * sum;
        n *= 2;

    } while (std::abs(integral_current - integral_prev) > epsilon);

    return integral_current;
}

// Метод Симпсона с критерием завершения
double simpsonIntegration(double a, double b, double epsilon) {
    int n = 1;
    double h = b - a;
    double integral_prev, integral_current = 0;

    do {
        integral_prev = integral_current;
        h /= 2.0;
        double sum = 0.0;

        for (int i = 0; i < n; ++i) {
            sum += (i % 2 == 0) ? 2 * f(a + i * h) : 4 * f(a + i * h);
        }

        integral_current = (h / 3.0) * (f(a) + sum + f(b));
        n *= 2;

    } while (std::abs(integral_current - integral_prev) > epsilon);

    return integral_current;
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
