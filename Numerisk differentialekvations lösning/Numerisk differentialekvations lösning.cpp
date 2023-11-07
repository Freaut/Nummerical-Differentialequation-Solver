#include <iostream>
#include <cmath>
#include <functional>
#include <fstream>
#include <locale>
#include <chrono>
#include "matplotlibcpp.h"

using namespace matplotlibcpp;

struct Comma final : std::numpunct<char>
{
    char do_decimal_point() const override { return ','; }
};

struct point {
    double x;
    double y;

    point(double _x, double _y) {
        x = _x;
        y = _y;
    }
};

struct result {
    std::vector<point> points;
    double y;

    result(std::vector<point> _points, double _y) {
        points = _points;
        y = _y;
    }
};

result solve_f_euler(std::function<double(double, double)> f, double x0, double y0, double t, double x1) {
    double x = x0;
    double y = y0;

    std::vector<point> values;
	
    while (x < x1) {
        double dy_dx = f(x, y);
        y = y + t * dy_dx;
        x += t;
        values.push_back(point(x, y));
    }

    return result(values, y);
}

result solve_b_euler(std::function<double(double, double)> f, double x0, double y0, double t, double x1) {
    double y_prev = y0;
    double y_next = y_prev + t * f(x0 + t, y_prev + t * f(x0, y_prev));
    double x = x0 + t;
	
    std::vector<point> values;
	
    while (x < x1) {
        y_prev = y_next;
        y_next = y_prev + t * f(x + t, y_prev + t * f(x, y_prev));
        x += t;
        values.push_back(point(x, y_next));
    }

    return result(values, y_next);
}

result runge_kutta(std::function<double(double, double)> f, double y0, double x0, double t, double x1) {
    double y = y0;
    double x = x0;

    std::vector<point> values;
	
    while (x < x1) {
        double k1 = t * f(x, y);
        double k2 = t * f(x + t / 2, y + k1 / 2);
        double k3 = t * f(x + t / 2, y + k2 / 2);
        double k4 = t * f(x + t, y + k3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x += t;
        values.push_back(point(x, y));
    }

    return result(values, y);
}


int main()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    double y0 = 0;
    double x0 = 0;
    double t = 0.01;
    double x1 = 1;
    bool write_to_file = false;
    const char* filename = "rk_points.txt";

    std::function<double(double, double)> f = [](double x, double y) {
        //return x * x;
        return ((x * x) - (2*x) + 1);
        //return (pow(Py_MATH_E, -y)) * sin(2 * Py_MATH_PI * x);
    };

    auto rk_t1 = high_resolution_clock::now();
    result r_k = runge_kutta(f, x0, y0, t, x1);
    auto rk_t2 = high_resolution_clock::now();

    auto fe_t1 = high_resolution_clock::now();
    result r_y1     = solve_f_euler(f, x0, y0, t, x1);
    auto fe_t2 = high_resolution_clock::now();

    auto be_t1 = high_resolution_clock::now();
    result r_y2     = solve_b_euler(f, x0, y0, t, x1);
    auto be_t2 = high_resolution_clock::now();


    duration<double, std::milli> rk_ms_double = rk_t2 - rk_t1;
    duration<double, std::milli> fe_ms_double = fe_t2 - fe_t1;
    duration<double, std::milli> be_ms_double = be_t2 - be_t1;
    std::cout << "Runge-kutta:  " << rk_ms_double.count() << "ms\n";
    std::cout << "Feuler:       " << fe_ms_double.count() << "ms\n";
    std::cout << "Beuler:       " << be_ms_double.count() << "ms\n";
	
    if (write_to_file) {
        std::ofstream outfile(filename);
        outfile.imbue(std::locale(std::locale::classic(), new Comma));

        for (const point& p : r_k.points) {
            outfile << p.x << " " << p.y << "\n";
        }

        outfile.close();
    }

    std::cout << "Runge-kutta x = " << x1 << " y = " << r_k.y << "\n";
    std::cout << "F euler x = " << x1 << " y = " << r_y1.y << "\n";
    std::cout << "B euler x = " << x1 << " y = " << r_y2.y << "\n";
}