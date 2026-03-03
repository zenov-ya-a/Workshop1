
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

#include <boost/math/distributions/students_t.hpp>

struct statistic {
  std::vector<double> data;

  void from_file(const std::string &file_name) {
    std::ifstream in(file_name);
    if (in.is_open()) {
      char buf[100];
      in.getline(buf, 100, ';');
      in.getline(buf, 100);

      size_t number;
      char delimer;
      double measure;
      for (; in >> number >> delimer >> measure;) {
        data.push_back(measure);
      }
    }
  }
  friend std::ostream &operator<<(std::ostream &out, const statistic &s) {
    for (auto x : s.data) {
      out << x << " ";
    }
    return out;
  }

private:
  static double student_coef(size_t n, double alpha) {
    if (n <= 1) {
      return 0.0;
    }
    boost::math::students_t dist(n - 1);
    /*
     *                     alpha
     <----|--------------------------------------|---->
        -t_stat              0                 +t_stat
    */
    auto t_stat = boost::math::quantile(dist, 0.5 + alpha / 2.0);
    return t_stat;
  }

public:
  double get_median() const {
    double res;
    for (auto x : data) {
      res += x;
    }
    if (!data.empty()) {
      res /= data.size();
    }
    return res;
  }

  std::vector<double> deviations_from_arithmetic_mean() const {
    std::vector<double> res;
    double mean = get_median();

    res.reserve(data.size());
    for (auto x : data) {
      res.push_back(x - mean);
    }
    return res;
  }
  std::vector<double> deviations_from_arithmetic_mean_square() const {
    std::vector<double> res;
    double mean = get_median();

    res.reserve(data.size());
    for (auto x : data) {
      res.push_back((x - mean) * (x - mean));
    }
    return res;
  }

  double sum_of_deviations_from_arithmetic_mean() const {
    double res = 0;
    double mean = get_median();

    for (auto x : data) {
      res += x - mean;
    }
    return res;
  }
  double sum_of_deviations_from_arithmetic_mean_square() const {
    double res = 0;
    double mean = get_median();

    for (auto x : data) {
      res += (x - mean) * (x - mean);
    }
    return res;
  }
  double biased_assessment() const {
    if (data.empty()) {
      return 0;
    }
    return sum_of_deviations_from_arithmetic_mean_square() / data.size();
  }

  double unbiased_assessment() const {
    if (data.size() < 2) {
      return 0;
    }
    return sum_of_deviations_from_arithmetic_mean_square() / (data.size() - 1);
  }
  double standard_deviation() { return sqrt(unbiased_assessment()); }
  double standard_error_of_average() const {
    if (data.size() < 2) {
      return 0;
    }
    auto n = data.size();
    return sqrt(sum_of_deviations_from_arithmetic_mean_square() /
                ((n - 1) * n));
  }
  // median + delta
  std::pair<double, double> confidence_interval_for_mean(double alpha) const {
    double delta =
        student_coef(data.size(), alpha) * standard_error_of_average();
    std::cout << "student: " << student_coef(data.size(), alpha) << std::endl
              << "eror of average: " << standard_error_of_average()
              << std::endl;
    return {get_median(), delta};
  }
  // gistogram, start, delta
  std::tuple<std::vector<size_t>, double, double>
  build_gistorgam(size_t step_count) const {
    auto sorted = data;
    std::sort(sorted.begin(), sorted.end());

    auto start = sorted.front();
    auto end = sorted.back();
    auto delta = (end - start) / step_count;
    /* |N1|N2|N3|N4|
     * 0  1  2  3
     * [ )[ )[ )[ ]
     * s 1d 2d 3d e
     */
    size_t j = 0;
    std::vector<size_t> res(step_count, 0);
    for (size_t i = 0; i < sorted.size() - 1; ++i) {
      if (sorted[i] >= start + delta * (j + 1)) {
        j += 1;
      }
      res[j] += 1;
    }
    res.back() += 1;

    return {res, start, delta};
  }
};

std::ostream &lab1_table(std::ostream &out, const statistic &s) {
  auto d = s.deviations_from_arithmetic_mean();
  out << "N;Source;d;d_squre" << std::endl;
  for (size_t i = 0; i < s.data.size(); ++i) {
    out << i + 1 << ";" << s.data[i] << ";" << d[i] << ";" << d[i] * d[i]
        << std::endl;
  }
  return out;
}

int main() {
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "lab: 1, file: data_raw.csv" << std::endl;

  statistic s;

  s.from_file("data_raw.csv");
  for (auto &x : s.data) {
    x *= 1000;
  }
  std::cout << "data (Hz): " << s << std::endl;

  auto median = s.get_median();
  std::cout << "median " << median << " Hz" << std::endl;

  auto [m, d] = s.confidence_interval_for_mean(0.98);
  std::cout << "X is " << m - d << " Hz" << " " << m + d << " Hz" << std::endl;
  std::cout << "d(50, 0,98) = " << d << " Hz" << std::endl;

  std::ofstream fout("lab_table.csv");
  lab1_table(fout, s);

  {
    std::ofstream fout("gistogram.csv");
    auto [res, start, delta] = s.build_gistorgam(20);
    fout << "X;dN" << std::endl;
    for (int i = 0; i < res.size(); ++i) {
      fout << start + delta * (i + 0.5) << ";"
           << (double)res[i] / (double)s.data.size() << std::endl;
    }
  }
  double max = 0;
  for (auto x : s.data) {
    max = std::max(max, 1 + 5 * 10e-7 * x * 1000);
  }
  std::cout << "max tool error = " << max / 1000 << " Hz" << std::endl;

  return 0;
}
