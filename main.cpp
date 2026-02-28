
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
      out << std::setprecision(4) << x << " ";
    }
    return out;
  }

private:
  static int student_coef(size_t n, double alpha) {
    boost::math::students_t dist(n);
    /*
     *                     alpha
     <----|--------------------------------------|---->
        -t_stat              0                 +t_stat
    */
    double t_stat = boost::math::quantile(dist, 0.5 + alpha / 2);
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
    double res;
    double mean = get_median();

    for (auto x : data) {
      res += x - mean;
    }
    return res;
  }
  double sum_of_deviations_from_arithmetic_mean_square() const {
    double res;
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
    return sqrt(sum_of_deviations_from_arithmetic_mean_square() /
                (data.size() - 1) / data.size());
  }
  // median + delta
  std::pair<double, double> confidence_interval_for_mean(double alpha) const {
    double delta =
        student_coef(data.size(), alpha) * standard_error_of_average();
    return {get_median(), delta};
  }
};

std::ostream &lab1_table(std::ostream &out, const statistic &s) {
  auto d = s.deviations_from_arithmetic_mean();
  out << "N;Source;d;d_squre";
  for (size_t i = 0; i < s.data.size(); ++i) {
    out << i << ";" << s.data[i] << ";" << d[i] << ";" << d[i] * d[i];
  }
  return out;
}

int main() {
  std::cout << "lab: 1, file: data_raw.csv" << std::endl;

  statistic s;

  s.from_file("data_raw.csv");
  std::cout << "data: " << s << std::endl;

  auto median = s.get_median();
  std::cout << "median " << median << std::endl;

  auto [m, d] = s.confidence_interval_for_mean(0.98);
  std::cout << "X is " << m - d << " " << m + d << std::endl;

  return 0;
}
