#pragma once
#include <algorithm>
#include <vector>

enum class PtType { lwr, est, upr };

struct DistPt {
  double val;  // Value of the point
  PtType type; // Type of the point

  DistPt();
  DistPt(double v, PtType t);
  DistPt(const DistPt &other);
};

class PtVector {
public:
  PtVector(const std::vector<double> &e, const std::vector<double> &l,
           const std::vector<double> &u);
  PtVector();

  DistPt &operator[](std::size_t index);
  const DistPt &operator[](std::size_t index) const;
  const std::size_t size() const;

private:
  std::vector<DistPt> pts;
};
