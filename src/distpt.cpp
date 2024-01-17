#include "distpt.h"

///////////////////////////// DistPt Implementation ////////////////////////////

/*
 * @brief Default constructor for DistPt.
 * @return A DistPt with value 0 and type PtType::lwr.
 */
DistPt::DistPt() : val(0), type(PtType::lwr) {}

/*
 * @brief Explicit constructor for DistPt.
 * @param v Value of the DistPt.
 * @param t Type of the DistPt.
 * @return A DistPt with value v and type t.
 */
DistPt::DistPt(double v, PtType t) : val(v), type(t) {}

/*
 * @brief Copy constructor for DistPt.
 * @param other DistPt to be copied.
 * @return A DistPt with value and type copied from other.
 */
DistPt::DistPt(const DistPt &other) : val(other.val), type(other.type) {}

/*
 * @brief Comparison operator for DistPts.
 * @param a first DistPt.
 * @param b second DistPt.
 * @return True if a.val < b.val or a.val == b.val and a.type < b.type.
 */
bool operator<(const DistPt &a, const DistPt &b) {
  return a.val == b.val ? a.type < b.type : a.val < b.val;
}

/*
 * @brief Default constructor for PtVector.
 * @return An empty PtVector.
 */
PtVector::PtVector(){};

/*
 * @brief Constructor for vectors of DistPts.
 * @param e Estimated risk.
 * @param l Lower bound on risk.
 * @param u Upper bound on risk.
 * @pre e.size() == l.size() == u.size();
 * @pre e.size() > 0;
 * @pre 0 <= e[i] <= l[i] <= u[i] <= 1 for all i.
 * @return A vector of DistPts generated from the input values sorted in
 *         increasing order of value.
 */
PtVector::PtVector(const std::vector<double> &e, const std::vector<double> &l,
                   const std::vector<double> &u) {
  // Add all points to the vector
  for (std::size_t i = 0; i < e.size(); ++i) {
    pts.push_back({e[i], PtType::est});
    pts.push_back({l[i], PtType::lwr});
    pts.push_back({u[i], PtType::upr});
  }

  // Sort the vector
  std::sort(pts.begin(), pts.end());
}

/*
 * @brief Subscript operator for PtVector.
 * @param index Index of the DistPt.
 * @pre 0 <= index < pts.size().
 * @return The DistPt at the given index.
 */
DistPt &PtVector::operator[](std::size_t index) { return pts[index]; }

/*
 * @brief Subscript operator for PtVector.
 * @param index Index of the DistPt.
 * @pre 0 <= index < pts.size().
 * @return The DistPt at the given index.
 */
const DistPt &PtVector::operator[](std::size_t index) const {
  return pts[index];
}

/*
 * @brief Get the size of the PtVector.
 * @return The size of the PtVector.
 */
const std::size_t PtVector::size() const { return pts.size(); }
