#pragma once
#include "cpp11.hpp"
#include "iter.h"
#include "max.h"
#include "min.h"
#include <thread>
#include <tuple>

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
sens(const std::vector<double> &e, const std::vector<double> &lwr,
     const std::vector<double> &upr, const std::vector<double> &lwr_iter,
     const std::vector<double> &upr_iter, const std::vector<int> &grp,
     const std::vector<bool> &a, const double epsilon, const double eta,
     const size_t m, const size_t chunk_size, const int num_threads);

cpp11::list sens_(const cpp11::doubles &e, const cpp11::doubles &lwr,
                  const cpp11::doubles &upr, const cpp11::doubles &lwr_iter,
                  const cpp11::doubles &upr_iter, const cpp11::integers &grp,
                  const cpp11::logicals &a, const double epsilon,
                  const double eta, const int m, const int chunk_size,
                  const int num_threads);
