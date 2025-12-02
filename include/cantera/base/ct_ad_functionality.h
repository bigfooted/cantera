#pragma once
#include "ct_typedefs.h"
#include "fmt.h"

#if AD_ENABLED


template <> struct fmt::formatter<CanteraDouble> : fmt::formatter<double> {

  auto format(const CanteraDouble& p, format_context& ctx) const {
    return fmt::formatter<double>::format(p.getValue(), ctx);
  }
};

CanteraDoublePassive castToPassiveDouble(CanteraDouble const& value) {
  return value.getValue();
}

std::vector<CanteraDouble> castToDoubleVector(std::vector<CanteraDoublePassive> const& vec) {
  std::vector<CanteraDouble> res(vec.size());
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = vec[i];
  }
  return res;
}

std::vector<std::vector<CanteraDouble>> castToDoubleVectorVector(std::vector<std::vector<CanteraDoublePassive>> const& vec) {
  std::vector<std::vector<CanteraDouble>> res(vec.size(), std::vector<CanteraDouble>(0));
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = castToDoubleVector(vec[i]);
  }
  return res;
}

std::vector<CanteraDoublePassive> castToPassiveDoubleVector(std::vector<CanteraDouble> const& vec) {
  std::vector<CanteraDoublePassive> res(vec.size());
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = vec[i].getValue();
  }
  return res;
}

std::vector<std::vector<CanteraDoublePassive>> castToPassiveDoubleVectorVector(std::vector<std::vector<CanteraDouble>> const& vec) {
  std::vector<std::vector<CanteraDoublePassive>> res(vec.size(), std::vector<CanteraDoublePassive>(0));
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = castToPassiveDoubleVector(vec[i]);
  }
  return res;
}
#else

CanteraDoublePassive castToPassiveDouble(CanteraDouble const& value) {
  return value;
}

std::vector<CanteraDouble> const& castToDoubleVector(std::vector<CanteraDoublePassive> const& vec) {
  return vec;
}

std::vector<std::vector<CanteraDouble>> const& castToDoubleVectorVector(std::vector<std::vector<CanteraDoublePassive>> const& vec) {
  return vec;
}
std::vector<CanteraDoublePassive> const& castToPassiveDoubleVector(std::vector<CanteraDouble> const& vec) {
  return vec;
}
std::vector<std::vector<CanteraDoublePassive>> const& castToPassiveDoubleVectorVector(std::vector<std::vector<CanteraDouble>> const& vec) {
  return vec;
}
#endif
