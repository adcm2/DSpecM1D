#ifndef DSPECM1D_REFERENCE_SERIES_IO_H
#define DSPECM1D_REFERENCE_SERIES_IO_H

#include <algorithm>
#include <string>
#include <Eigen/Core>
#include "ReadYSpec.h"
#include "ReadSpecnm.h"
#include "ReadMineos.h"

namespace DSpecM {

struct ReferenceTimeSeries {
  Eigen::MatrixXd yspec_t;
  Eigen::MatrixXd mineos_t;
};

inline Eigen::MatrixXd
loadYSpecTimeSeries(const std::string &yspecPath, int ncols) {
  YSPECREADER::DataColumns yspecData(yspecPath);
  const std::size_t maxCols = std::max(0, ncols);
  Eigen::MatrixXd yspecT =
      Eigen::MatrixXd::Zero(3, static_cast<Eigen::Index>(maxCols));

  const std::size_t count = std::min(maxCols, yspecData.getColumn1().size());
  for (std::size_t idx = 0; idx < count; ++idx) {
    yspecT(0, static_cast<Eigen::Index>(idx)) = yspecData.getColumn2()[idx];
    yspecT(1, static_cast<Eigen::Index>(idx)) = yspecData.getColumn3()[idx];
    yspecT(2, static_cast<Eigen::Index>(idx)) = yspecData.getColumn4()[idx];
  }

  return yspecT;
}

inline Eigen::MatrixXd
loadSpecnmTimeSeries(const std::string &specnmPath, int ncols) {
  SPECNMREADER::DataColumns specnmData(specnmPath);
  const std::size_t maxCols = std::max(0, ncols);
  Eigen::MatrixXd specnmT =
      Eigen::MatrixXd::Zero(3, static_cast<Eigen::Index>(maxCols));

  const std::size_t count = std::min(maxCols, specnmData.getColumn1().size());
  for (std::size_t idx = 0; idx < count; ++idx) {
    specnmT(0, static_cast<Eigen::Index>(idx)) = specnmData.getColumn2()[idx];
    specnmT(1, static_cast<Eigen::Index>(idx)) = specnmData.getColumn3()[idx];
    specnmT(2, static_cast<Eigen::Index>(idx)) = specnmData.getColumn4()[idx];
  }

  return specnmT;
}

inline Eigen::MatrixXd
loadMineosTimeSeries(const std::string &mineosZPath,
                     const std::string &mineosNPath,
                     const std::string &mineosEPath, int ncols,
                     double scale = 1e-9) {
  MINEOSREADER::DataColumns mineosZ(mineosZPath);
  MINEOSREADER::DataColumns mineosN(mineosNPath);
  MINEOSREADER::DataColumns mineosE(mineosEPath);

  const std::size_t maxCols = std::max(0, ncols);
  Eigen::MatrixXd mineosT =
      Eigen::MatrixXd::Zero(3, static_cast<Eigen::Index>(maxCols));

  const std::size_t count =
      std::min({maxCols, mineosZ.getColumn1().size(),
                mineosN.getColumn1().size(), mineosE.getColumn1().size()});
  for (std::size_t idx = 0; idx < count; ++idx) {
    mineosT(0, static_cast<Eigen::Index>(idx)) =
        mineosZ.getColumn2()[idx] * scale;
    mineosT(1, static_cast<Eigen::Index>(idx)) =
        mineosN.getColumn2()[idx] * scale;
    mineosT(2, static_cast<Eigen::Index>(idx)) =
        mineosE.getColumn2()[idx] * scale;
  }

  return mineosT;
}

inline ReferenceTimeSeries
loadReferenceTimeSeries(const std::string &yspecPath,
                        const std::string &mineosZPath,
                        const std::string &mineosNPath,
                        const std::string &mineosEPath, int ncols,
                        double mineosScale = 1e-9) {
  ReferenceTimeSeries out;
  out.yspec_t = loadYSpecTimeSeries(yspecPath, ncols);
  out.mineos_t = loadMineosTimeSeries(mineosZPath, mineosNPath, mineosEPath,
                                      ncols, mineosScale);
  return out;
}

}   // namespace DSpecM

#endif   // DSPECM1D_REFERENCE_SERIES_IO_H
