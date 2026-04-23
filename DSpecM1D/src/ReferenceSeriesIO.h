#ifndef DSPECM1D_REFERENCE_SERIES_IO_H
#define DSPECM1D_REFERENCE_SERIES_IO_H

#include <algorithm>
#include <string>
#include <Eigen/Core>
#include "ReadYSpec.h"
#include "ReadSpecnm.h"
#include "ReadMineos.h"

namespace DSpecM {

/**
 * @brief Bundle of reference time-series arrays used for comparison workflows.
 */
struct ReferenceTimeSeries {
  Eigen::MatrixXd yspecTime;
  Eigen::MatrixXd mineosTime;
};

/**
 * @brief Loads a three-component YSpec time series into a fixed-size matrix.
 *
 * @param yspecPath Path to the YSpec four-column text file.
 * @param ncols Number of output samples to allocate.
 * @return Matrix with rows `Z`, `N`, `E` and `ncols` columns, padded with
 * zeros when the file is shorter than requested.
 */
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

/**
 * @brief Loads a three-component SpecNM time series from semicolon-delimited
 * text.
 *
 * @param specnmPath Path to the SpecNM four-column text file.
 * @param ncols Number of output samples to allocate.
 * @return Matrix with rows `Z`, `N`, `E` and `ncols` columns, padded with
 * zeros when the file is shorter than requested.
 */
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

/**
 * @brief Loads three Mineos component files into a single matrix.
 *
 * @param mineosZPath Path to the vertical component file.
 * @param mineosNPath Path to the north component file.
 * @param mineosEPath Path to the east component file.
 * @param ncols Number of output samples to allocate.
 * @param scale Multiplicative scale factor applied to all loaded amplitudes.
 * @return Matrix with rows `Z`, `N`, `E` and `ncols` columns, padded with
 * zeros when the files are shorter than requested.
 */
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

/**
 * @brief Convenience loader that returns both YSpec and Mineos reference
 * traces.
 *
 * @param yspecPath Path to the YSpec time-series file.
 * @param mineosZPath Path to the Mineos vertical component file.
 * @param mineosNPath Path to the Mineos north component file.
 * @param mineosEPath Path to the Mineos east component file.
 * @param ncols Number of output samples to allocate.
 * @param mineosScale Multiplicative scale factor applied to Mineos amplitudes.
 * @return Populated `ReferenceTimeSeries` bundle.
 */
inline ReferenceTimeSeries
loadReferenceTimeSeries(const std::string &yspecPath,
                        const std::string &mineosZPath,
                        const std::string &mineosNPath,
                        const std::string &mineosEPath, int ncols,
                        double mineosScale = 1e-9) {
  ReferenceTimeSeries out;
  out.yspecTime = loadYSpecTimeSeries(yspecPath, ncols);
  out.mineosTime = loadMineosTimeSeries(mineosZPath, mineosNPath, mineosEPath,
                                        ncols, mineosScale);
  return out;
}

}   // namespace DSpecM

#endif   // DSPECM1D_REFERENCE_SERIES_IO_H
