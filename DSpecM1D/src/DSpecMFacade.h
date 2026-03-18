#ifndef DSPECM1D_FACADE_H
#define DSPECM1D_FACADE_H

#include "Timer.h"
#include "ReadStation.h"
#include "InputParser.h"
#include "InputParametersNew.h"
#include "SRInfo.h"
#include "SourceInfo.h"
#include "StartElement.h"
#include "SpecHelpers.h"
#include "SpectraRunContext.h"
#include "ReadYSpec.h"
#include "ReadMineos.h"
#include "ReferenceSeriesIO.h"
#include "SignalFiltering.h"
#include "OutputWriters.h"
#include "SEM/SEM.h"
#include "FullSpec.h"

namespace DSpecM {

// Namespace-level aliases for all major library modules.
namespace Full1D = ::Full1D;
namespace SparseSpec = ::SPARSESPEC;
namespace Source = ::SourceInfo;
namespace YSpecReader = ::YSPECREADER;
namespace MineosReader = ::MINEOSREADER;
namespace SpectralTools = ::SpectralTools;

using InputParameters = ::InputParameters;
using InputParametersNew = ::InputParametersNew;
using SRInfo = ::SRInfo;
using EarthquakeCMT = SourceInfo::EarthquakeCMT;
using Timer = ::Timer;
using SiteChanEntry = ::SiteChanEntry;
using SpectraRunContext = SPARSESPEC::SpectraRunContext;

using SEM = Full1D::SEM;
using SparseFSpec = SPARSESPEC::SparseFSpec;

using YSpecDataColumns = YSPECREADER::DataColumns;
using MineosDataColumns = MINEOSREADER::DataColumns;

using ::read_full_sitechan_file;
using SPARSESPEC::attenFactor;
using SPARSESPEC::factorizeOrCompute;
using SPARSESPEC::outputFactor;
using SPARSESPEC::resolveModeFlags;

}   // namespace DSpecM

#endif   // DSPECM1D_FACADE_H
