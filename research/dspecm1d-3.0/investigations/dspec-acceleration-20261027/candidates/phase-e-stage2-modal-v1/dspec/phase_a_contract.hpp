#ifndef DSPECM1D_PHASE_A_CONTRACT_HPP
#define DSPECM1D_PHASE_A_CONTRACT_HPP

#include <algorithm>
#include <array>
#include <iterator>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace phasea {

using Complex = std::complex<double>;

enum class FieldRole : std::int32_t { ForwardSource = 0, AdjointDual = 1 };
enum class Family : std::uint8_t { Radial = 1, Spheroidal = 2, Toroidal = 4 };
inline constexpr std::size_t SourceOrderSlots = 5;
inline constexpr std::size_t AdjointOrderSlots = 3;
inline constexpr std::size_t ModalOrderSlots = SourceOrderSlots + AdjointOrderSlots;

struct ComplexPair {
  double real = 0.0;
  double imag = 0.0;
  Complex value() const noexcept { return {real, imag}; }
  static ComplexPair from(Complex z) noexcept { return {z.real(), z.imag()}; }
};

/** One native side and degree, padded to five source orders and three adjoint
 * orders. Slots 0..4 are forward source m=-2..2; slots 5..7 are adjoint-dual
 * m=-1..1 for the single receiver 0. Each mode is [U,U',V,V',W,W']; primes
 * are d/dr_nd. Spheroidal U/V and toroidal W have already been merged at
 * identical (l,m). */
struct RadialSideDegreeRecord {
  std::int32_t bin = 0;
  std::int32_t receiver = 0;
  std::int32_t degree = 0;
  std::int32_t layer = -1;
  std::int32_t element = -1;
  std::int32_t localNode = -1;
  std::uint8_t fluid = 0;
  std::uint8_t adjointFrameIsIsotropic = 0;
  std::uint8_t sourceFrameSet = 0;
  std::uint8_t adjointFrameSet = 0;
  double radiusNormalized = 0.0;
  double radiusMetres = 0.0;
  double elementLowerNormalized = 0.0;
  double elementUpperNormalized = 0.0;
  double elementJacobianNormalized = 0.0;
  std::array<double,9> sourceFrameAxes{};
  std::array<double,9> adjointFrameAxes{};
  std::array<std::uint8_t, ModalOrderSlots> familyMask{};
  std::array<std::uint8_t, ModalOrderSlots> sourceLiftApplied{};
  std::array<ComplexPair, ModalOrderSlots * 6> radial{};
};

inline constexpr int sourceOrder(std::size_t slot) {
  return static_cast<int>(slot) - 2;
}
inline constexpr int adjointOrder(std::size_t slot) {
  return static_cast<int>(slot - SourceOrderSlots) - 1;
}

/** Accepted DS packed/tensor reference at one native side and section point. */
struct ReferenceFieldRecord {
  std::int32_t bin = 0;
  FieldRole role = FieldRole::ForwardSource;
  std::int32_t receiver = -1;
  std::int32_t section = 0;
  std::int32_t layer = -1;
  std::int32_t element = -1;
  std::int32_t localNode = -1;
  std::int32_t sectionPoint = -1;
  std::uint8_t fluid = 0;
  std::uint8_t reserved[7]{};
  double radiusNormalized = 0.0;
  double thetaRadians = 0.0;
  double phiRadians = 0.0;
  std::array<ComplexPair, 6> strain{};
};

struct PacketHeader {
  char magic[8] = {'D','S','P','H','A','R','1','\0'};
  std::uint32_t schema = 1;
  std::uint32_t recordBytes = 0;
  std::uint64_t recordCount = 0;
  std::int32_t bin = 0;
  std::uint32_t reserved = 0;
  double zReal = 0.0;
  double zImag = 0.0;
  double lengthScaleMetres = 0.0;
  double timeScaleSeconds = 0.0;
};

static_assert(std::is_trivially_copyable_v<ComplexPair>);
static_assert(std::is_trivially_copyable_v<RadialSideDegreeRecord>);
static_assert(std::is_trivially_copyable_v<ReferenceFieldRecord>);
static_assert(std::is_trivially_copyable_v<PacketHeader>);

template<class Record>
void writePacket(const std::filesystem::path &path, PacketHeader header,
                 const std::vector<Record> &records) {
  static_assert(std::is_trivially_copyable_v<Record>);
  if (sizeof(Record)>std::numeric_limits<std::uint32_t>::max())
    throw std::overflow_error("Phase A packet record size exceeds schema limit");
  if (std::filesystem::exists(path))
    throw std::runtime_error("Phase A packet output already exists: " + path.string());
  if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
    if (records.size() > std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("Phase A packet record count overflows");
  }
  if (records.size() > static_cast<std::size_t>(
          std::numeric_limits<std::streamsize>::max() / sizeof(Record)))
    throw std::overflow_error("Phase A packet byte count exceeds streamsize");
  header.recordCount = static_cast<std::uint64_t>(records.size());
  header.recordBytes = static_cast<std::uint32_t>(sizeof(Record));
  std::ofstream out(path, std::ios::binary | std::ios::out);
  if (!out) throw std::runtime_error("cannot create Phase A packet: " + path.string());
  out.write(reinterpret_cast<const char *>(&header), sizeof(header));
  if (!records.empty())
    out.write(reinterpret_cast<const char *>(records.data()),
              static_cast<std::streamsize>(records.size() * sizeof(Record)));
  out.flush();
  if (!out) throw std::runtime_error("failed writing Phase A packet: " + path.string());
}

template<class Record>
std::vector<Record> readPacket(const std::filesystem::path &path,
                               PacketHeader &header) {
  static_assert(std::is_trivially_copyable_v<Record>);
  std::ifstream in(path, std::ios::binary | std::ios::in);
  if (!in) throw std::runtime_error("cannot open Phase A packet: " + path.string());
  in.read(reinterpret_cast<char *>(&header), sizeof(header));
  const char expected[8] = {'D','S','P','H','A','R','1','\0'};
  if (!in || !std::equal(std::begin(header.magic), std::end(header.magic), expected) ||
      header.schema != 1 || header.recordBytes != sizeof(Record))
    throw std::runtime_error("invalid Phase A packet header: " + path.string());
  if (header.recordCount > std::numeric_limits<std::size_t>::max() / sizeof(Record))
    throw std::overflow_error("Phase A packet size overflows address space");
  if (header.recordCount > static_cast<std::uint64_t>(
          std::numeric_limits<std::streamsize>::max() / sizeof(Record)))
    throw std::overflow_error("Phase A packet byte count exceeds streamsize");
  std::vector<Record> records(static_cast<std::size_t>(header.recordCount));
  if (!records.empty())
    in.read(reinterpret_cast<char *>(records.data()),
            static_cast<std::streamsize>(records.size() * sizeof(Record)));
  if (!in || in.peek() != std::char_traits<char>::eof())
    throw std::runtime_error("truncated or trailing bytes in Phase A packet: " + path.string());
  return records;
}

} // namespace phasea

#endif
