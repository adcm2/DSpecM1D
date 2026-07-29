#ifndef DSPECM1D_LOVE_NUMBERS_DOF_MAP_H
#define DSPECM1D_LOVE_NUMBERS_DOF_MAP_H

#include <optional>
#include <stdexcept>
#include <vector>

#include <DSpecM1D/EarthMesh>

namespace DSpecM1D::LoveNumbers::detail {

enum class Field { U, V, P };

class LoveDofMap {
public:
  LoveDofMap(const EarthMesh::RadialMesh &mesh, int degree) {
    if (degree < 0) {
      throw std::invalid_argument("degree must be non-negative.");
    }

    indices_.resize(
        mesh.NE(), std::vector<NodeIndices>(mesh.NN()));

    for (int element = 0; element < mesh.NE(); ++element) {
      for (int node = 0; node < mesh.NN(); ++node) {
        auto &current = indices_[element][node];
        const bool shared_node = element > 0 && node == 0;
        const bool surface_node =
            element == mesh.NE() - 1 && node == mesh.NN() - 1;

        if (degree == 0) {
          if (shared_node) {
            current = indices_[element - 1].back();
          } else {
            current.u = size_++;
            current.p = size_++;
          }
          continue;
        }

        if (shared_node) {
          current.p = indices_[element - 1].back().p;
        }

        if (mesh.IsSolid(element)) {
          if (shared_node && mesh.IsSolid(element - 1)) {
            current.u = indices_[element - 1].back().u;
            current.v = indices_[element - 1].back().v;
          } else {
            current.u = size_++;
            current.v = size_++;
          }
        }

        if (!shared_node && !(degree == 1 && surface_node)) {
          current.p = size_++;
        }
      }
    }
  }

  std::optional<int> index(Field field, int element, int node) const {
    const auto &indices = indices_.at(element).at(node);
    switch (field) {
    case Field::U:
      return indices.u;
    case Field::V:
      return indices.v;
    case Field::P:
      return indices.p;
    }
    return std::nullopt;
  }

  int size() const noexcept { return size_; }

private:
  struct NodeIndices {
    std::optional<int> u;
    std::optional<int> v;
    std::optional<int> p;
  };

  std::vector<std::vector<NodeIndices>> indices_;
  int size_ = 0;
};

}   // namespace DSpecM1D::LoveNumbers::detail

#endif   // DSPECM1D_LOVE_NUMBERS_DOF_MAP_H
