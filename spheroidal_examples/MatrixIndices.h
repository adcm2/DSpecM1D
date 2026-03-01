#ifndef MATRIX_INDICES_GUARD_H
#define MATRIX_INDICES_GUARD_H

namespace Spheroidal {

class MatrixIndices {
  // this class is used to index the matrix
  // it is used to convert the indices of the matrix to a single index
  // the ordering is via element, node, l, m
  // l is the degree of the spherical harmonic, m is the order

public:
  // constructors
  MatrixIndices() {};
  MatrixIndices(int, int, int, int, int);

  // number of independent spherical harmonics
  int N_SH() const;

  // full index
  std::size_t mat_index(int, int, int, int);

  // total size
  std::size_t mat_size() const;

  // get the APPROXIMATE number of non-zero elements in the matrix
  std::size_t mat_nonzero() const;

  // matrix if just at a given l
  std::size_t mat_index_l(int idxelem, int idxnode) const {
    // this gives the index for a given element, node and l
    std::size_t tmpidx = (_NQ - 1) * (idxelem - _el) + idxnode;
    if (_el == 0) {
      --tmpidx;
    }
    return tmpidx;
  };

  std::size_t mat_size_l() const {
    // this gives the size of the matrix for a given l
    std::size_t tmpidx = (_NQ - 1) * _en + 1;
    if (_el == 0) {
      --tmpidx;
    }
    return tmpidx;
  };

  std::size_t mat_nonzero_l() const {
    // this gives APPROXIMATELY the number of non-zero elements in the matrix
    // for a given l
    std::size_t tmpidx = _NQ * _NQ * _en;
    if (_el == 0) {
      tmpidx -= _NQ;
    }
    return tmpidx;
  };

private:
  int _lMax, _lMin, _linlen, _NQ, _elemlen, _el, _eu, _en;
  std::size_t _numlen, _totnonzero;
  // _numlen is the number of elements in the matrix
};

// constructor for matrix indices
MatrixIndices::MatrixIndices(int el, int eu, int lMin, int lMax, int NQ)
    : _lMax{lMax}, _lMin{lMin}, _NQ{NQ}, _el{el}, _eu{eu}, _en{eu - el} {
  // total number of lm from lmin to lmax
  _linlen = (lMax + 1 + lMin) * (lMax + 1 - lMin);

  // number of points in element excluding the shared layer
  _elemlen = (NQ - 1) * _linlen;

  // number of elements in the matrix
  _numlen = _elemlen * _en + _linlen;

  // if the first element is the first element of the mesh, we need to remove
  //  the first element from the number of elements in the matrix
  if (el == 0) {
    _numlen -= _linlen;
  }

  // total number of non-zero elements in the matrix
  _totnonzero = NQ * NQ * _en * _linlen;

  // if the first element is the first element of the mesh, we need to
  // remove the first element from the number of non-zero elements in the
  // matrix
  if (el == 0) {
    _totnonzero -= 2 * NQ * _linlen;
  }
  //   std::cout << "MatrixIndices: "
  //             << "lMin: " << _lMin << ", lMax: " << _lMax
  //             << ", linlen: " << _linlen << ", elemlen: " << _elemlen
  //             << ", numlen: " << _numlen << ", totnonzero: " << _totnonzero
  //             << "\n";
};

// number of independent spherical harmonics
int
MatrixIndices::N_SH() const {
  return this->_linlen;
};   // number of independent spherical harmonics

// full index
// this is the index of the matrix element corresponding to the element, node,
// spherical harmonic degree and order
std::size_t
MatrixIndices::mat_index(int idxelem, int idxnode, int idxl, int idxm) {
  // initialise the index at start of the element
  // add the index of the node, then the index of the spherical harmonic
  // degree and order
  // the index of the spherical harmonic is given by the formula:
  // l(l+1) + m - lmin^2
  std::size_t tmpidx = this->_elemlen * (idxelem - _el);
  //   std::cout << "\ntmpidx: " << tmpidx << ", elemlen: " << this->_elemlen
  //             << ", idxelem - el: " << idxelem - _el << "\n";
  tmpidx += this->_linlen * idxnode;
  //   std::cout << "tmpidx: " << tmpidx << "\n";
  tmpidx += idxl * idxl + idxl + idxm - this->_lMin * this->_lMin;
  //   std::cout << "tmpidx: " << tmpidx << "\n";

  // if the first element is the first element of the mesh, we need to
  //  remove the first element from the index
  if (_el == 0) {
    tmpidx -= this->_linlen;
    //   std::cout << "tmpidx: " << tmpidx << "\n";
  }
  return tmpidx;
};

// total size
std::size_t
MatrixIndices::mat_size() const {
  return this->_numlen;
};

// get the number of non-zero elements in the matrix
std::size_t
MatrixIndices::mat_nonzero() const {
  return _totnonzero;
};

}   // namespace Spheroidal
#endif