#ifndef FILE_MATRIX
#define FILE_MATRIX

#include "matexpr.hpp"
#include "vector.hpp"
#include <algorithm>


namespace nanoblas
{
  
  // enum ORDERING { RowMajor, ColMajor };

  template <typename T, ORDERING ORD>
  class MatrixView : public MatExpr<MatrixView<T,ORD>>
  {
  protected:
    T* m_data;
    size_t m_rows, m_cols;
    size_t m_dist;

    size_t index(size_t i, size_t j) const
    {
      return (ORD==RowMajor) ? i*m_dist+j : j*m_dist+i;
    }
  public:
    MatrixView() = default;
    MatrixView(const MatrixView &) = default;
        
    MatrixView (size_t rows, size_t cols, T* data)
      : m_data(data), m_rows(rows), m_cols(cols)
    {
      m_dist = (ORD==RowMajor) ? m_cols : m_rows;
    }
    
    MatrixView (size_t rows, size_t cols, size_t dist, T* data)
      : m_data(data), m_rows(rows), m_cols(cols), m_dist(dist) { }
    
    template <typename TB, ORDERING ORD2>
    MatrixView (const MatrixView<TB,ORD2>& m2)
      : m_data(m2.data()), m_rows(m2.rows()), m_cols(m2.cols()), m_dist(m2.dist()) { }


    MatrixView& operator= (const MatrixView& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }
    
    template <typename TB>
    MatrixView& operator= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }
        
    MatrixView& operator= (T scal)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = scal;
      return *this;
    }
        
    T* data() const { return m_data; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    size_t dist() const { return m_dist; }
    auto shape() const { return std::array<size_t,2>{m_rows, m_cols}; }

    T& operator()(size_t i, size_t j) 
    { 
      return m_data[index(i,j)];
    }
        
    const T& operator()(size_t i, size_t j) const 
    {
      return m_data[index(i,j)];
    }


    template <typename TB>
    MatrixView& operator+= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) += m2(i,j);
      return *this;
    }
    
    template <typename TB>
    MatrixView& operator-= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) -= m2(i,j);
      return *this;
    }

    MatrixView& operator*= (T scal)
    {
      for (size_t i = 0; i < rows(); i++)
        for (size_t j = 0; j < cols(); j++)
          (*this)(i,j) *= scal;
      return *this;
    }


    auto row(size_t i) const 
    {
      if constexpr (ORD == RowMajor)
        return VectorView<T>(m_cols, m_data + i*m_dist);
      else
        return VectorView<T,size_t>(m_cols, m_dist, m_data + i);
    } 

    auto col(size_t j) const 
    {
      if constexpr (ORD == RowMajor)
        return VectorView<T,size_t>(m_rows, m_dist, m_data + j);
      else
        return VectorView<T>(m_rows, m_data + j*m_dist);
    }

    auto diag() const 
    {
      return VectorView<T,size_t>(std::min(m_rows, m_cols), m_dist+1, m_data);
    } 

    auto rows(size_t first, size_t next) const 
    {
      return MatrixView<T, ORD>(next - first, m_cols, m_dist, m_data+index(first,0));
    }

    auto cols(size_t first, size_t next) const 
    {
      return MatrixView<T, ORD>(m_rows, next - first, m_dist, m_data+index(0,first));
    }
  };


  template <typename T, ORDERING ORD>
  auto trans (MatrixView<T,ORD> mat)
  {
    if constexpr (ORD==RowMajor)
      return MatrixView<T,ColMajor>(mat.cols(), mat.rows(), mat.dist(), mat.data());
    else
      return MatrixView<T,RowMajor>(mat.cols(), mat.rows(), mat.dist(), mat.data());
  }
  
  template <typename T=double, ORDERING ORD=RowMajor>
  class Matrix : public MatrixView<T,ORD>
  {
    typedef MatrixView<T,ORD> BASE;
    using BASE::m_cols;
    using BASE::m_data;
    using BASE::m_rows;

  public:
    Matrix (size_t rows, size_t cols)
      : BASE(rows, cols, new T[rows*cols]) { }
          
    Matrix (const Matrix& m2)
      : BASE(m2.rows(), m2.cols(), new T[m2.rows()*m2.cols()])
    {
      *this = m2;
    }
          
    Matrix (std::initializer_list<std::initializer_list<T>> list)
      : BASE(list.size(), list.begin()->size(), new T[list.size()*list.begin()->size()])
    {
      size_t i = 0;
      for (auto row : list)
        {
          size_t j = 0;
          for (auto val : row)
            {
              (*this)(i,j) = val;
              j++;
            }
          i++;
        }
    } 

    ~Matrix() { delete[] m_data; }

    using BASE::operator=;
    Matrix& operator= (const Matrix& m2)
    {
      if (this != &m2)
        {
          assert(m_rows==m2.m_rows && m_cols==m2.m_cols);
          BASE::operator=(m2);
        }
      return *this;
    }
                    
  };

// Annahme: MatrixView-Schnittstelle wie im Skript:
// - size_t rows() const;
// - size_t cols() const;
// - size_t dist() const;      // leading dimension
// - T& operator()(size_t i, size_t j);
// - MatrixView rows(i1,i2), cols(j1,j2), etc.

template <typename T = double, ORDERING ORD = ColMajor>
void addMatMat (MatrixView<T,ORD> A,
                MatrixView<T,ORD> B,
                MatrixView<T,ORD> C)
{
  constexpr size_t BH = 96;
  constexpr size_t BW = 96;

  static_assert(ORD == ColMajor, "addMatMat ist für ColMajor ausgelegt");

  alignas(64) T memBA[BH*BW];

  for (size_t i1 = 0; i1 < A.rows(); i1 += BH)
    for (size_t j1 = 0; j1 < A.cols(); j1 += BW) {
      size_t i2 = std::min(A.rows(), i1 + BH);
      size_t j2 = std::min(A.cols(), j1 + BW);

      MatrixView<T,ORD> Ablock(i2 - i1, j2 - j1, BW, memBA);
      Ablock = A.rows(i1, i2).cols(j1, j2);

      // B-Block: (j1..j2) Zeilen, alle Spalten
      // C-Block: (i1..i2) Zeilen, alle Spalten
      addMatMat2(Ablock, B.rows(j1, j2), C.rows(i1, i2));
    }
}


/// Mikro-Kernel für einen H×W-Block von C.
/// A: Zeiger auf A(i,0) (H Zeilen, K Spalten, ColMajor, leading dimension a_dist)
/// B: Zeiger auf B(0,j) (K Zeilen, W Spalten, ColMajor, leading dimension b_dist)
/// C: Zeiger auf C(i,j) (H Zeilen, W Spalten, ColMajor, leading dimension c_dist)
template <size_t H, size_t W, typename T = double>
void AddMatMatKernel(size_t K,
                     const T* A, size_t a_dist,
                     const T* B, size_t b_dist,
                     T*       C, size_t c_dist)
{
  // Akkumulatoren im Register/Stack
  T acc[H][W];

  // C-Block laden
  for (size_t h = 0; h < H; ++h)
    for (size_t w = 0; w < W; ++w)
      acc[h][w] = C[h + w * c_dist];

  // Spaltenzeiger für B vorbereiten: B(0, j+w)
  const T* b_cols[W];
  for (size_t w = 0; w < W; ++w)
    b_cols[w] = B + w * b_dist;

  // K-Schleife
  for (size_t k = 0; k < K; ++k) {
    const T* a_col_k = A + k * a_dist;  // Zeiger auf Spalte k von A (ab Zeile i)
    T a_vals[H];
    for (size_t h = 0; h < H; ++h)
      a_vals[h] = a_col_k[h];

    for (size_t w = 0; w < W; ++w) {
      T b_kw = b_cols[w][k];            // B(k, j+w)
      for (size_t h = 0; h < H; ++h)
        acc[h][w] += a_vals[h] * b_kw;
    }
  }

  // Zurück nach C schreiben
  for (size_t h = 0; h < H; ++h)
    for (size_t w = 0; w < W; ++w)
      C[h + w * c_dist] = acc[h][w];
}


template <typename T = double, ORDERING ORD = ColMajor>
void addMatMat2 (MatrixView<T,ORD> A,
                 MatrixView<T,ORD> B,
                 MatrixView<T,ORD> C)
{
  constexpr size_t H = 4;
  constexpr size_t W = 12;

  static_assert(ORD == ColMajor, "addMatMat2 ist für ColMajor ausgelegt");

  const size_t M = C.rows();
  const size_t N = C.cols();
  const size_t K = A.cols();   // = B.rows()

  // Vollständige H×W-Blöcke
  size_t j = 0;
  for (; j + W <= N; j += W) {
    size_t i = 0;
    for (; i + H <= M; i += H) {
      AddMatMatKernel<H, W, T>(
        K,
        &A(i, 0), A.dist(),
        &B(0, j), B.dist(),
        &C(i, j), C.dist()
      );
    }

    // Restzeilen für diese W-Spalten
    for (; i < M; ++i) {
      for (size_t k = 0; k < K; ++k) {
        T aik = A(i, k);
        for (size_t jj = 0; jj < W; ++jj)
          C(i, j + jj) += aik * B(k, j + jj);
      }
    }
  }

  // Restspalten (Spalten < W)
  for (; j < N; ++j) {
    for (size_t i = 0; i < M; ++i) {
      T sum = C(i, j);
      for (size_t k = 0; k < K; ++k)
        sum += A(i, k) * B(k, j);
      C(i, j) = sum;
    }
  }
}



}


#endif
