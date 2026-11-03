///|
#derive(Debug)
pub(all) suberror MatrixError {
  InvalidModulus(Int64)
  InvalidDimensions(Int32, Int32)
  InvalidEntryCount(Int32, Int32)
  DimensionMismatch(Int32, Int32)
  NotSquare(Int32, Int32)
  Singular
  NegativeExponent(Int64)
}

///|
pub extend MatrixError with @debug.Debug::{to_repr}

///|
/// Exact dense linear algebra over F_p. The bound p <= 2^31-1 makes every
/// product of two residues fit in signed Int64 before modular reduction.
#derive(Eq, Debug)
pub struct Matrix {
  rows : Int32
  columns : Int32
  modulus : Int64
  entries : FixedArray[Int64]
}

///|
pub extend Matrix with Eq::{equal, not_equal}

///|
pub extend Matrix with @debug.Debug::{to_repr}

///|
pub fn is_prime(modulus : Int64) -> Bool {
  if modulus < 2 || modulus > 2147483647L {
    return false
  }
  if modulus % 2L == 0L {
    return modulus == 2L
  }
  let mut divisor : Int64 = 3L
  while divisor <= modulus / divisor {
    if modulus % divisor == 0L {
      return false
    }
    divisor += 2L
  }
  true
}

///|
pub fn residue(value : Int64, modulus : Int64) -> Int64 {
  let remainder = value % modulus
  if remainder < 0L {
    remainder + modulus
  } else {
    remainder
  }
}

///|
/// Constant-memory extended Euclid inversion in a prime field.
pub fn inverse(value : Int64, modulus : Int64) -> Int64 raise MatrixError {
  if !is_prime(modulus) {
    raise MatrixError::InvalidModulus(modulus)
  }
  let mut remainder = residue(value, modulus)
  if remainder == 0L {
    raise MatrixError::Singular
  }
  let mut old_remainder = modulus
  let mut coefficient : Int64 = 1L
  let mut old_coefficient : Int64 = 0L
  while remainder != 0L {
    let quotient = old_remainder / remainder
    let next = old_remainder - quotient * remainder
    old_remainder = remainder
    remainder = next
    let next_coefficient = old_coefficient - quotient * coefficient
    old_coefficient = coefficient
    coefficient = next_coefficient
  }
  residue(old_coefficient, modulus)
}

///|
pub fn Matrix::new(
  rows : Int32,
  columns : Int32,
  modulus : Int64,
  entries : Array[Int64],
) -> Matrix raise MatrixError {
  if !is_prime(modulus) {
    raise MatrixError::InvalidModulus(modulus)
  }
  if rows < 0 || columns < 0 || (rows != 0 && columns > 2147483647 / rows) {
    raise MatrixError::InvalidDimensions(rows, columns)
  }
  if entries.length() != rows * columns {
    raise MatrixError::InvalidEntryCount(entries.length(), rows * columns)
  }
  {
    rows,
    columns,
    modulus,
    entries: FixedArray::makei(entries.length(), fn(i) {
      residue(entries[i], modulus)
    }),
  }
}

///|
pub fn Matrix::zero(
  rows : Int32,
  columns : Int32,
  modulus : Int64,
) -> Matrix raise MatrixError {
  Matrix::new(rows, columns, modulus, Array::make(rows * columns, 0L))
}

///|
pub fn Matrix::identity(
  size : Int32,
  modulus : Int64,
) -> Matrix raise MatrixError {
  Matrix::new(
    size,
    size,
    modulus,
    Array::makei(size * size, fn(i) {
      if i / size == i % size {
        1L
      } else {
        0L
      }
    }),
  )
}

///|
pub fn Matrix::rows(self : Matrix) -> Int32 {
  self.rows
}

///|
pub fn Matrix::columns(self : Matrix) -> Int32 {
  self.columns
}

///|
pub fn Matrix::modulus(self : Matrix) -> Int64 {
  self.modulus
}

///|
pub fn Matrix::get(
  self : Matrix,
  row : Int32,
  column : Int32,
) -> Int64 raise MatrixError {
  if row < 0 || row >= self.rows || column < 0 || column >= self.columns {
    raise MatrixError::InvalidDimensions(row, column)
  }
  self.entries[row * self.columns + column]
}

///|
pub fn Matrix::to_array(self : Matrix) -> Array[Array[Int64]] {
  Array::makei(self.rows, fn(row) {
    Array::makei(self.columns, fn(column) {
      self.entries[row * self.columns + column]
    })
  })
}

///|
fn swap_rows(matrix : Matrix, left : Int32, right : Int32) -> Unit {
  if left == right {
    return
  }
  for column in 0..<matrix.columns {
    let a = left * matrix.columns + column
    let b = right * matrix.columns + column
    let saved = matrix.entries[a]
    matrix.entries[a] = matrix.entries[b]
    matrix.entries[b] = saved
  }
}

///|
#derive(Debug)
pub struct RowEchelon {
  matrix : Matrix
  pivot_columns : Array[Int32]
}

///|
pub extend RowEchelon with @debug.Debug::{to_repr}

///|
/// Gauss–Jordan elimination with normalized pivots and elimination above
/// and below each pivot. The input matrix is left unchanged.
pub fn Matrix::rref(self : Matrix) -> RowEchelon raise MatrixError {
  let matrix = { ..self, entries: self.entries.copy(), }
  let pivots : Array[Int32] = []
  let mut row = 0
  for column in 0..<matrix.columns {
    if row == matrix.rows {
      break
    }
    let mut pivot = row
    while pivot < matrix.rows &&
          matrix.entries[pivot * matrix.columns + column] == 0L {
      pivot += 1
    }
    if pivot == matrix.rows {
      continue
    }
    swap_rows(matrix, row, pivot)
    let scale = inverse(
      matrix.entries[row * matrix.columns + column],
      matrix.modulus,
    ) catch {
      _ => raise MatrixError::Singular
    }
    for j in 0..<matrix.columns {
      let index = row * matrix.columns + j
      matrix.entries[index] = matrix.entries[index] * scale % matrix.modulus
    }
    for other in 0..<matrix.rows {
      if other == row {
        continue
      }
      let factor = matrix.entries[other * matrix.columns + column]
      if factor == 0L {
        continue
      }
      for j in 0..<matrix.columns {
        let index = other * matrix.columns + j
        matrix.entries[index] = residue(
          matrix.entries[index] -
          factor * matrix.entries[row * matrix.columns + j],
          matrix.modulus,
        )
      }
    }
    pivots.push(column)
    row += 1
  }
  { matrix, pivot_columns: pivots, }
}

///|
pub fn RowEchelon::matrix(self : RowEchelon) -> Matrix {
  self.matrix
}

///|
pub fn RowEchelon::pivot_columns(self : RowEchelon) -> Array[Int32] {
  self.pivot_columns.copy()
}

///|
pub fn RowEchelon::rank(self : RowEchelon) -> Int32 {
  self.pivot_columns.length()
}

///|
#derive(Debug)
pub struct AffineSolution {
  particular : FixedArray[Int64]
  basis : Array[FixedArray[Int64]]
}

///|
pub extend AffineSolution with @debug.Debug::{to_repr}

///|
pub fn AffineSolution::particular(self : AffineSolution) -> FixedArray[Int64] {
  self.particular.copy()
}

///|
pub fn AffineSolution::basis(self : AffineSolution) -> Array[FixedArray[Int64]] {
  self.basis.map(fn(vector) { vector.copy() })
}

///|
/// Return one solution and a basis of the homogeneous solution space.
/// None means inconsistency, not an empty homogeneous solution space.
pub fn Matrix::solve(
  self : Matrix,
  right : Array[Int64],
) -> AffineSolution? raise MatrixError {
  if right.length() != self.rows {
    raise MatrixError::DimensionMismatch(right.length(), self.rows)
  }
  let augmented = Array::makei(self.rows * (self.columns + 1), fn(i) {
    let row = i / (self.columns + 1)
    let column = i % (self.columns + 1)
    if column == self.columns {
      right[row]
    } else {
      self.entries[row * self.columns + column]
    }
  })
  let echelon = (try! Matrix::new(
    self.rows,
    self.columns + 1,
    self.modulus,
    augmented,
  )).rref()
  let reduced = echelon.matrix
  let pivots = echelon.pivot_columns
  for row in echelon.rank()..<self.rows {
    if reduced.entries[row * reduced.columns + self.columns] != 0L {
      return None
    }
  }
  for column in pivots {
    if column == self.columns {
      return None
    }
  }
  let is_pivot = FixedArray::make(self.columns, false)
  let particular = FixedArray::make(self.columns, 0L)
  for row, column in pivots {
    is_pivot[column] = true
    particular[column] = reduced.entries[row * reduced.columns + self.columns]
  }
  let basis : Array[FixedArray[Int64]] = []
  for free in 0..<self.columns {
    if is_pivot[free] {
      continue
    }
    let vector = FixedArray::make(self.columns, 0L)
    vector[free] = 1L
    for row, column in pivots {
      vector[column] = residue(
        -reduced.entries[row * reduced.columns + free],
        self.modulus,
      )
    }
    basis.push(vector)
  }
  Some({ particular, basis, })
}

///|
pub fn Matrix::kernel(
  self : Matrix,
) -> Array[FixedArray[Int64]] raise MatrixError {
  let solution = try! self.solve(Array::make(self.rows, 0L))
  match solution {
    Some(value) => value.basis
    None => raise MatrixError::Singular
  }
}

///|
pub fn Matrix::determinant(self : Matrix) -> Int64 raise MatrixError {
  if self.rows != self.columns {
    raise MatrixError::NotSquare(self.rows, self.columns)
  }
  let matrix = { ..self, entries: self.entries.copy(), }
  let mut result : Int64 = 1L
  for column in 0..<self.columns {
    let mut pivot = column
    while pivot < self.rows &&
          matrix.entries[pivot * self.columns + column] == 0L {
      pivot += 1
    }
    if pivot == self.rows {
      return 0L
    }
    if pivot != column {
      swap_rows(matrix, pivot, column)
      result = residue(-result, self.modulus)
    }
    let value = matrix.entries[column * self.columns + column]
    result = result * value % self.modulus
    let scale = inverse(value, self.modulus)
    for row in (column + 1)..<self.rows {
      let factor = matrix.entries[row * self.columns + column] *
        scale %
        self.modulus
      for j in column..<self.columns {
        let index = row * self.columns + j
        matrix.entries[index] = residue(
          matrix.entries[index] -
          factor * matrix.entries[column * self.columns + j],
          self.modulus,
        )
      }
    }
  }
  result
}

///|
pub fn Matrix::multiply(
  self : Matrix,
  right : Matrix,
) -> Matrix raise MatrixError {
  if self.modulus != right.modulus {
    raise MatrixError::InvalidModulus(right.modulus)
  }
  if self.columns != right.rows {
    raise MatrixError::DimensionMismatch(self.columns, right.rows)
  }
  let entries = Array::make(self.rows * right.columns, 0L)
  for row in 0..<self.rows {
    for inner in 0..<self.columns {
      let factor = self.entries[row * self.columns + inner]
      for column in 0..<right.columns {
        let index = row * right.columns + column
        entries[index] = (
            entries[index] +
            factor * right.entries[inner * right.columns + column]
          ) %
          self.modulus
      }
    }
  }
  Matrix::new(self.rows, right.columns, self.modulus, entries)
}

///|
pub fn Matrix::apply(
  self : Matrix,
  vector : FixedArray[Int64],
) -> FixedArray[Int64] raise MatrixError {
  if vector.length() != self.columns {
    raise MatrixError::DimensionMismatch(vector.length(), self.columns)
  }
  FixedArray::makei(self.rows, fn(row) {
    let mut sum : Int64 = 0L
    for column, value in vector {
      sum = (sum + self.entries[row * self.columns + column] * value) %
        self.modulus
    }
    sum
  })
}

///|
pub fn Matrix::power(
  self : Matrix,
  exponent : Int64,
) -> Matrix raise MatrixError {
  if self.rows != self.columns {
    raise MatrixError::NotSquare(self.rows, self.columns)
  }
  if exponent < 0L {
    raise MatrixError::NegativeExponent(exponent)
  }
  let mut result = try! Matrix::identity(self.rows, self.modulus)
  let mut base = self
  let mut remaining = exponent
  while remaining > 0L {
    if remaining % 2L == 1L {
      result = try! result.multiply(base)
    }
    remaining /= 2L
    if remaining > 0L {
      base = try! base.multiply(base)
    }
  }
  result
}

///|
pub fn Matrix::subtract(
  self : Matrix,
  right : Matrix,
) -> Matrix raise MatrixError {
  if self.rows != right.rows || self.columns != right.columns {
    raise MatrixError::DimensionMismatch(self.columns, right.columns)
  }
  if self.modulus != right.modulus {
    raise MatrixError::InvalidModulus(right.modulus)
  }
  Matrix::new(
    self.rows,
    self.columns,
    self.modulus,
    Array::makei(self.entries.length(), fn(i) {
      residue(self.entries[i] - right.entries[i], self.modulus)
    }),
  )
}
