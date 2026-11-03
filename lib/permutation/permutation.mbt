///|
/// Errors at the public finite-action boundary.
#derive(Debug)
pub(all) suberror DomainError {
  NegativeDegree(Int32)
  PointOutOfRange(Int32, Int32)
  NotBijection
  DegreeMismatch(Int32, Int32)
  InvalidCycle
}

///|
pub extend DomainError with @debug.Debug::{to_repr}

///|
/// A bijection of {0, ..., degree-1}. The representation is private so that
/// every value reaching another package has passed the bijectivity check.
#derive(Eq, Debug, Hash)
pub struct Permutation {
  images : FixedArray[Int32]
}

///|
pub extend Permutation with Eq::{equal, not_equal}

///|
pub extend Permutation with Hash::{hash, hash_combine}

///|
pub extend Permutation with @debug.Debug::{to_repr}

///|
pub fn Permutation::identity(degree : Int32) -> Permutation raise DomainError {
  if degree < 0 {
    raise DomainError::NegativeDegree(degree)
  }
  { images: FixedArray::makei(degree, fn(i) { i }), }
}

///|
pub fn Permutation::from_images(
  images : Array[Int32],
) -> Permutation raise DomainError {
  let n = images.length()
  let seen = FixedArray::make(n, false)
  for image in images {
    if image < 0 || image >= n {
      raise DomainError::PointOutOfRange(image, n)
    }
    if seen[image] {
      raise DomainError::NotBijection
    }
    seen[image] = true
  }
  { images: FixedArray::from_array(images[:]), }
}

///|
pub fn Permutation::degree(self : Permutation) -> Int32 {
  self.images.length()
}

///|
/// Defensive copy: changing the result cannot invalidate a permutation.
pub fn Permutation::to_images(self : Permutation) -> FixedArray[Int32] {
  self.images.copy()
}

///|
pub fn Permutation::image(
  self : Permutation,
  point : Int32,
) -> Int32 raise DomainError {
  let n = self.degree()
  if point < 0 || point >= n {
    raise DomainError::PointOutOfRange(point, n)
  }
  self.images[point]
}

///|
/// Composition self ∘ right; the right permutation acts first.
pub fn Permutation::compose(
  self : Permutation,
  right : Permutation,
) -> Permutation raise DomainError {
  let n = self.degree()
  if n != right.degree() {
    raise DomainError::DegreeMismatch(n, right.degree())
  }
  { images: FixedArray::makei(n, fn(i) { self.images[right.images[i]] }), }
}

///|
pub fn Permutation::inverse(self : Permutation) -> Permutation {
  let result = FixedArray::make(self.degree(), 0)
  for i, image in self.images {
    result[image] = i
  }
  { images: result, }
}

///|
/// Canonical cycles: each starts at its least point; fixed points are omitted.
pub fn Permutation::cycles(self : Permutation) -> Array[Array[Int32]] {
  let seen = FixedArray::make(self.degree(), false)
  let result : Array[Array[Int32]] = []
  for start in 0..<self.degree() {
    if seen[start] {
      continue
    }
    let cycle : Array[Int32] = []
    let mut point = start
    while !seen[point] {
      seen[point] = true
      cycle.push(point)
      point = self.images[point]
    }
    if cycle.length() > 1 {
      result.push(cycle)
    }
  }
  result
}

///|
/// Cycles are applied from left to right. Repeated points within one cycle
/// are rejected; points may occur in different cycles.
pub fn Permutation::from_cycles(
  degree : Int32,
  cycles : Array[Array[Int32]],
) -> Permutation raise DomainError {
  let mut result = Permutation::identity(degree)
  for cycle in cycles {
    if cycle.length() < 2 {
      continue
    }
    let images = FixedArray::makei(degree, fn(i) { i })
    let seen = FixedArray::make(degree, false)
    for point in cycle {
      if point < 0 || point >= degree || seen[point] {
        raise DomainError::InvalidCycle
      }
      seen[point] = true
    }
    for i, point in cycle {
      images[point] = cycle[(i + 1) % cycle.length()]
    }
    result = { images, }.compose(result)
  }
  result
}

///|
pub fn Permutation::is_identity(self : Permutation) -> Bool {
  for i, image in self.images {
    if i != image {
      return false
    }
  }
  true
}

///|
pub fn Permutation::fixed_points(self : Permutation) -> Array[Int32] {
  let result : Array[Int32] = []
  for i, image in self.images {
    if i == image {
      result.push(i)
    }
  }
  result
}

///|
/// The order is the lcm of cycle lengths, represented without word overflow.
pub fn Permutation::order(self : Permutation) -> @bigint.BigInt {
  let mut result : @bigint.BigInt = 1N
  for cycle in self.cycles() {
    let length = @bigint.BigInt::from_int(cycle.length())
    let mut a = result
    let mut b = length
    while !b.is_zero() {
      let remainder = a % b
      a = b
      b = remainder
    }
    result = result / a * length
  }
  result
}

///|
pub fn Permutation::sign(self : Permutation) -> Int32 {
  let mut transpositions = 0
  for cycle in self.cycles() {
    transpositions += cycle.length() - 1
  }
  if transpositions % 2 == 0 {
    1
  } else {
    -1
  }
}

///|
/// Rotate each cycle by exponent modulo its length. Negative exponents use
/// the inverse. This also handles the minimum signed Int64 without negating it.
pub fn Permutation::power(self : Permutation, exponent : Int64) -> Permutation {
  let images = FixedArray::makei(self.degree(), fn(i) { i })
  for cycle in self.cycles() {
    let length = cycle.length().to_int64()
    let shift = ((exponent % length + length) % length).to_int()
    for i, point in cycle {
      images[point] = cycle[(i + shift) % cycle.length()]
    }
  }
  { images, }
}
