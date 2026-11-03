///|
#derive(Debug)
pub(all) suberror FieldError {
  InvalidCharacteristic(Int64)
  EmptyModulus
  NotMonic
  ReducibleModulus
  FieldMismatch
  ZeroInverse
  InvalidSubfieldDegree(Int32, Int32)
  DimensionMismatch(Int32, Int32)
}

///|
pub extend FieldError with @debug.Debug::{to_repr}

///|
/// Dense coefficients in ascending degree, with no trailing zeroes.
priv struct Polynomial {
  coefficients : Array[Int64]
}

///|
fn normalize(coefficients : Array[Int64], p : Int64) -> Polynomial {
  let values = coefficients.map(fn(value) { @linalg.residue(value, p) })
  while values.length() > 0 && values[values.length() - 1] == 0L {
    ignore(values.pop())
  }
  { coefficients: values, }
}

///|
fn poly_subtract(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> Polynomial {
  let size = if left.coefficients.length() > right.coefficients.length() {
    left.coefficients.length()
  } else {
    right.coefficients.length()
  }
  normalize(
    Array::makei(size, fn(i) {
      let a = if i < left.coefficients.length() {
        left.coefficients[i]
      } else {
        0L
      }
      let b = if i < right.coefficients.length() {
        right.coefficients[i]
      } else {
        0L
      }
      a - b
    }),
    p,
  )
}

///|
fn poly_remainder(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> Polynomial raise {
  let remainder = left.coefficients.copy()
  let degree = right.coefficients.length() - 1
  let leading_inverse = @linalg.inverse(right.coefficients[degree], p)
  while remainder.length() > degree {
    let factor = remainder[remainder.length() - 1] * leading_inverse % p
    let shift = remainder.length() - 1 - degree
    for i in 0..<=degree {
      let index = shift + i
      remainder[index] = @linalg.residue(
        remainder[index] - factor * right.coefficients[i],
        p,
      )
    }
    while remainder.length() > 0 && remainder[remainder.length() - 1] == 0L {
      ignore(remainder.pop())
    }
  }
  normalize(remainder, p)
}

///|
fn poly_multiply(
  left : Polynomial,
  right : Polynomial,
  modulus : Polynomial,
  p : Int64,
) -> Polynomial raise {
  if left.coefficients.length() == 0 || right.coefficients.length() == 0 {
    return { coefficients: [], }
  }
  let product = Array::make(
    left.coefficients.length() + right.coefficients.length() - 1,
    0L,
  )
  for i, a in left.coefficients {
    for j, b in right.coefficients {
      product[i + j] = (product[i + j] + a * b) % p
    }
  }
  poly_remainder(normalize(product, p), modulus, p)
}

///|
fn poly_power(
  base : Polynomial,
  exponent : UInt64,
  modulus : Polynomial,
  p : Int64,
) -> Polynomial raise {
  let mut result = normalize([1L], p)
  let mut factor = poly_remainder(base, modulus, p)
  let mut remaining = exponent
  while remaining > 0UL {
    if remaining % 2UL == 1UL {
      result = poly_multiply(result, factor, modulus, p)
    }
    remaining /= 2UL
    if remaining > 0UL {
      factor = poly_multiply(factor, factor, modulus, p)
    }
  }
  result
}

///|
fn poly_gcd(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> Polynomial raise {
  let mut a = left
  let mut b = right
  while b.coefficients.length() > 0 {
    let remainder = poly_remainder(a, b, p)
    a = b
    b = remainder
  }
  if a.coefficients.length() == 0 {
    return a
  }
  let scale = @linalg.inverse(a.coefficients[a.coefficients.length() - 1], p)
  normalize(a.coefficients.map(fn(value) { value * scale % p }), p)
}

///|
/// Extended Euclid in F_p[x], returning (g, s, t) with s*a+t*b=g.
fn poly_extended_gcd(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> (Polynomial, Polynomial, Polynomial) raise {
  let mut old_r = left
  let mut r = right
  let mut old_s = normalize([1L], p)
  let mut s = normalize([], p)
  let mut old_t = normalize([], p)
  let mut t = normalize([1L], p)
  while r.coefficients.length() > 0 {
    let quotient = poly_quotient(old_r, r, p)
    let next_r = poly_subtract(old_r, poly_plain_multiply(quotient, r, p), p)
    old_r = r
    r = next_r
    let next_s = poly_subtract(old_s, poly_plain_multiply(quotient, s, p), p)
    old_s = s
    s = next_s
    let next_t = poly_subtract(old_t, poly_plain_multiply(quotient, t, p), p)
    old_t = t
    t = next_t
  }
  if old_r.coefficients.length() > 0 {
    let scale = @linalg.inverse(
      old_r.coefficients[old_r.coefficients.length() - 1],
      p,
    )
    old_r = normalize(old_r.coefficients.map(fn(x) { x * scale % p }), p)
    old_s = normalize(old_s.coefficients.map(fn(x) { x * scale % p }), p)
    old_t = normalize(old_t.coefficients.map(fn(x) { x * scale % p }), p)
  }
  (old_r, old_s, old_t)
}

///|
fn poly_quotient(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> Polynomial raise {
  let remainder = left.coefficients.copy()
  let degree = right.coefficients.length() - 1
  let quotient = Array::make(
    if remainder.length() > degree {
      remainder.length() - degree
    } else {
      0
    },
    0L,
  )
  let leading_inverse = @linalg.inverse(right.coefficients[degree], p)
  while remainder.length() > degree {
    let shift = remainder.length() - 1 - degree
    let factor = remainder[remainder.length() - 1] * leading_inverse % p
    quotient[shift] = factor
    for i in 0..<=degree {
      let index = shift + i
      remainder[index] = @linalg.residue(
        remainder[index] - factor * right.coefficients[i],
        p,
      )
    }
    while remainder.length() > 0 && remainder[remainder.length() - 1] == 0L {
      ignore(remainder.pop())
    }
  }
  normalize(quotient, p)
}

///|
fn poly_plain_multiply(
  left : Polynomial,
  right : Polynomial,
  p : Int64,
) -> Polynomial {
  if left.coefficients.length() == 0 || right.coefficients.length() == 0 {
    return { coefficients: [], }
  }
  let product = Array::make(
    left.coefficients.length() + right.coefficients.length() - 1,
    0L,
  )
  for i, a in left.coefficients {
    for j, b in right.coefficients {
      product[i + j] = (product[i + j] + a * b) % p
    }
  }
  normalize(product, p)
}

///|
fn distinct_prime_factors(value : Int32) -> Array[Int32] {
  let result : Array[Int32] = []
  let mut remaining = value
  let mut divisor = 2
  while divisor <= remaining / divisor {
    if remaining % divisor == 0 {
      result.push(divisor)
      while remaining % divisor == 0 {
        remaining /= divisor
      }
    }
    divisor += 1
  }
  if remaining > 1 {
    result.push(remaining)
  }
  result
}

///|
/// Rabin's deterministic criterion: x^(p^n)=x and for every prime divisor
/// l of n, gcd(x^(p^(n/l))-x, f)=1. Frobenius powers reuse one matrix.
pub fn is_irreducible(
  characteristic : Int64,
  coefficients : Array[Int64],
) -> Bool raise {
  if !@linalg.is_prime(characteristic) {
    raise FieldError::InvalidCharacteristic(characteristic)
  }
  let polynomial = normalize(coefficients, characteristic)
  let n = polynomial.coefficients.length() - 1
  if n < 1 {
    raise FieldError::EmptyModulus
  }
  if polynomial.coefficients[n] != 1L {
    raise FieldError::NotMonic
  }
  if n == 1 {
    return true
  }
  let x = normalize([0L, 1L], characteristic)
  let frobenius_x = poly_power(
    x,
    characteristic.reinterpret_as_uint64(),
    polynomial,
    characteristic,
  )
  let columns = Array::make(n, normalize([1L], characteristic))
  for index in 1..<n {
    columns[index] = poly_multiply(
      columns[index - 1],
      frobenius_x,
      polynomial,
      characteristic,
    )
  }
  let entries : Array[Int64] = []
  for row in 0..<n {
    for column in columns {
      entries.push(
        if row < column.coefficients.length() {
          column.coefficients[row]
        } else {
          0L
        },
      )
    }
  }
  let matrix = try! @linalg.Matrix::new(n, n, characteristic, entries)
  let x_vector = FixedArray::makei(n, fn(i) { if i == 1 { 1L } else { 0L } })
  for divisor in distinct_prime_factors(n) {
    let power = try! matrix.power((n / divisor).to_int64())
    let image = try! power.apply(x_vector)
    let difference = poly_subtract(
      { coefficients: Array::from_iter(image.iter()), },
      x,
      characteristic,
    )
    if poly_gcd(difference, polynomial, characteristic).coefficients.length() >
      1 {
      return false
    }
  }
  let full = try! matrix.power(n.to_int64())
  try! (full.apply(x_vector) == x_vector)
}
