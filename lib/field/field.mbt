///|
/// The validated quotient F_p[x]/(modulus), with a cached Frobenius matrix.
#derive(Eq, Debug)
pub struct FiniteField {
  characteristic : Int64
  modulus : FixedArray[Int64]
  frobenius : @linalg.Matrix
}

///|
pub extend FiniteField with Eq::{equal, not_equal}

///|
pub extend FiniteField with @debug.Debug::{to_repr}

///|
pub fn FiniteField::new(
  characteristic : Int64,
  modulus : Array[Int64],
) -> FiniteField raise {
  let irreducible = try! is_irreducible(characteristic, modulus)
  if !irreducible {
    raise FieldError::ReducibleModulus
  }
  let polynomial = normalize(modulus, characteristic)
  let n = polynomial.coefficients.length() - 1
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
  {
    characteristic,
    modulus: FixedArray::from_array(polynomial.coefficients[:]),
    frobenius: try! @linalg.Matrix::new(n, n, characteristic, entries),
  }
}

///|
pub fn FiniteField::characteristic(self : FiniteField) -> Int64 {
  self.characteristic
}

///|
pub fn FiniteField::degree(self : FiniteField) -> Int32 {
  self.modulus.length() - 1
}

///|
pub fn FiniteField::modulus(self : FiniteField) -> FixedArray[Int64] {
  self.modulus.copy()
}

///|
pub fn FiniteField::frobenius_matrix(self : FiniteField) -> @linalg.Matrix {
  self.frobenius
}

///|
#derive(Eq, Debug)
pub struct Element {
  field : FiniteField
  coefficients : FixedArray[Int64]
}

///|
pub extend Element with Eq::{equal, not_equal}

///|
pub extend Element with @debug.Debug::{to_repr}

///|
pub fn FiniteField::element(
  self : FiniteField,
  coefficients : Array[Int64],
) -> Element raise {
  let reduced = poly_remainder(
    normalize(coefficients, self.characteristic),
    { coefficients: Array::from_iter(self.modulus.iter()), },
    self.characteristic,
  )
  {
    field: self,
    coefficients: FixedArray::makei(self.degree(), fn(i) {
      if i < reduced.coefficients.length() {
        reduced.coefficients[i]
      } else {
        0L
      }
    }),
  }
}

///|
pub fn FiniteField::zero(self : FiniteField) -> Element raise {
  self.element([])
}

///|
pub fn FiniteField::one(self : FiniteField) -> Element raise {
  self.element([1L])
}

///|
pub fn FiniteField::generator(self : FiniteField) -> Element raise {
  self.element([0L, 1L])
}

///|
pub fn Element::field(self : Element) -> FiniteField {
  self.field
}

///|
pub fn Element::coefficients(self : Element) -> FixedArray[Int64] {
  self.coefficients.copy()
}

///|
fn require_same(left : Element, right : Element) -> Unit raise {
  if left.field != right.field {
    raise FieldError::FieldMismatch
  }
}

///|
pub fn Element::add(self : Element, right : Element) -> Element raise {
  try! require_same(self, right)
  self.field.element(
    Array::makei(self.coefficients.length(), fn(i) {
      self.coefficients[i] + right.coefficients[i]
    }),
  )
}

///|
pub fn Element::subtract(self : Element, right : Element) -> Element raise {
  try! require_same(self, right)
  self.field.element(
    Array::makei(self.coefficients.length(), fn(i) {
      self.coefficients[i] - right.coefficients[i]
    }),
  )
}

///|
pub fn Element::negate(self : Element) -> Element raise {
  self.field.element(
    Array::from_iter(self.coefficients.iter()).map(fn(value) { -value }),
  )
}

///|
pub fn Element::multiply(self : Element, right : Element) -> Element raise {
  try! require_same(self, right)
  let p = self.field.characteristic
  let modulus = { coefficients: Array::from_iter(self.field.modulus.iter()), }
  let product = poly_multiply(
    { coefficients: Array::from_iter(self.coefficients.iter()), },
    { coefficients: Array::from_iter(right.coefficients.iter()), },
    modulus,
    p,
  )
  self.field.element(product.coefficients)
}

///|
pub fn Element::inverse(self : Element) -> Element raise {
  let p = self.field.characteristic
  let (gcd, coefficient, _) = poly_extended_gcd(
    { coefficients: Array::from_iter(self.coefficients.iter()), },
    { coefficients: Array::from_iter(self.field.modulus.iter()), },
    p,
  )
  if gcd.coefficients.length() == 0 || gcd.coefficients[0] != 1L {
    raise FieldError::ZeroInverse
  }
  self.field.element(coefficient.coefficients)
}

///|
pub fn Element::power(self : Element, exponent : Int64) -> Element raise {
  let p = self.field.characteristic
  let modulus = { coefficients: Array::from_iter(self.field.modulus.iter()), }
  let mut base = { coefficients: Array::from_iter(self.coefficients.iter()), }
  let mut remaining = exponent.reinterpret_as_uint64()
  if exponent < 0L {
    base = {
      coefficients: Array::from_iter((try! self.inverse()).coefficients.iter()),
    }
    // UInt64 can represent abs(Int64::minimum), unlike signed negation.
    remaining = (-(exponent + 1L)).reinterpret_as_uint64() + 1UL
  }
  self.field.element(poly_power(base, remaining, modulus, p).coefficients)
}

///|
pub fn Element::frobenius_power(
  self : Element,
  exponent : Int64,
) -> Element raise {
  let matrix = try! self.field.frobenius.power(exponent)
  self.field.element(
    Array::from_iter(try! matrix.apply(self.coefficients).iter()),
  )
}

///|
pub fn Element::is_zero(self : Element) -> Bool {
  for value in self.coefficients {
    if value != 0L {
      return false
    }
  }
  true
}

///|
pub fn Element::trace(self : Element) -> Element {
  let mut result = try! self.field.zero()
  let mut conjugate = self
  for _ in 0..<self.field.degree() {
    result = try! result.add(conjugate)
    conjugate = try! conjugate.frobenius_power(1L)
  }
  result
}

///|
pub fn Element::norm(self : Element) -> Element {
  let mut result = try! self.field.one()
  let mut conjugate = self
  for _ in 0..<self.field.degree() {
    result = try! result.multiply(conjugate)
    conjugate = try! conjugate.frobenius_power(1L)
  }
  result
}

///|
/// First linear dependence among 1,a,...,a^d gives the monic minimal
/// polynomial. Searching d in increasing order guarantees minimality.
pub fn Element::minimal_polynomial(self : Element) -> FixedArray[Int64] raise {
  let n = self.field.degree()
  let powers = [try! self.field.one()]
  for degree in 1..<=n {
    powers.push(try! powers[degree - 1].multiply(self))
    let entries : Array[Int64] = []
    for row in 0..<n {
      for column in 0..<degree {
        entries.push(powers[column].coefficients[row])
      }
    }
    let matrix = try! @linalg.Matrix::new(
      n,
      degree,
      self.field.characteristic,
      entries,
    )
    let right = powers[degree].coefficients.map(fn(x) { -x })
    let solution = try! matrix.solve(Array::from_iter(right.iter()))
    match solution {
      Some(value) => {
        let coefficients = value.particular()
        return FixedArray::makei(degree + 1, fn(i) {
          if i == degree {
            1L
          } else {
            coefficients[i]
          }
        })
      }
      None => ()
    }
  }
  raise FieldError::DimensionMismatch(n, n + 1)
}
