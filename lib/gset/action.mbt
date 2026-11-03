///|
#derive(Debug)
pub(all) suberror ActionError {
  GeneratorCountMismatch(Int32, Int32)
  InvalidMapLength(Int32, Int32)
  InvalidMapImage(Int32, Int32)
  NotEquivariant(Int32, Int32)
  ConflictingSeed(Int32, Int32, Int32)
  MissingSeed(Int32)
  SizeOverflow(Int32, Int32)
}

///|
pub extend ActionError with @debug.Debug::{to_repr}

///|
/// A finite action specified by images of one common abstract generating
/// tuple. Two actions are comparable only with equal generator-list length.
#derive(Eq)
#derive(Debug)
pub struct Action {
  degree : Int32
  generators : Array[@permutation.Permutation]
}

///|
pub extend Action with Eq::{equal, not_equal}

///|
pub extend Action with @debug.Debug::{to_repr}

///|
pub fn Action::new(
  degree : Int32,
  generators : Array[@permutation.Permutation],
) -> Action raise @permutation.DomainError {
  if degree < 0 {
    raise @permutation.DomainError::NegativeDegree(degree)
  }
  for generator in generators {
    if generator.degree() != degree {
      raise @permutation.DomainError::DegreeMismatch(generator.degree(), degree)
    }
  }
  { degree, generators: generators.copy(), }
}

///|
pub fn Action::degree(self : Action) -> Int32 {
  self.degree
}

///|
pub fn Action::generators(self : Action) -> Array[@permutation.Permutation] {
  self.generators.copy()
}

///|
fn checked_product(a : Int32, b : Int32) -> Int32 raise ActionError {
  let limit : Int64 = 2147483647L
  if a.to_int64() * b.to_int64() > limit {
    raise ActionError::SizeOverflow(a, b)
  }
  a * b
}

///|
/// Diagonal action on X × Y, indexed by x*|Y| + y.
pub fn Action::product(self : Action, other : Action) -> Action raise {
  if self.generators.length() != other.generators.length() {
    raise ActionError::GeneratorCountMismatch(
      self.generators.length(),
      other.generators.length(),
    )
  }
  let degree = checked_product(self.degree, other.degree)
  let generators : Array[@permutation.Permutation] = []
  for index, left in self.generators {
    let right = other.generators[index]
    let images = Array::makei(degree, fn(pair) {
      let x = pair / other.degree
      let y = pair % other.degree
      let next_x = try! left.image(x)
      let next_y = try! right.image(y)
      next_x * other.degree + next_y
    })
    generators.push(try! @permutation.Permutation::from_images(images))
  }
  { degree, generators, }
}

///|
/// Disjoint union action, with the second summand shifted by |X|.
pub fn Action::coproduct(self : Action, other : Action) -> Action raise {
  if self.generators.length() != other.generators.length() {
    raise ActionError::GeneratorCountMismatch(
      self.generators.length(),
      other.generators.length(),
    )
  }
  let limit : Int64 = 2147483647L
  if self.degree.to_int64() + other.degree.to_int64() > limit {
    raise ActionError::SizeOverflow(self.degree, other.degree)
  }
  let degree = self.degree + other.degree
  let generators : Array[@permutation.Permutation] = []
  for index, left in self.generators {
    let right = other.generators[index]
    let images = Array::makei(degree, fn(point) {
      if point < self.degree {
        try! left.image(point)
      } else {
        let next_y = try! right.image(point - self.degree)
        self.degree + next_y
      }
    })
    generators.push(try! @permutation.Permutation::from_images(images))
  }
  { degree, generators, }
}

///|
#derive(Debug)
pub struct EquivariantMap {
  source : Action
  target : Action
  images : FixedArray[Int32]
}

///|
pub extend EquivariantMap with @debug.Debug::{to_repr}

///|
/// Generator equivariance implies equivariance for every generated word.
pub fn EquivariantMap::new(
  source : Action,
  target : Action,
  images : Array[Int32],
) -> EquivariantMap raise {
  if source.generators.length() != target.generators.length() {
    raise ActionError::GeneratorCountMismatch(
      source.generators.length(),
      target.generators.length(),
    )
  }
  if images.length() != source.degree {
    raise ActionError::InvalidMapLength(images.length(), source.degree)
  }
  for image in images {
    if image < 0 || image >= target.degree {
      raise ActionError::InvalidMapImage(image, target.degree)
    }
  }
  for index, generator in source.generators {
    let target_generator = target.generators[index]
    for point, image in images {
      let left = images[try! generator.image(point)]
      let right = try! target_generator.image(image)
      if left != right {
        raise ActionError::NotEquivariant(point, index)
      }
    }
  }
  { source, target, images: FixedArray::from_array(images[:]), }
}

///|
pub fn EquivariantMap::image(
  self : EquivariantMap,
  point : Int32,
) -> Int32 raise {
  if point < 0 || point >= self.source.degree {
    raise ActionError::InvalidMapImage(point, self.source.degree)
  }
  self.images[point]
}

///|
pub fn EquivariantMap::to_images(self : EquivariantMap) -> FixedArray[Int32] {
  self.images.copy()
}

///|
/// Unique equivariant extension of compatible seeds, requiring at least one
/// seed in every source orbit. Breadth-first propagation checks consistency.
pub fn extend_map(
  source : Action,
  target : Action,
  seeds : Array[(Int32, Int32)],
) -> EquivariantMap raise {
  if source.generators.length() != target.generators.length() {
    raise ActionError::GeneratorCountMismatch(
      source.generators.length(),
      target.generators.length(),
    )
  }
  let images = FixedArray::make(source.degree, -1)
  for (point, image) in seeds {
    if point < 0 || point >= source.degree {
      raise ActionError::InvalidMapImage(point, source.degree)
    }
    if image < 0 || image >= target.degree {
      raise ActionError::InvalidMapImage(image, target.degree)
    }
    if images[point] != -1 && images[point] != image {
      raise ActionError::ConflictingSeed(point, images[point], image)
    }
    images[point] = image
  }
  let queue : Array[Int32] = []
  for point in 0..<source.degree {
    if images[point] != -1 {
      queue.push(point)
    }
  }
  let mut head = 0
  while head < queue.length() {
    let point = queue[head]
    head += 1
    for index, generator in source.generators {
      let next = try! generator.image(point)
      let next_image = try! target.generators[index].image(images[point])
      if images[next] == -1 {
        images[next] = next_image
        queue.push(next)
      } else if images[next] != next_image {
        raise ActionError::NotEquivariant(point, index)
      }
    }
  }
  for point, image in images {
    if image == -1 {
      raise ActionError::MissingSeed(point)
    }
  }
  EquivariantMap::new(source, target, Array::from_iter(images.iter()))
}
