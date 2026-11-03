///|
/// The pullback X ×_Z Y in finite G-sets, with its two projections.
#derive(Debug)
pub struct Pullback {
  left_target : Action
  right_target : Action
  object : Action
  left_projection : FixedArray[Int32]
  right_projection : FixedArray[Int32]
}

///|
pub extend Pullback with @debug.Debug::{to_repr}

///|
pub fn Pullback::new(
  left : EquivariantMap,
  right : EquivariantMap,
) -> Pullback raise {
  if left.target != right.target {
    raise ActionError::GeneratorCountMismatch(
      left.target.degree,
      right.target.degree,
    )
  }
  if left.source.generators.length() != right.source.generators.length() {
    raise ActionError::GeneratorCountMismatch(
      left.source.generators.length(),
      right.source.generators.length(),
    )
  }
  let pairs : Array[(Int32, Int32)] = []
  let index = FixedArray::make(
    try! checked_product(left.source.degree, right.source.degree),
    -1,
  )
  for x in 0..<left.source.degree {
    for y in 0..<right.source.degree {
      if left.images[x] == right.images[y] {
        index[x * right.source.degree + y] = pairs.length()
        pairs.push((x, y))
      }
    }
  }
  let generators : Array[@permutation.Permutation] = []
  for generator_index in 0..<left.source.generators.length() {
    let images = Array::makei(pairs.length(), fn(point) {
      let (x, y) = pairs[point]
      let next_x = try! left.source.generators[generator_index].image(x)
      let next_y = try! right.source.generators[generator_index].image(y)
      index[next_x * right.source.degree + next_y]
    })
    generators.push(try! @permutation.Permutation::from_images(images))
  }
  {
    left_target: left.source,
    right_target: right.source,
    object: { degree: pairs.length(), generators, },
    left_projection: FixedArray::makei(pairs.length(), fn(i) { pairs[i].0 }),
    right_projection: FixedArray::makei(pairs.length(), fn(i) { pairs[i].1 }),
  }
}

///|
pub fn Pullback::action(self : Pullback) -> Action {
  self.object
}

///|
pub fn Pullback::projections(
  self : Pullback,
) -> (EquivariantMap, EquivariantMap) {
  (
    try! EquivariantMap::new(
      self.object,
      self.left_target,
      Array::from_iter(self.left_projection.iter()),
    ),
    try! EquivariantMap::new(
      self.object,
      self.right_target,
      Array::from_iter(self.right_projection.iter()),
    ),
  )
}
