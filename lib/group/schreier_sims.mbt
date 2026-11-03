///|
/// Deterministic incremental Schreier–Sims with the natural fixed base.
/// No randomized verification, group enumeration, or probabilistic membership.
#derive(Debug)
struct Level {
  generators : Array[@permutation.Permutation]
  generator_keys : Map[FixedArray[Int32], Unit]
  representatives : Array[@permutation.Permutation?]
  points : Array[Int32]
}

///|
pub extend Level with @debug.Debug::{to_repr}

///|
#derive(Debug)
pub struct PermutationGroup {
  degree : Int32
  levels : Array[Level]
}

///|
pub extend PermutationGroup with @debug.Debug::{to_repr}

///|
fn make_level(
  degree : Int32,
  base_point : Int32,
) -> Level raise @permutation.DomainError {
  let representatives : Array[@permutation.Permutation?] = Array::make(
    degree,
    None,
  )
  representatives[base_point] = Some(@permutation.Permutation::identity(degree))
  {
    generators: [],
    generator_keys: Map([]),
    representatives,
    points: [base_point],
  }
}

///|
fn sift(
  group : PermutationGroup,
  element : @permutation.Permutation,
  start : Int32,
) -> (@permutation.Permutation, Int32)? raise @permutation.DomainError {
  let mut residue = element
  for level_index in start..<group.levels.length() {
    if residue.is_identity() {
      return None
    }
    let level = group.levels[level_index]
    let image = residue.image(level_index)
    match level.representatives[image] {
      None => return Some((residue, level_index))
      Some(representative) =>
        residue = representative.inverse().compose(residue)
    }
  }
  if residue.is_identity() {
    None
  } else {
    Some((residue, group.levels.length()))
  }
}

///|
fn close_chain(
  group : PermutationGroup,
  initial : Array[@permutation.Permutation],
) -> Unit raise @permutation.DomainError {
  let elements = initial.copy()
  let pairs : Array[(Int32, @permutation.Permutation, @permutation.Permutation)] = []
  let mut element_head = 0
  let mut pair_head = 0
  while element_head < elements.length() || pair_head < pairs.length() {
    if element_head < elements.length() {
      let element = elements[element_head]
      element_head += 1
      let (residue, failure) = match sift(group, element, 0) {
        Some(value) => value
        None => continue
      }
      // A new strong generator may enlarge every prefix orbit, even when
      // it fixes that level's base point. Close each affected orbit now.
      for index in 0..<=failure {
        let level = group.levels[index]
        let key = residue.to_images()
        if level.generator_keys.contains(key) {
          continue
        }
        level.generator_keys.set(key, ())
        level.generators.push(residue)
        for point in level.points {
          let representative = match level.representatives[point] {
            Some(value) => value
            None => raise @permutation.DomainError::NotBijection
          }
          pairs.push((index, representative, residue))
        }
        let mut point_head = 0
        while point_head < level.points.length() {
          let point = level.points[point_head]
          point_head += 1
          let prefix = match level.representatives[point] {
            Some(value) => value
            None => raise @permutation.DomainError::NotBijection
          }
          let mut generator_index = 0
          while generator_index < level.generators.length() {
            let generator = level.generators[generator_index]
            generator_index += 1
            let next = generator.image(point)
            if level.representatives[next] is None {
              let transport = generator.compose(prefix)
              level.representatives[next] = Some(transport)
              level.points.push(next)
              let mut strong_index = 0
              while strong_index < level.generators.length() {
                pairs.push((index, transport, level.generators[strong_index]))
                strong_index += 1
              }
            }
          }
        }
      }
    } else {
      let (index, transport, generator) = pairs[pair_head]
      pair_head += 1
      let image = generator.image(transport.image(index))
      let target = match group.levels[index].representatives[image] {
        Some(value) => value
        None => raise @permutation.DomainError::NotBijection
      }
      let schreier = target.inverse().compose(generator).compose(transport)
      if !schreier.is_identity() {
        elements.push(schreier)
      }
    }
  }
}

///|
pub fn PermutationGroup::new(
  degree : Int32,
  generators : Array[@permutation.Permutation],
) -> PermutationGroup raise @permutation.DomainError {
  let _ = @permutation.Permutation::identity(degree)
  for generator in generators {
    if generator.degree() != degree {
      raise @permutation.DomainError::DegreeMismatch(generator.degree(), degree)
    }
  }
  let levels : Array[Level] = []
  for point in 0..<degree {
    levels.push(make_level(degree, point))
  }
  let group = { degree, levels, }
  close_chain(group, generators)
  group
}

///|
pub fn PermutationGroup::degree(self : PermutationGroup) -> Int32 {
  self.degree
}

///|
pub fn PermutationGroup::base(self : PermutationGroup) -> Array[Int32] {
  Array::makei(self.degree, fn(i) { i })
}

///|
pub fn PermutationGroup::contains(
  self : PermutationGroup,
  element : @permutation.Permutation,
) -> Bool raise @permutation.DomainError {
  if element.degree() != self.degree {
    raise @permutation.DomainError::DegreeMismatch(
      element.degree(),
      self.degree,
    )
  }
  sift(self, element, 0) is None
}

///|
pub fn PermutationGroup::order(self : PermutationGroup) -> @bigint.BigInt {
  let mut result : @bigint.BigInt = 1N
  for level in self.levels {
    result = result * @bigint.BigInt::from_int(level.points.length())
  }
  result
}

///|
pub fn PermutationGroup::orbit_sizes(self : PermutationGroup) -> Array[Int32] {
  self.levels.map(fn(level) { level.points.length() })
}

///|
/// Union of the chain's level generators, a strong generating set.
pub fn PermutationGroup::strong_generators(
  self : PermutationGroup,
) -> Array[@permutation.Permutation] {
  let result : Array[@permutation.Permutation] = []
  let seen : Map[FixedArray[Int32], Unit] = Map([])
  for level in self.levels {
    for generator in level.generators {
      let key = generator.to_images()
      if !seen.contains(key) {
        seen.set(key, ())
        result.push(generator)
      }
    }
  }
  result
}

///|
/// Generators of the pointwise stabilizer of the first count base points.
pub fn PermutationGroup::prefix_stabilizer_generators(
  self : PermutationGroup,
  count : Int32,
) -> Array[@permutation.Permutation] raise @permutation.DomainError {
  if count < 0 || count > self.degree {
    raise @permutation.DomainError::PointOutOfRange(count, self.degree)
  }
  let result : Array[@permutation.Permutation] = []
  let seen : Map[FixedArray[Int32], Unit] = Map([])
  for index in count..<self.degree {
    for generator in self.levels[index].generators {
      let key = generator.to_images()
      if !seen.contains(key) {
        seen.set(key, ())
        result.push(generator)
      }
    }
  }
  result
}

///|
/// Largest k for which the action is k-transitive. Natural-base orbit sizes
/// must be n, n-1, ..., n-k+1. The answer is between 0 and degree.
pub fn PermutationGroup::transitivity_degree(self : PermutationGroup) -> Int32 {
  let mut result = 0
  for index, level in self.levels {
    if level.points.length() != self.degree - index {
      break
    }
    result += 1
  }
  result
}

///|
pub fn PermutationGroup::is_transitive(self : PermutationGroup) -> Bool {
  self.degree == 0 ||
  (self.levels.length() > 0 && self.levels[0].points.length() == self.degree)
}
