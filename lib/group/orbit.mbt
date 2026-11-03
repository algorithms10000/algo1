///|
/// Breadth-first orbit tree with factored transport maps. Parent links store
/// one generator index per point instead of a full permutation per point.
#derive(Debug)
pub struct OrbitTree {
  degree : Int32
  root_point : Int32
  generators : Array[@permutation.Permutation]
  points : Array[Int32]
  parent : FixedArray[Int32]
  edge : FixedArray[Int32]
  position : FixedArray[Int32]
}

///|
pub extend OrbitTree with @debug.Debug::{to_repr}

///|
pub fn OrbitTree::build(
  degree : Int32,
  generators : Array[@permutation.Permutation],
  root : Int32,
) -> OrbitTree raise @permutation.DomainError {
  if degree < 0 {
    raise @permutation.DomainError::NegativeDegree(degree)
  }
  if root < 0 || root >= degree {
    raise @permutation.DomainError::PointOutOfRange(root, degree)
  }
  for generator in generators {
    if generator.degree() != degree {
      raise @permutation.DomainError::DegreeMismatch(generator.degree(), degree)
    }
  }
  let parent = FixedArray::make(degree, -1)
  let edge = FixedArray::make(degree, -1)
  let points = [root]
  let position = FixedArray::make(degree, -1)
  position[root] = 0
  parent[root] = root
  let mut head = 0
  while head < points.length() {
    let point = points[head]
    head += 1
    for index, generator in generators {
      let next = generator.image(point)
      if parent[next] == -1 {
        parent[next] = point
        edge[next] = index
        position[next] = points.length()
        points.push(next)
      }
    }
  }
  { degree, root_point: root, generators, points, parent, edge, position, }
}

///|
pub fn OrbitTree::points(self : OrbitTree) -> Array[Int32] {
  self.points.copy()
}

///|
pub fn OrbitTree::contains(self : OrbitTree, point : Int32) -> Bool {
  point >= 0 && point < self.degree && self.parent[point] != -1
}

///|
/// A group element taking root to point, reconstructed in O(depth * degree).
pub fn OrbitTree::transport(
  self : OrbitTree,
  point : Int32,
) -> @permutation.Permutation raise @permutation.DomainError {
  if !self.contains(point) {
    raise @permutation.DomainError::PointOutOfRange(point, self.degree)
  }
  let path : Array[Int32] = []
  let mut cursor = point
  while cursor != self.root_point {
    path.push(self.edge[cursor])
    cursor = self.parent[cursor]
  }
  let mut result = @permutation.Permutation::identity(self.degree)
  let mut i = path.length() - 1
  while i >= 0 {
    result = self.generators[path[i]].compose(result)
    i -= 1
  }
  result
}

///|
/// Schreier generators t_(s(x))^(-1) ∘ s ∘ t_x for the left action.
pub fn OrbitTree::stabilizer_generators(
  self : OrbitTree,
) -> Array[@permutation.Permutation] raise @permutation.DomainError {
  let transports : Array[@permutation.Permutation] = []
  for point in self.points {
    transports.push(self.transport(point))
  }
  let seen : Map[FixedArray[Int32], Unit] = Map([])
  let result : Array[@permutation.Permutation] = []
  for i, point in self.points {
    for generator in self.generators {
      let next = generator.image(point)
      let target_index = self.position[next]
      let element = transports[target_index]
        .inverse()
        .compose(generator)
        .compose(transports[i])
      let key = element.to_images()
      if !element.is_identity() && !seen.contains(key) {
        seen.set(key, ())
        result.push(element)
      }
    }
  }
  result
}

///|
pub fn orbit_partition(
  degree : Int32,
  generators : Array[@permutation.Permutation],
) -> Array[Array[Int32]] raise @permutation.DomainError {
  let seen = FixedArray::make(degree, false)
  let result : Array[Array[Int32]] = []
  for point in 0..<degree {
    if seen[point] {
      continue
    }
    let tree = OrbitTree::build(degree, generators, point)
    let orbit = tree.points()
    for point_in_orbit in orbit {
      seen[point_in_orbit] = true
    }
    result.push(orbit)
  }
  result
}
