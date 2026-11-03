///|
priv struct Partition {
  parent : FixedArray[Int32]
  size : FixedArray[Int32]
}

///|
fn partition_root(partition : Partition, point : Int32) -> Int32 {
  let mut root = point
  while partition.parent[root] != root {
    root = partition.parent[root]
  }
  let mut cursor = point
  while cursor != root {
    let next = partition.parent[cursor]
    partition.parent[cursor] = root
    cursor = next
  }
  root
}

///|
fn partition_blocks(partition : Partition) -> Array[Array[Int32]] {
  let index = FixedArray::make(partition.parent.length(), -1)
  let blocks : Array[Array[Int32]] = []
  for point in 0..<partition.parent.length() {
    let root = partition_root(partition, point)
    if index[root] == -1 {
      index[root] = blocks.length()
      blocks.push([])
    }
    blocks[index[root]].push(point)
  }
  blocks
}

///|
/// Least invariant equivalence relation containing a~b. Each successful
/// union propagates all generator images. Forward closure suffices because
/// each finite permutation's inverse is a positive power of itself.
pub fn invariant_closure(
  action : Action,
  a : Int32,
  b : Int32,
) -> Array[Array[Int32]] raise {
  if a < 0 || a >= action.degree || b < 0 || b >= action.degree {
    raise ActionError::InvalidMapImage(a, action.degree)
  }
  let partition = {
    parent: FixedArray::makei(action.degree, fn(i) { i }),
    size: FixedArray::make(action.degree, 1),
  }
  let queue : Array[(Int32, Int32)] = [(a, b)]
  let mut head = 0
  while head < queue.length() {
    let (left, right) = queue[head]
    head += 1
    let mut left_root = partition_root(partition, left)
    let mut right_root = partition_root(partition, right)
    if left_root == right_root {
      continue
    }
    if partition.size[left_root] < partition.size[right_root] {
      let swap = left_root
      left_root = right_root
      right_root = swap
    }
    partition.parent[right_root] = left_root
    partition.size[left_root] += partition.size[right_root]
    for generator in action.generators {
      queue.push((try! generator.image(left), try! generator.image(right)))
    }
  }
  partition_blocks(partition)
}

///|
/// A transitive action of degree at least two is primitive exactly when
/// every pair of distinct points generates the universal invariant relation.
pub fn is_primitive(action : Action) -> Bool {
  if action.degree < 2 {
    return false
  }
  let orbits = try! @group.orbit_partition(action.degree, action.generators)
  if orbits.length() != 1 {
    return false
  }
  for point in 1..<action.degree {
    let blocks = try! invariant_closure(action, 0, point)
    for block in blocks {
      if block.length() > 1 && block.length() < action.degree {
        return false
      }
    }
  }
  true
}
