# permutation

## 定義與局部公理

Ω_n={0,…,n−1}，Permutation 是 Ω_n 到自身的雙射。構造器檢查範圍及單射性；有限等勢集合上的單射即雙射。內部陣列私有，公開 `to_images` 回傳拷貝。合成 `p.compose(q)` 表示 p∘q，右邊先作用。

## 已獲得定理列表

- **T01.1 循環分解：**每個置換唯一分解為不交循環（忽略循環次序及旋轉）。逐點沿有限雙射迭代必回到起點，visited 標記保證線性掃描。
- **T01.2 階與符號：**ord(g)=lcm(循環長度)，sgn(g)=(−1)^(Σ(ℓ−1))。不交循環獨立作用；長度 ℓ 的循環可分解為 ℓ−1 個對換。
- **T01.3 負冪邊界：**g^k 在每個循環上按 k mod ℓ 位移；因此 Int64 最小負值亦無需取絕對值。

## 已實作演算法列表

A01：`cycles`、`order`、`sign`、`power`，以及支援其閉合運算的 `compose`、`inverse`、`from_cycles`。`from_cycles` 的列表由左至右複合，允許不同循環共享點，拒絕同一循環內重複點。

## 複雜度

循環分解 O(n)；m 個循環的階為 O(n) 次機器操作加 O(m) 次 BigInt gcd/乘除；冪 O(n)。 BigInt 位元成本另計。構造器 O(n)，合成／逆元 O(n)。

## 依賴圖、動線與接口

```text
Array[Int] → from_images → Permutation → compose/inverse/power
                                  └→ cycles → order/sign
```

無上游本地包；下游為 group、gset。輸入點必須在 [0,degree)，錯誤為 `DomainError`。空次數 0 有唯一恒等置換。輸出包括 `FixedArray[Int]` 拷貝、`Array[Array[Int]]` 循環及 BigInt 階。

## 驗證邊界

已核對重疊循環次序、逆元、−20…20 的冪、Int64 兩端、非法非雙射輸入。數學證明為紙本推導，非 `moon prove` 證書。
