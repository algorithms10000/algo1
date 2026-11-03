# field

## 定義與局部公理

E=F_p[x]/(f)，f 單首、不可約、deg f=n≥1。元素以次數<n 的升冪係數表示。模數經 Rabin 判準驗證後才構造 FiniteField；不同模數的同階體不自動視為同一個體。

## 已獲得定理列表

- **T08 Rabin 判準：**f 不可約 iff x^(p^n)≡x mod f，且對 n 的每個質因數 ℓ，gcd(x^(p^(n/ℓ))−x,f)=1。不可約因子的次數整除 n；真因子的次數必除某個 n/ℓ。
- **T08.1 逆元：**a≠0 時 gcd(a,f)=1，擴展 Euclid 給 sa+tf=1，所以 s 是 a 的模 f 逆元。
- **T09 Frobenius：**φ(a)=a^p 是 F_p-線性自同構，階為 n；φ^k 由快取矩陣的 k 次冪計算。
- **T09.1 跡／範數：**Tr(a)=Σ_(i=0)^(n−1)φ^i(a)，N(a)=∏φ^i(a)，皆在 F_p；由循環置換求和／求積即得。
- **T09.2 最小多項式：**1,a,…,a^n 線性相關；最小 d 的首個關係给出 d 次單首消去多項式。

## 已實作演算法列表

A13 `is_irreducible`；A14 `Element::add/subtract/multiply/inverse/power`；A15 `frobenius_power/trace/norm`；A16 `minimal_polynomial`。多項式冪用 UInt64 指數，負元素冪先求逆，安全處理 Int64 最小值。

## 複雜度

稠密教科書乘法 O(n²)；Rabin 的 Frobenius 矩陣構造 O(n³+n² log p)，各次矩陣冪 O(n³ log n)，gcd 另計多項式 Euclid 成本。 inverse 為 O(n²) 次模算術的擴展 Euclid；跡／範數 O(n³) 量級；最小多項式逐次求解，保守 O(n⁵) 量級。

## 依賴圖、動線與接口

```text
p,f → Rabin → FiniteField → Element → inverse/power/trace/norm
                       └→ Frobenius Matrix → power
Element powers → linalg.solve(first dependence) → minimal_polynomial
```

上游 linalg；下游 galois。輸入 p 不超過 2^31−1，f 必須單首不可約。輸出為元素、固定長度係數或 Frobenius 矩陣。零逆元、可約模數、不同體運算均顯式報錯。

## 驗證邊界

核對 F16 全 16 元素的交換律、逆元消去、a^15、φ^4、跡／範數落於基體；另核對 F9、極端 signed exponent、可約模數拒絕及最小多項式。例子 x⁴+x+1 的係數是 [1,1,0,0,1]。
