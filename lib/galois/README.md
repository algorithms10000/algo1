# galois

## 定義與全局對應的嚴格實例

E/F 有限 Galois 時，H↦E^H、K↦Gal(E/K) 是子群格與中間體格的反序雙射。本包實作 E=F_(p^n)、F=F_p 的情形；Gal(E/F)=⟨φ⟩≅C_n。這不把任意有限群作用冒充域擴張。

## 已獲得定理列表

- **T10.1 循環子群：**每個子群唯一為 H_d=⟨φ^d⟩，d|n，|H_d|=n/d。φ^k 的階為 n/gcd(n,k)。
- **T10 固定子體：**E^(H_d)=ker(F^d−I)=F_(p^d)，維數 d。固定點構成子體；有限體唯一子體階數及 Galois 基本定理給出等號。
- **T10.2 反序格：**E_d∩E_e=E_gcd(d,e)，E_dE_e=E_lcm(d,e)；對應子群運算分別為 join=H_gcd、intersection=H_lcm。
- **T10.3 自同構合成：**φ^a∘φ^b=φ^(a+b mod n)，逆元 φ^(−a mod n)。

## 已實作演算法列表

A17 `Automorphism::new/apply/compose/inverse/order`；A18 `SubfieldCorrespondence::new/contains`、交／合成次數、子群 join/intersection 指數、`subfield_lattice`。子體基由實際 kernel 計算，不按次數填入預設答案。

## 複雜度

列舉因數 O(√n)，每個 d 的固定基需矩陣冪 O(n³ log d) 加消去 O(n³)。完整格共 τ(n) 項；自同構合成 O(1)，作用 O(n³ log exponent)（稠密矩陣冪）。

## 依賴圖、動線與接口

```text
field.FiniteField → Frobenius matrix → power(d) → subtract(I)
                                      └→ linalg.kernel → fixed-field basis
n → divisors → H_d ↔ E_d → intersection/compositum
```

上游 field、linalg；目前無本地依賴者。輸入 d 必須整除 n；輸出基、子群生成元、子群階及格運算指數。非因子次數與不同 ambient field 明確報錯。

## 驗證邊界

F16 的子體次數 [1,2,4] 對應子群階 [4,2,1]；逐元素核對固定點數 2^d、基對乘法封閉。F64 核對 [1,2,3,6]；F9 核對奇特徵 [1,2]。亦核對自同構合成、逆元、階及非法子體次數。
