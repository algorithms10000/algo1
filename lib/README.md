# lib：庫群與依賴視角

## 已獲得定理

T01 軌道—穩定子；T02 Schreier 引理；T03 BSGS 完備性；T04 不變等價關係與本原性；T05 等變映射延拓；T06 纖維積泛性質；T07 秩—零化度；T08 Rabin 判準；T09 有限體 Frobenius；T10 有限體 Galois 對應。各子庫 README 給出假設、證明及對應實作。

## 已實作演算法

| ID | 演算法／族 | 套件 |
|---|---|---|
| A01 | 循環分解、階、符號與整數冪 | permutation |
| A02 | BFS 軌道及因式分解傳輸 | group |
| A03 | Schreier 穩定子生成 | group |
| A04 | 確定性增量 Schreier–Sims、篩分、群階 | group |
| A05 | 自然基穩定子鏈與 k-傳遞度 | group |
| A06 | Atkinson 式不變閉包與本原性 | gset |
| A07 | 對角乘積與不交餘積作用 | gset |
| A08 | 軌道種子等變延拓 | gset |
| A09 | 有限 G-集合纖維積 | gset |
| A10 | 模質數 RREF、核與仿射求解 | linalg |
| A11 | 模質數行列式 | linalg |
| A12 | 二進位矩陣冪與線性映射複合 | linalg |
| A13 | Rabin 不可約判定 | field |
| A14 | 有限體多項式運算與擴展 Euclid 逆元 | field |
| A15 | Frobenius 矩陣、跡與範數 | field |
| A16 | 首個線性相關求最小多項式 | field |
| A17 | 循環自同構合成與階 | galois |
| A18 | 固定子體基及反序格對應 | galois |

## 依賴圖

```text
permutation → group → gset
permutation ─────────→ gset
linalg → field → galois
linalg ──────────────→ galois
```

## 動線及接口

群作用路徑：validated images → Permutation → PermutationGroup / Action → OrbitTree / EquivariantMap / Pullback。

有限體路徑：prime p + monic f → Rabin → FiniteField → Element → Frobenius matrix → kernel(F^d−I) → SubfieldCorrespondence。

輸入輸出細節見各包 README。臨時 inspect 不屬庫接口，也不保留在倉庫。
