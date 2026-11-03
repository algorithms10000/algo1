# linalg

## 定義與局部公理

矩陣元素在 F_p，p 為質數且 2≤p≤2^31−1。此範圍使兩個殘基的乘積小於 2^62，可在 Int64 內精確計算。質數性用確定性試除驗證；這是代數前提，不以模數名稱假定可逆。

## 已獲得定理列表

- **T07.1 RREF 不變量：**列交換、非零倍乘、倍加保持解集；主元正規化並上下消去後得到唯一列簡化階梯形。
- **T07 秩—零化度：**dim ker A=columns−rank A。每個自由變數取 1、其餘自由變數取 0，構造線性無關核基。
- **T07.2 仿射解：**Ax=b 的解集若非空，為一個特解加 ker A；增廣矩陣出現 [0…0|1] 則無解。
- **T07.3 行列式：**列交換變號、倍加不改變行列式；消去後對角乘積配對交換符號即 det。

## 已實作演算法列表

A10 `rref/kernel/solve`；A11 `determinant`；A12 `multiply/apply/power/subtract`。`inverse` 以常數額外記憶體的擴展 Euclid 計算；矩陣冪用二進位分解。

## 複雜度

rref/求解 O(r²c)，核 O(rc²) 量級；行列式與 n×n 乘法 O(n³)；矩陣冪 O(n³ log e) 次模算術。輸入拷貝保證原矩陣不變。質數檢查 O(√p)，每個公開算術入口只在需要時執行。

## 依賴圖、動線與接口

```text
rows,columns,p,entries → Matrix → rref → rank/pivots
                              └→ solve → particular + kernel basis
                              └→ determinant / multiply / power
```

無本地上游；下游 field、galois。輸入列主序 Array[Int64]，輸出防御拷貝。無解返回 None；零元素求逆、非方陣行列式、負矩陣冪等返回 `MatrixError`。

## 驗證邊界

核對奇異／長方形矩陣、四種質數上的秩—零化度、無解系統、行列式、矩陣冪與直接連乘一致。沒有浮點容差或機率性秩判定。
