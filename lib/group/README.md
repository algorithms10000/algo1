# group

## 定義與局部公理

G=⟨S⟩≤Sym(Ω)，自然基 B=(0,…,n−1)，G_i 固定所有 j<i。左作用記作 s(x)。每個軌道點儲存一個傳輸元素 t_x，使 t_x(i)=x。

## 已獲得定理列表

- **T01 軌道—穩定子：**|G_i|=|i^G_i|·|G_{i+1}|，故 |G|=∏_i |i^G_i|。映射 g↦(g(i), 對應陪集) 給出雙射。
- **T02 Schreier 引理：**若 S 生成 G_i，則 {t_(s(x))⁻¹∘s∘t_x : x∈i^G_i,s∈S} 生成 G_(i+1)。將任意穩定子元素拆成傳輸與生成元字，再伸縮消去中間軌道代表。
- **T03 BSGS 完備性：**增量加入篩分餘式、閉合所有受影響軌道及所有新傳輸×強生成元對；終止時所有 Schreier 生成元可篩為恒等，歸納得到各層穩定子生成集。每個失敗層至少增加一個軌道點，保證終止。
- **T03.1 成員判定：**完備鏈上逐層消去像；到恒等即成員，遇到缺少的軌道代表即非成員。
- **T03.2 傳遞度：**作用為 k-傳遞 iff 前 k 個自然基點的自然穩定子在其餘點上傳遞，即軌道大小依次 n,n−1,…,n−k+1。

## 已實作演算法列表

A02 `OrbitTree::build/transport`；A03 `stabilizer_generators`；A04 `PermutationGroup::new/contains/order/strong_generators`；A05 `prefix_stabilizer_generators/transitivity_degree`。採確定性增量閉合，不枚舉整個群作成員判定，不使用隨機驗證。

## 複雜度與實作界限

BFS 軌道 O(n|S|) 次像查詢、O(n) .parent/edge 儲存，傳輸重建 O(depth·n)。本批 BSGS 採保守實作：至多 n² 個失敗餘式、每個加入至多 n 個前綴層；以稠密置換及顯式工作隊列計，最壞時間上界 O(n⁶+m·n²)、儲存 O(n⁴)（不計 BigInt 位元）。這不是聲稱達到 nearly-linear 或文獻最佳 Schreier–Sims 界。

## 依賴圖、動線與接口

```text
permutation.Permutation[] → close_chain(sift + orbit + Schreier pairs)
                         → PermutationGroup → order/contains/stabilizer
                         → OrbitTree → transport
```

上游 permutation；下游 gset。群階輸出 BigInt，避免 n! 溢位。 degree 0 的群階為 1、傳遞度為 0。所有生成元必須同次數。

## 驗證邊界

小次數與獨立生成元閉包比較群階及成員；72 組確定性生成 family 核對各前綴軌道—穩定子；另核對 A5 階 60／3-傳遞、D8 階 16／1-傳遞、S20 階 20!／20-傳遞。臨時題已刪除。
