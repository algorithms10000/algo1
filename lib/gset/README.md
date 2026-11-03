# gset

## 定義與局部公理

兩個 Action 可組合的前提是生成元列表長度相同，第 j 對代表同一抽象生成元在不同集合上的作用。等變映射滿足 f(s_jx)=s_jf(x)。

## 已獲得定理列表

- **T04 最小不變閉包：**由 a~b 出發，每次合併等價類便對所有生成元傳播像，得到最小 G-不變等價關係。有限置換的逆是正冪，正向生成元閉合已足夠。
- **T04.1 本原判準：**次數≥2 的傳遞作用本原， iff 任意兩不同點生成的不變等價關係只有單點類或全集合；等價類即區塊。
- **T05 等變延拓：**一個軌道上的等變映射由一個種子唯一決定；傳播一致 iff 種子被該軌道穩定子固定。每軌道至少一個相容種子時得到唯一全域延拓。
- **T06 纖維積：**X×_Z Y={(x,y):f(x)=g(y)} 在對角作用下封閉，兩個投影等變；任何使 f∘u=g∘v 的 W 與等變映射唯一分解到該配對集合。

## 已實作演算法列表

A06 `invariant_closure/is_primitive`；A07 `Action::product/coproduct`；A08 `extend_map`；A09 `Pullback::new/projections`。乘積索引 x|Y|+y，餘積第二段平移 |X|；纖維積只保留相容配對並重新編號。

## 複雜度

不變閉包 O(n|S|·α(n))；本原性逐點呼叫，O(n²|S|·α(n))。乘積／餘積 O((n+m)|S|) 或 O(nm|S|)。等變延拓 O(n|S|)。纖維積 O(nm+k|S|)，k 為相容配對數；目前索引空間 O(nm)。32 位次數乘法先在 Int64 檢查。

## 依賴圖、動線與接口

```text
permutation → Action → product/coproduct/invariant_closure
                    → EquivariantMap → Pullback → projections
group ─────────────────────────────→ orbit/primitivity checks
```

輸入為 Action、種子對或兩個共享目標的 EquivariantMap；輸出為 Action、完整映射表或 Pullback。錯誤區分生成元數、目標、映射範圍、種子衝突、缺種子及尺寸溢位。

## 驗證邊界

核對 C4 的 {0,2}/{1,3} 區塊、本原 S4、乘積軌道、餘積、等變平移、衝突種子，以及纖維積的交換方形。未把任意函數裝作等變映射。
