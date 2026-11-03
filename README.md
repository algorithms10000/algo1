# algo1：精確數學演算法庫

MoonBit 實作；首個開發批次把**有限群作用**與**有限體 Galois 對應**接成可編譯、可呼叫的演算法鏈。

## 本批成果

- 6 個 library packages、18 個演算法／演算法族登記項。
- 確定性增量 Schreier–Sims：強生成集、精確群階、成員判定、穩定子、傳遞度。
- 有限 G-集合：軌道傳輸、最小不變等價關係、本原性、乘積／餘積、等變映射延拓、纖維積。
- 模質數線性代數：RREF、秩、核、仿射解、行列式、矩陣快速冪。
- 有限體：Rabin 不可約判準、擴展 Euclid 逆元、Frobenius、跡／範數、最小多項式。
- 顯式 Galois 對應：循環自同構群、固定子體基、子體格與子群格的反序運算。

`algorithm.xlsx` 是定理、演算法、等式及驗證索引；`Dev.txt` 是逐回合流水帳。10,000 條演算法及 5,000,000 行是累計目標，**不是本批宣稱完成的數字**。

## 數學範圍

全局對應命題採有限 Galois 擴張基本定理：對有限 Galois 擴張 E/F，子群 H≤Gal(E/F) 與中間體 F⊆K⊆E 透過固定體及固定子群形成反序一一對應。本庫的可計算實例是 F_(p^n)/F_p；不聲稱任意群作用都有這種對應。

所有定理均附 人類可讀 推導及適用條件；沒有把編譯成功當作機器形式化證明，亦沒有使用隨機成員判定代替精確結果。

## 依賴與動線

```text
permutation ──> group ──> gset
      └───────────────────> gset
linalg ──> field ──> galois
   └────────────────> galois
```

輸入採零基點集合、升冪多項式係數及列主序矩陣。公開構造器驗證雙射、次數、模數、不可約性及映射等變性；回傳可變陣列時防禦複製。

## 已獲得定理列表

軌道—穩定子、Schreier 引理、強生成集完備性、有限群作用區塊判準、等變映射由軌道種子決定、G-集合纖維積、秩—零化度、Rabin 判準、有限體 Frobenius／固定子體、有限體 Galois 對應。完整陳述、證明概要、接口及反例邊界見各資料夾 README。

## 已實作演算法列表

見 `lib/README.md` 的 A01–A18；每項均對應原始碼及 `algorithm.xlsx`，不把構造器或 CI 設定充當演算法條目。

## 編譯與驗證

```sh
curl -fsSL https://cli.moonbitlang.com/install/unix.sh | bash
export PATH="$HOME/.moon/bin:$PATH"
moon check --deny-warn --frozen
moon build --target wasm --deny-warn --frozen
```

本批開發使用官方穩定通道 `moon 0.1.20261102`、`moonc v0.10.19+bee578715`。本地採局部編譯及臨時 `debug_inspect`（現行穩定版中 `inspect` 的正式名稱）；驗證題已在提交前刪除。CI 對 wasm、wasm-gc、js、native 做全局 check/build，不使用 `--release`，不執行測試或基準測試。

## 接入與輸出接口

- `algorithms10000/algo1/lib/permutation`：`Permutation`。
- `algorithms10000/algo1/lib/group`：`PermutationGroup`、`OrbitTree`。
- `algorithms10000/algo1/lib/gset`：`Action`、`EquivariantMap`、`Pullback`。
- `algorithms10000/algo1/lib/linalg`：`Matrix`、`RowEchelon`、`AffineSolution`。
- `algorithms10000/algo1/lib/field`：`FiniteField`、`Element`。
- `algorithms10000/algo1/lib/galois`：`Automorphism`、`SubfieldCorrespondence`。

目前沒有 CLI 或網絡服務；輸出是具型別的庫值及可恢復錯誤。各層 README 給出呼叫動線。
