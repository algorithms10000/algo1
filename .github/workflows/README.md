# 全局編譯工作流程

`compile.yml` 在 main、dev/** 推送及手動觸發時執行。矩陣後端：wasm、wasm-gc、js、native。

動線：checkout → 官方 latest 穩定安裝 → 記錄版本 → moon check --deny-warn --frozen → moon build --deny-warn --frozen。

輸入為倉庫程式碼；輸出為 GitHub Actions 狀態。依賴官方 MoonBit 發行服務及標準庫。此資料夾沒有數學定理或獨立演算法實作。
