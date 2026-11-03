# 持續整合

此資料夾只容納 GitHub Actions 設定。從官方穩定通道安裝 MoonBit，對支援的四種後端做全局檢查與編譯。沒有測試、release 或基準測試步驟。

依賴：workflows/compile.yml。輸入：push 至 main 或 dev/**；輸出：各後端編譯結果。
已證定理：不適用。已實作演算法：不適用（CI 設定由 workflows 管轄）。
