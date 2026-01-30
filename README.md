# Schaefer / HCR / MSE Workshop（00_workshop_run.Rmd）

## 1. このリポジトリは何か（概要）
このリポジトリは、**Schaefer モデル**を用いた資源量のシミュレーション（Operating Model）、
CPUE の観測モデル、**最尤推定（ML）**によるパラメータ推定、
そして **HCR（Hockey-stick / Type2=2k）**の比較と **MSE（Management Strategy Evaluation）**の流れを
一連で体験できるワークショップ用デモです。

実行入口は **`00_workshop_run.Rmd`** です。ここから `R/` 配下のスクリプトを `source()` して動作します。

## 2. クイックスタート（最短で動かす）
### 推奨：RStudio で Knit
1. `260128_LectureDemo.Rproj` を RStudio で開く
2. `00_workshop_run.Rmd` を開き、Knit を実行

### 代替：R コンソールから実行
```r
rmarkdown::render("00_workshop_run.Rmd")
```

### Windows での注意
- 作業ディレクトリは **リポジトリのルート**にしてください。
- RStudio の Project を開くと自動的にルートが作業ディレクトリになります。
- コンソールから実行する場合は、必要に応じて `setwd("path/to/repo")` を使ってください。

## 3. 必要環境
- R（明示はないため **R 4.x 系を推奨**）
- 必須パッケージ：`ggplot2`, `rmarkdown`, `knitr`
- 任意（高速化）：`Rcpp`（C++ 版の最尤推定を使う場合）

インストール例：
```r
install.packages(c("ggplot2", "rmarkdown", "knitr"))
# 任意（高速化）
install.packages("Rcpp")
```

※ `renv` は導入されていません。必要なら後から追加できます。

## 4. 実行すると何が起きるか（処理フロー）
`00_workshop_run.Rmd` の構成に沿って、概ね次の流れで処理が進みます：

1. 乱数シード設定、必要パッケージの読み込み
2. Schaefer Operating Model による資源量・漁獲のシミュレーション（単一ケース）
3. CPUE 観測モデル（対数正規）による観測値生成
4. 最尤推定（ML）で r, K, q, B1, σ を推定し、真値との比較
5. 感度解析（TTYear, σ, 漁獲の過小報告の影響）
6. HCR の形状確認（Hockey-stick / Type2=2k）
7. MSE ループで HCR の性能比較、指標計算、時系列プロット

## 5. 出力物（生成ファイル）
- **HTML**：`00_workshop_run.html`（Knit/Render の出力。現在は同名ファイルが同梱されています）
- Rmd 内の図表は HTML に埋め込まれます（保存先を明示するコードはありません）

※ 計算負荷は軽量ですが、`use_cpp <- TRUE` の場合は初回に C++ コンパイルが走ります。

## 6. ディレクトリ構成
```
.
├─00_workshop_run.Rmd
├─00_workshop_run.html
├─260128_LectureDemo.Rproj
├─R/
│  ├─fit_schaefer_ml.R          # ML 推定（R / Rcpp）
│  ├─hcr_hockey_stick.R         # Hockey-stick HCR
│  ├─hcr_new2kei_2k.R           # Type2 (2k) HCR
│  ├─metrics_plots.R            # 指標計算・可視化ヘルパ
│  ├─mse_loop.R                 # MSE ループ本体
│  ├─obs_cpue_lognormal.R       # CPUE 観測モデル
│  └─om_schaefer.R              # Schaefer OM（資源動態）
├─output/
│  ├─fig/                       # 参考図（PNG）
│  └─tab/                       # 参考表（CSV）
├─FRA-SA2020-ABCWG01-01.pdf
├─AGENTS.md
├─.gitignore
└─.gitattributes
```

## 7. `R/` スクリプト一覧（最重要）

| File | Role | Key functions / objects |
|---|---|---|
| `R/om_schaefer.R` | Schaefer OM の資源動態シミュレーション | `om_schaefer_sim`, `om_step_schaefer`, `predict_B_schaefer` |
| `R/obs_cpue_lognormal.R` | CPUE 観測モデル（対数正規） | `obs_cpue_lognormal`, `check_cpue_lognormal_mean` |
| `R/fit_schaefer_ml.R` | Schaefer モデルの最尤推定（R / Rcpp） | `fit_schaefer_ml`, `fit_schaefer_ml_cpp` |
| `R/hcr_hockey_stick.R` | Hockey-stick 型 HCR と形状計算 | `hcr_hockey_stick`, `hcr_hockey_shape` |
| `R/hcr_new2kei_2k.R` | Type2 (2k) HCR の ABC 計算 | `hcr_new2kei_2k`, `hcr_new2kei_2k_from_D` |
| `R/metrics_plots.R` | MSE 指標・プロットの共通関数 | `calc_mse_metrics`, `plot_ts_case`, `plot_true_vs_hat`, `plot_hockey_shape`, `plot_2k_shape` |
| `R/mse_loop.R` | MSE ループ（OM→観測→推定→HCR→更新） | `mse_loop` |

**簡易依存関係（読み順）**
- `00_workshop_run.Rmd` → `R/om_schaefer.R`, `R/obs_cpue_lognormal.R`, `R/fit_schaefer_ml.R`, `R/hcr_hockey_stick.R`, `R/hcr_new2kei_2k.R`, `R/metrics_plots.R`, `R/mse_loop.R`
- `R/mse_loop.R` は `om_step_schaefer`, `obs_cpue_lognormal`, `fit_schaefer_ml(_cpp)`, `hcr_*` を内部で利用します。

## 8. よくあるエラーと対処
1. **パッケージが無い**（`ggplot2` / `rmarkdown` / `knitr` / `Rcpp`）
   - `install.packages(...)` で導入してください（Rcpp は任意）。
2. **`source("R/...")` が失敗する**
   - 作業ディレクトリが repo ルートでない可能性があります。`getwd()` を確認し、`setwd()` で修正してください。
3. **Rcpp コンパイルでエラーが出る**
   - `use_cpp <- FALSE` に変更して R 版を使うか、Rcpp のビルド環境を整えてください。
4. **`CPUE_obs has too few positive values` が出る**
   - 入力データやシミュレーション条件を変更した場合に起こります。`TTYear` や `sigma_CPUE` を見直してください。
5. **`Blim must be > Bban` など HCR パラメータエラー**
   - HCR の設定値を確認し、`Bban_ratio < Blim_ratio` になるよう調整してください。
