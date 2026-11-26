# OCFEM_Develop

3次元OCFEMによる非定常Navier–Stokes計算のソース一式です。2026-11-26に、ローカル原本との照合を行って整理しました。

## 原本の扱い

- `src/`と既存の回帰試験のFortranファイルは、原本と**同じファイル名・同じバイト列**です。コメント、計算式、改行コードも変更していません。
- `source_manifest.json`に役割、原本名、原本と登録対象のSHA-256を記録しています。
- `.gitattributes`でFortranと入力データの自動改行変換を無効化しています。
- GitHub Actionsは照合するだけです。ソースを書き換えたり、自動コミットしたりする権限はありません。
- 原本中の歴史的なコメントや内部変数名はそのまま保持しています。実際の役割は下の一覧で確認してください。

## 計算に使用するファイル

| 場所・ファイル名 | 内容 |
| --- | --- |
| `src/solvers/UnSteady_Navire_Stokes_Straight_Duct_OCFEM_Marquardt_PARDISO_Extrap_Newton_lmn_nCSR.f90` | 非定常3次元NSの主プログラム。時間微分作用素を保持し、Marquardt補正をPARDISOで計算 |
| `src/common/marquardt_pardiso_solver.f90` | 通常・拡大系・キャッシュ付き拡大系のPARDISO補正ソルバー |
| `src/common/sparse_triplet_builder.f90` | 疎行列のCOO/CSR組立、転置、補償積和など |
| `src/initial/Steady_Stokes_Straight_Duct_OCFEM_takaki_intel_LU_3D_MM_Bound2_CSR_lmn.f90` | 線形定常Stokesの初期場作成。NSの非線形ソルバーとは別プログラム |
| `src/preprocess/calcollQ_TAKAKI_newton_q_xy.f90` | 1・2方向の内部選点生成。時間選点の入力にも使用 |
| `src/preprocess/calcollQ_TAKAKI_newton_q_xyz.f90` | 3方向の内部選点生成 |
| `src/preprocess/COMV1T_TAKAKI_q.f90` | 1次元微分作用素。時間方向の作用素作成 |
| `src/preprocess/COMV3d_TAKAKI_q_Mxyz.f90` | 3次元の微分作用素作成 |
| `src/preprocess/calweight_TAKAKI_2D_q.f90` | 断面流量積分用の重み作成 |
| `src/preprocess/Straight_Duct_gene_TAKAKI_3D_LMN.f90` | 既存5要素管路の選点座標、接続法線、境界値を生成する前処理。ファイル名は原本を維持 |

旧アップロードの取り違えは`docs/previous_upload_audit.json`に記録しました。旧CG系など、現行非定常主プログラムとは時間項の扱いが異なる版は今回の実行対象に混在させていません。旧履歴は削除せず、通常の変更コミットで置換します。

## 必要環境

- Intel Fortran（ifxまたはifort）
- Intel oneMKL（PARDISO、BLAS/LAPACK）
- Python 3。ビルド補助・照合は標準ライブラリだけで動作します。

現在のビルドスクリプトはデフォルトINTEGER 32bitとMKL LP64を対応させています。HPCでILP64へ変更する場合は、Fortran側の整数幅とMKLのインターフェースを一緒に変更してください。

配列境界チェック付き試験で、原本の`src/common/sparse_triplet_builder.f90`の216行に範囲外参照が検出されています。また、大容量配列を持つ原本のStokesソースは、Windowsの通常ビルドではNTCOFFの2GB制限に当たります。原本を無断修正せず、`docs/VALIDATION.md`に区別して記録しています。現段階を「全試験合格・本計算の収束確認済み」とは扱わないでください。

## Windowsでのビルド

Intel oneAPIとVisual Studio 2022の環境で実行します。

```bat
python tools\verify_sources.py
scripts\build_windows.cmd --mode release
```

コンパイラーやMKLの場所が標準と異なる場合は、`OCFEM_INTEL_COMPILER_ROOT`と`MKLROOT`を指定してください。Pythonの実行ファイルは`OCFEM_PYTHON`で指定できます。ビルド結果は`build/release/bin/`に入ります。デフォルトの`runtime`対象は非定常NS主プログラムと前処理6本です。Stokes初期場作成・試験を含めた対象は`--target all`です。ビルドだけで計算は開始しません。

## Linux/HPCでのビルド

Intel oneAPIコンパイラーとMKLの環境をロードしてから実行します。

```sh
bash scripts/build_linux.sh --mode release
```

このLinux用手順は同じLP64構成を指定しています。Windowsでの確認結果をLinuxでの実行確認と同一視しないでください。

## 小規模確認

```bat
scripts\build_windows.cmd --mode smoke --target all
scripts\test_windows.cmd
```

`smoke`ビルドは、生成する検証用コピーの`KNAA=10000`だけを`KNAA=256`に変え、2×2×2内部選点の小規模確認で大容量配列を確保しない構成です。`src/`内の原本コピーは変更しません。検証用実行ファイルを大きな本計算に使用しないでください。

`examples/smoke_2_2_2/`には、原本名の微分作用素・幾何・境界・保存場と、新しく作成した入力ファイルをまとめています。入力例は全長時間0.01、時間区間0.01の1区間です。公称Re500の旧ファイル名を維持しており、指定値から計算するReは約498です。

`input.in`は非定常NS用、`stokes.in`は初期場作成用です。主プログラムは標準入力からファイル名等を読みます。実行時の作業ディレクトリは入力ファイルのあるディレクトリにしてください。結果はその下の`results/`へ出力します。

検証の実施範囲・成否は`docs/VALIDATION.md`に記録します。今回の整理は計算式の改変ではなく、原本の正しい配置と、依存ファイル・ビルド手順の整備です。
