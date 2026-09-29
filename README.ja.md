# neo-smo

[English](README.md) | [简体中文](README.zh-CN.md) | **日本語**

neo-smo は [zsmooth](https://github.com/adworacz/zsmooth) から直接移植した VapourSynth と AviSynth 用の映像フィルタープラグインです。空間・時間方向のノイズ除去、近傍を用いた修復、色ノイズ除去、ブロック DCT フィルタリングを提供します。

元の Zig ソースコードを基に、C++17 と Google Highway で実装・最適化し、ブロック DCT には neo-libdct を使用します。DualSynth2 が計算コアを両ホストに接続します。VapourSynth では `core.neo_smo`、AviSynth では `neo_smo_` 接頭辞の関数を使用します。

## プロジェクトの由来

zsmooth は Austin Dworaczyk Wiltshire（adworacz）が開発しました。neo-smo のアルゴリズムとインターフェースは zsmooth を直接継承しています。zsmooth 自体も、映像フィルターのコミュニティが蓄積してきたアルゴリズムと実装を基にしています。そのドキュメントには、以下の出典や参考実装が挙げられています。

- [AviSynth RgTools](https://github.com/pinterf/RgTools) と [VapourSynth RemoveGrain](https://github.com/vapoursynth/vs-removegrain)：RemoveGrain 系列の従来実装と動作の参考。
- [VapourSynth TemporalSoften](https://github.com/dubhater/vapoursynth-temporalsoften2)、[VapourSynth TemporalMedian](https://github.com/dubhater/vapoursynth-temporalmedian)、[Neo Temporal Median](https://github.com/HomeOfAviSynthPlusEvolution/neo_TMedian)：時間方向の平滑化とメディアンフィルター。
- [VapourSynth FluxSmooth](https://github.com/dubhater/vapoursynth-fluxsmooth)：FluxSmooth 系列。
- [Dogway の AviSynth Scripts](https://github.com/Dogway/Avisynth-Scripts)：`ex_median`、四分位平均、SmartMedian などのアルゴリズムの着想。
- [End-of-Eternity の CCD](https://github.com/End-of-Eternity/vs-ccd) と [vs-jetpack の CCD](https://github.com/Jaded-Encoding-Thaumaturgy/vs-jetpack/blob/e0f47d86930150fd0bf92b0845ccc2b0491f7807/vsdenoise/ccd.py#L95)：zsmooth CCD の主な参考実装。CCD は当初 Sergey Stolyarevsky が VirtualDub 用に作成しました。
- [Cnr2](http://avisynth.nl/index.php/Cnr2)：zsmooth Cnr4 のアルゴリズムの出典の一つ。

zsmooth とこれらの上流プロジェクトの作者・貢献者に感謝します。

## 設計

neo-smo はフィルターの計算とホストのフレーム管理を分離しています。コアは画像プレーンのサンプリング、ソート、重み付け、変換を担当し、ホスト層はパラメーター、フレーム要求、プロパティ、出力の割り当てを管理します。コアは単独でビルドおよびテストできます。DCTFilter は neo-libdct の固定 8×8 変換を使用し、FFTW ランタイムに依存しません。

移植では zsmooth の公開関数名とパラメーター構造を維持しつつ、見つかった問題を修正し、計算経路を最適化しています。上流の不具合の修正、FMA、浮動小数点の丸めは結果に影響する場合があり、過去のすべての zsmooth バージョンとの画素単位の一致は保証しません。呼び出し方と数値規則は[移行ガイド（英語）](docs/api/en/migration.md)を参照してください。

フィルターと DCT はホストの呼び出しスレッド上で処理を行い、ワーカースレッドやスレッドプールを作成しません。ホストは複数フレームを並行して要求できます。SIMD やバッチ DCT は内部のマルチスレッド処理を意味しません。

## 対応する処理

| 分類 | 関数と用途 |
|---|---|
| 空間方向の順序統計とクリーニング | `Median`、`InterQuartileMean`、`SmartMedian`、`VerticalCleaner`：中央値、四分位平均、条件付きメディアン補正、垂直方向のクリーニング。 |
| 近傍を用いたノイズ除去と修復 | `RemoveGrain`、`Repair`：複数モードの空間処理と参照近傍による値の制限。 |
| 時間方向の値の制限と修復 | `Clense`、`ForwardClense`、`BackwardClense`、`TemporalRepair`：隣接フレームを用いて現在の画素値を制約。 |
| 時間・時空間ノイズ除去 | `TemporalMedian`、`TemporalSoften`、`DegrainMedian`、`FluxSmoothT`、`FluxSmoothST`、`TTempSmooth`：時間方向の順序統計、しきい値付き平均、方向選択、重み付き平滑化。 |
| しきい値付き近傍平滑化 | [Deen、MiniDeen](docs/api/en/deen.md)：空間・3 フレームのノイズ除去。固定分母、適応型、距離重み付きの平均。 |
| 色ノイズ除去 | `CCD`、`Cnr4`：色差で誘導する空間サンプリングと時間方向の色度ノイズ除去。 |
| ブロック周波数フィルタリング | `DCTFilter`：重ならない 8×8 ブロックの DCT 係数に重み付け。 |

全部で 21 関数です。通常は、形式とサイズが固定のプレーナー GRAY/YUV/RGB に対応し、8–16 ビット整数、16 ビット浮動小数点（F16）、32 ビット浮動小数点（F32）を扱います。F16 は VapourSynth のみで利用できます。AviSynth+ はアルファを持たないプレーナー形式に対応します。TTempSmooth は F16、CCD は GRAY に対応せず、Cnr4 は整数 YUV のみ、MiniDeen は整数入力のみを受け付けます。

出力は入力の形式、サイズ、フレーム数、フレームレートを維持します。`planes` を持つ関数では、省略時に全プレーンを処理し、`planes=[0]` で最初のプレーンだけを処理します。空のプレーン配列は受け付けません。一部の関数は `mode` やしきい値でプレーン処理を制御するため、すべての関数に一律に `planes` を渡すことはできません。これらのフィルターは動きベクトルの推定や動き補償を行いません。

## ドキュメントと使用方法

API リファレンスでは関数の呼び出し方を、ナレッジベースではサンプリング、ソート、重み付け、変換によって入力から出力を求める過程を解説します。英語版を参照してください。

- [API reference (English)](docs/api/en/README.md)：関数シグネチャ、パラメーター、既定値、使用例。
- [Knowledge base (English)](docs/knowledge/en/README.md)：データ表現、数式、演算順序、境界、精度。
- [zsmooth からの移行（英語）](docs/api/en/migration.md)：名前空間、しきい値の単位、シーンチェンジ、精度の規則。

ビルドしたプラグインを明示的に読み込むか、VapourSynth のプラグイン自動読み込みディレクトリに配置してください。以下は Windows のファイル名を使用しています。Linux では `neo-smo.so` を使用し、他のプラットフォームでも実際のプラグインパスに置き換えてください。VapourSynth のプラグイン識別子は `org.neofilters.neo_smo` です。

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-smo.dll")

clip = core.std.BlankClip(width=640, height=360, format=vs.YUV420P8, length=24)
output = core.neo_smo.Median(clip, radius=[1], planes=[0])
output.set_output()
```

この最小例は合成クリップを使って輝度のメディアン処理を示します。実際の使用時は `clip` を入力映像に置き換えてください。`scalep` を持つ関数は、明示されたしきい値を既定で入力形式の単位として解釈します。省略時の内部既定値は、すでにビット深度に応じてスケーリングされている場合があります。スクリプトを移行する際は各関数のパラメーター説明を確認してください。

同じプラグインファイルに AviSynth C++ インターフェースも含まれます。インターフェースのバージョン 11 に対応したホストで `LoadPlugin` を使用してください。

```avs
LoadPlugin("/path/to/neo-smo.dll")
clip = BlankClip(width=640, height=360, length=24, pixel_type="YV12")
return neo_smo_Median(clip, radius=[1], planes=[0])
```

関数名と引数の順序は API リファレンスに対応し、関数名に `neo_smo_` 接頭辞を付けます。配列パラメーターは `[0, 1]` などのネイティブ配列を受け取り、単一の値も 1 要素の配列として扱います。DCTFilter の `factors` には引き続き正確に 8 要素が必要です。Cnr4 の `sense`、`str`、`pow` はそれぞれ正確に 3 要素が必要で、単一の値では代用できません。真偽値には `true`/`false` を使用します。音声とフィールドパリティは主入力クリップから引き継ぎ、出力フレームのプロパティは対応する入力フレームから取得します。[AviSynth インターフェース（英語）](docs/api/en/README.md#avisynth-calls-and-builds)を参照してください。

パラメーターの省略を明示するには、VapourSynth では `None`、AviSynth では `Undefined()` を使用するか、引数自体を省略します。省略時は既定の動作が選ばれ、ゼロや空配列を渡すこととは異なります。例えば `planes=None` または `planes=Undefined()` は既定のプレーン選択を使用しますが、neo-smo は `planes=[]` を受け付けません。

## SIMD と CPU 選択

Highway は、ビルドに含まれ、実行中の CPU が対応している SIMD ターゲットを自動選択します。公開の `opt`、SIMD OFF、`KernelInfo` インターフェースはありません。

半精度演算を使用する経路では、コンパイル対象と CPU が対応していれば F16 で直接計算し、それ以外では F32 で計算して F16 に書き戻します。アルゴリズム上必要な場合は、より広い中間型も使用します。DCTFilter は常に F32 で変換し、F16 は入出力の格納にのみ使用します。

FMA、演算順序、中間精度は、丸めやしきい値判定に影響することがあります。より広い SIMD や F16 格納が、すべてのフィルターで高速になるとは限りません。[サンプル、しきい値、精度（英語）](docs/knowledge/en/shared/sample-and-precision.md)を参照してください。

## ビルドとテスト

CMake 3.24 以降、Git、C++17 対応コンパイラーが必要です。テストには Python 3 も必要です。CMake は固定バージョンの DualSynth2、Highway、neo-libdct を取得します。両ホストの SDK はローカルで検出するか、自動取得します。

以下の手順は、映像ホストをインストールせずに両ホスト用プラグインとコアテストをビルドします。Windows では clang-cl を優先し、MinGW ビルドでは AVS インターフェースを無効にしてください。Linux では Clang または GCC を使用できます。ネイティブ半精度経路にはコンパイラーと CPU の両方の対応が必要で、Clang 22 以降を推奨します。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

Windows で clang-cl を使用する場合は、Visual Studio の開発者コマンドプロンプトから Ninja ビルドを構成します。

```sh
cmake -S . -B build/clang-cl -G Ninja -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/clang-cl --parallel 4
ctest --test-dir build/clang-cl --output-on-failure
```

| オプション | 用途 |
|---|---|
| `NEO_SMO_BUILD_VAPOURSYNTH=OFF` | VapourSynth インターフェースを無効化。 |
| `NEO_SMO_BUILD_AVISYNTH=OFF` | AviSynth インターフェースを無効化。両ホストを無効化すると、コアとそのテストのみをビルド。 |
| `BUILD_TESTING=OFF` | テストを無効化。 |
| `NEO_SMO_TEST_VAPOURSYNTH=ON` | VapourSynth ホストテストを有効化。既定では無効。指定した Python 環境に VapourSynth、NumPy、同じアーキテクチャのランタイムが必要。 |
| `Python3_EXECUTABLE=/path/to/python` | ホストテスト用の Python を指定。 |
| `NEO_SMO_VS_SDK=/path/to/sdk` | ローカルの VapourSynth SDK を指定。 |
| `NEO_SMO_AVS_SDK=/path/to/sdk` | ローカルの AviSynth SDK を指定。 |
| `NEO_SMO_TEST_AVISYNTH=ON` | AviSynth ホストテストを有効化。既定では無効。Python と同じアーキテクチャのランタイムが必要。 |
| `NEO_SMO_AVISYNTH_RUNTIME=/path/to/avisynth.dll` | AviSynth ホストテスト用のランタイムライブラリを指定。 |
| `NEO_SMO_TEST_CROSS_HOST=ON` | AVS/VS の出力比較を有効化。既定では無効。両ホスト、NumPy、VapourSynth の Python 環境が必要。 |
| `NEO_SMO_TEST_BLACKBOX=ON` | zsmooth とのブラックボックス比較を有効化。既定では無効。VapourSynth の Python 環境、NumPy、参照プラグインが必要。 |
| `NEO_SMO_REFERENCE_PLUGIN=/path/to/zsmooth` | ブラックボックステストで使う参照プラグインファイルを指定。 |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | 固定バージョンのダウンロードの代わりにローカルの DualSynth2 ソースを使用。 |
| `FETCHCONTENT_SOURCE_DIR_NEO_LIBDCT=/path/to/neo-libdct` | 固定バージョンのダウンロードの代わりにローカルの neo-libdct ソースを使用。 |

プラグインのビルドターゲットは `neo_smo`、出力ファイルのベース名は `neo-smo` で、既定では両ホストのインターフェースを含みます。テストはサンプリング境界、整数・浮動小数点演算、半精度経路、DCT、ホストの動作を対象とします。スカラーの参照計算との照合に加え、ビルドに含まれ、実行環境が対応する各 Highway ターゲット間の数値比較も行います。比較は各経路の精度規則に従い、すべての浮動小数点経路でビット単位の一致を要求するものではありません。zsmooth とのブラックボックス比較には、別途参照プラグインが必要です。

CI は Windows x64、Linux x64/ARM64、macOS ARM64、Linux ASan/UBSan を対象とします。Linux runner は Ubuntu 26.04、システム既定の GCC、Clang 22 を使用します。リリースワークフローは Windows、Linux、macOS の x64/ARM64 パッケージを作成します。VapourSynth ホストテストは現在 Windows x64 で実行し、AviSynth ホストテストは上記オプションで個別に有効化します。実際のビルド環境とテスト範囲はパッケージ内に記録され、すべての Linux ディストリビューションへの対応を意味するものではありません。

## 性能

以下は VapourSynth で測定した構成について、**neo-smo のフレーム生成時間 / zsmooth のフレーム生成時間**をフィルター別にまとめたものです。**1 未満は neo-smo の方が高速であることを示します。** 以下の範囲は、計測した 8/16 ビット整数、F16/F32、カラー形式、パラメーターの各設定をまとめたもので、平均値や信頼区間ではありません。各フィルターには、そのフィルターが対応する計測済みの設定のみを含めています。

| フィルター | 所要時間比の範囲 |
|---|---:|
| Median | 0.40–1.21× |
| InterQuartileMean | 0.31–1.23× |
| SmartMedian | 0.70–1.16× |
| VerticalCleaner | 0.50–1.00× |
| RemoveGrain | 0.42–1.14× |
| Repair | 0.48–1.02× |
| Clense | 0.50–1.00× |
| ForwardClense | 0.49–0.98× |
| BackwardClense | 0.49–0.98× |
| TemporalMedian | 0.49–1.02× |
| TemporalSoften | 0.66–1.19× |
| TemporalRepair | 0.27–1.20× |
| DegrainMedian | 0.29–0.68× |
| FluxSmoothT | 0.53–1.30× |
| FluxSmoothST | 0.76–1.18× |
| TTempSmooth | 0.55–0.58× |
| CCD | 0.69–1.37× |
| Cnr4 | 0.68–1.08× |
| DCTFilter | 0.23–0.28× |

これらの数値は個々のフィルターのフレーム生成時間を測定したもので、処理チェーン全体のスループットを示すものではありません。F16 の結果は、対応する上流の F32 経路を比較基準としています。実際の速度は入力、パラメーター、ハードウェア、スレッド数によって変わります。

## 開発と貢献

メンテナーが技術方針、変更のレビュー、リリースを担当します。不具合報告、提案、貢献を歓迎します。数値動作、公開インターフェース、重要な設計変更については、実装前に目的と方針を相談してください。

本プロジェクトでは実装、テスト、レビューに AI を活用します。貢献には問題、方法、検証内容、AI の関与を記載してください。報告にはバージョン、OS、CPU、コンパイラー、ビルド設定、入出力形式、最小再現例を含めてください。数値差分の報告には参照バージョン、パラメーター、要求順序を、性能報告には計測範囲とスレッド設定も記載してください。

## 謝辞とライセンス

neo-smo は以下のライブラリも使用しています。

- [Google Highway](https://github.com/google/highway)：クロスプラットフォームの SIMD を提供します。
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2)：VapourSynth、AviSynth と共有計算コアを接続します。
- [neo-libdct](https://github.com/HomeOfAviSynthPlusEvolution/neo-libdct)：DCTFilter の固定 8×8 ブロック変換に使用します。

テスト、問題報告、改善に協力する開発者とユーザーの皆様に感謝します。

開発に使用する LLM サブスクリプションをご支援いただいた [SB.SB](https://sb.sb) に感謝します。

zsmooth の元コードの MIT 著作権表示とライセンス全文は、[zsmooth ライセンスファイル](LICENSES/zsmooth-MIT.txt)に保持しています。

neo-smo は GNU General Public License バージョン 2 以降（`GPL-2.0-or-later`）で提供します。全文は [LICENSE](LICENSE) を参照してください。第三者コンポーネントはそれぞれの著作権表示とライセンス条項を維持します。
