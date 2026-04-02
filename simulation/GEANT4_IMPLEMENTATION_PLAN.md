# Geant4統合実装計画（FAIR向け）

## 目的
- FAIR内部にGeant4ベースのシミュレーション実行系を実装し、`EventStore`に`vector<AHCALSimHit>`と`SimData`を直接投入できるようにする。
- **まずはSimHit出力を最優先**とし、digitization（`DigitizationAlg`）へそのまま接続できる状態を作る。
- 幾何は当面 **HCAL option3** を固定採用し、将来切替可能な設計にする。
- 参照実装 `AHCAL-simulation`（branch: `2026`）と物理的・入出力的に同等化できるよう、比較用の受け入れ基準を明確化する。

---

## 現状整理（FAIR側）

1. FAIRにはすでに`AHCALSimHit`/`SimData`のEDMがあり、`RootWriterAlg`で出力可能。
2. `DigitizationAlg`は`EventStore`の`vector<AHCALSimHit>`を入力に取るため、Geant4統合後の最短パスは「Geant4→SimHit→digitization」。
3. 既存の`SimHitReader`は外部ROOT（AHCAL-simulation形式）を読む設計で、これを内部生成に置き換えるか、並行運用できる。
4. top-level CMakeではGeant4/Pythia8依存がまだ未定義。

---

## 全体アーキテクチャ方針

### A. `IReader`ベースの新入力モジュールとして実装
- 新規: `IO/reader/Geant4SimReader.{hpp,cpp}`
- 役割:
  - 初回`initialize`相当でGeant4 RunManager・Detector・PhysicsList・PrimaryGenerator・Action群を構築
  - `next(...)`呼び出しごとに1イベント実行し、`EventStore`へ
    - `vector<AHCALSimHit>`（例: key=`SimHits`）
    - `SimData`（例: key=`SimData`）
    を投入
- 利点:
  - 既存の`MultiInOut`/`MultiInOneOut`の「reader→algs→writer」流儀を維持
  - 既存の`DigitizationAlg`や`RootWriterAlg`を無改造で再利用可能

### B. Geant4実装本体は`simulation/geant4/`へ分離
推奨ディレクトリ:
- `simulation/geant4/core/` : RunManager, Config bridge, event buffer
- `simulation/geant4/detector/` : option3 geometry construction
- `simulation/geant4/physics/` : physics list factory
- `simulation/geant4/generator/` : primary generator（初期はAHCAL-simulation互換）
- `simulation/geant4/actions/` : run/event/stepping/sensitive detector actions
- `simulation/geant4/decayer/` : Pythia8 decayer adapter

**ポイント**: FAIRの`IO/reader`からは「1イベント実行してSimHitを返す」薄いAPIだけ見せる。

---

## 実装フェーズ（優先順）

## Phase 0: インターフェース固定（最短1週間目標）

### 0-1. コンフィグ仕様を先に固定
`reader.type: Geant4SimReader`を追加し、最低限以下を持つ:

```yaml
reader:
  type: Geant4SimReader
  cfg:
    out_simhits_key: SimHits
    out_simdata_key: SimData
    detector:
      model: hcal_option3
    generator:
      mode: ahcal_compatible
      particle: "mu-"
      energy_GeV: 100.0
      position_mm: [0.0, 0.0, -2000.0]
      direction: [0.0, 0.0, 1.0]
    physics:
      list: FTFP_BERT
    decayer:
      enable_pythia8: true
      tune: "default"
    random:
      seed: 12345
```

### 0-2. `EventStore`出力契約
- 毎イベントで必ず
  - `vector<AHCALSimHit>`
  - `SimData`
  を`put/set`する。
- `SimData`に入らない項目は暫定値（既存`SimHitReader`の流儀）を採用し、将来段階で埋める。

---

## Phase 1: Geant4最小実行（SimHitのみ）

### 1-1. Detector Construction（option3固定）
- まずはAHCAL-simulation `2026`と同じ寸法・材質・layer/chip/channel対応を移植。
- `cellID = layer*100000 + asic*10000 + channel`規約に一致させる。

### 1-2. Sensitive Detector / Hit集約
- stepレベルで同一cellへ加算し、イベント終端で1cell=1`AHCALSimHit`に集約。
- 埋める項目:
  - `cellID`
  - `Edep`（visible energyの定義を参照実装に合わせる）
  - `Nmip = Edep / AHCALGeometry::MIPEnergy`
  - `HitTime`（max deposit time等、参照実装と一致）
  - `TimeOfArrival`

### 1-3. FAIR接続
- `Geant4SimReader::next`で1イベント実行→上記ベクタを返す。
- 既存`RootWriterAlg`で`vector<AHCALSimHit>`をROOT出力できることを確認。

**このPhase完了条件**
- `config/SimTest.yaml`相当の設定で、外部ROOT入力無しに`vector<AHCALSimHit>`が生成される。

---

## Phase 2: AHCAL-simulation同等化

### 2-1. 入射条件同等
- 初期粒子注入（粒子種・エネルギー・位置・角度・time）をAHCAL-simulationと同一実装へ。
- 将来拡張を見据え、generator modeをstrategy化（`ahcal_compatible`, `particle_gun`, `file_input`など）。

### 2-2. Truth/secondaryの取り扱い
- `SimData`の
  - `injected_*`
  - neutrino interaction系
  - `Secondary*`
  を段階的に埋める。
- branch互換の観点で、既存`SimHitReader`が読むROOT相当情報との差分表を作る。

### 2-3. 物理同等性検証
- 同じseed/入射条件で、少なくとも以下を比較:
  - `Nhits/event`
  - `ΣEdep`
  - layer別Edep分布
  - time分布
- 一致基準（例）: 平均値差 < 2–5%、主要分布のKS検定で許容範囲。

---

## Phase 3: PhysicsListのコンフィグ化

### 3-1. PhysicsList factory
- 文字列→G4VModularPhysicsListを生成するFactoryを実装。
- 例: `FTFP_BERT`, `QGSP_BERT`, `QGSP_BIC`, `Shielding`。

### 3-2. Config適用
- `physics.list`をYAMLで選択可能に。
- 不正値時はエラー終了（候補一覧をログ出力）。

### 3-3. 将来拡張
- カット値/ステップ制御/EMオプション等を`physics`セクションで拡張可能に。

---

## Phase 4: Pythia8 decayer導入

### 4-1. Decayer adapter層
- Geant4側から呼べる`Pythia8DecayerService`を独立実装。
- ON/OFFを`decayer.enable_pythia8`で切替。

### 4-2. Link/Build
- CMakeにPythia8検出（`find_package` or custom finder）を追加。
- 未導入環境では`enable_pythia8: false`のみ許容、`true`なら明示エラー。

### 4-3. バリデーション
- 崩壊生成物のPDG/運動量保存のsanity check。
- 参照実装と主要チャネルの分岐・スペクトル比較。

---

## 受け入れテスト計画

1. **Smoke test**: 10イベントで`SimHits`が非空、クラッシュ無し。
2. **Schema test**: `vector<AHCALSimHit>` + `SimData`が`RootWriterAlg`で書ける。
3. **Compatibility test**: AHCAL-simulation同条件比較（統計量差分）。
4. **Pipeline test**: `Geant4SimReader -> DigitizationAlg -> RootWriterAlg`の一気通し。

---

## 実装時の注意点

- Geant4オブジェクト寿命（RunManager/Detector/Actions）をreader寿命に束縛し、イベントごと再初期化は避ける。
- MT対応は後段。初期は単一スレッドで再現性を優先。
- seed管理は`runNumber + event_counter`等で追跡可能に。
- `EventStore` key名はconfig可変にして既存ジョブと衝突しないようにする。

---

## 具体的な着手順（短期）

1. `simulation/geant4/`骨格 + CMakeオプション（`FAIR_WITH_GEANT4`, `FAIR_WITH_PYTHIA8`）
2. `Geant4SimReader`雛形追加（まず固定1粒子・固定physics）
3. option3 detector + sensitive detectorで`AHCALSimHit`生成
4. `config/SimTest_geant4.yaml`追加（SimHit出力のみ）
5. `DigitizationAlg`接続確認
6. AHCAL-simulation `2026`との差分詰め（generator/physics/decayer）

---

## Deliverable定義（第一段階）

- FAIR単体でGeant4イベント生成し、`vector<AHCALSimHit>`をROOT出力できる。
- 出力`SimHit`が既存`DigitizationAlg`入力としてそのまま動作する。
- configで少なくとも `physics.list` と `decayer.enable_pythia8` を設定可能。

---

## AHCAL-simulation（2026）実ファイル確認メモ

以下は、`AHCAL-simulation` branch `2026` の実装を実際に確認した上で、FAIRへ移植対象として優先度を付けた対応表。

### 1) ビルド・依存関係（Geant4/ROOT/Pythia8）
- `CMakeLists.txt`
  - Geant4（UI/Vis付きオプション）を`find_package(Geant4 ...)`で取得。
  - ROOTを`find_package(ROOT REQUIRED)`でリンク。
  - Pythia8は`$PYTHIA8/bin/pythia8-config` + `$PYTHIA8/lib/libpythia8.so`前提でリンク。
- FAIR側反映:
  - `FAIR_WITH_GEANT4` / `FAIR_WITH_PYTHIA8`の2段階オプション化。
  - Pythia8未導入時は`decayer.enable_pythia8: false`のみ許容するガードを実装。

### 2) Detector（幾何・材質）
- `src/DetectorConstruction.cc`
  - HCAL: 18x18セル、40層、`SensitiveX=40.0 mm`, `House_X=40.3 mm`などの主要寸法が定義。
  - 材質/Birks設定（PlasticSciHCAL等）を内部定義。
- FAIR側反映:
  - Phase 1でoption3（当面固定）を移植し、将来`detector.model`で差し替え可能にする。
  - `AHCALGeometry`との整合（cell ID規約・座標系）を最優先確認項目に置く。

### 3) Primary generator（初期粒子入射）
- `src/PrimaryGenerator.cc`
  - `gFaser` TTreeを読み、`pdgc/px/py/pz/E/status/firstMother`から一次粒子を生成。
  - neutrino interactionラベル（`ftagNulabel`）やsecondary情報をイベントごとに抽出。
- FAIR側反映:
  - Phase 2で`generator.mode: ahcal_compatible`としてまず同等実装。
  - その後`particle_gun`等を追加するstrategy設計で拡張。

### 4) Event集約（SimHit相当の元データ）
- `src/EventAction.cc`, `src/RunAction.cc`
  - `vecHcalCellID`, `vecHcalVisibleEdepCell`, `vecHcalHitTimeCell`, `vecHcalToaCell`をROOT treeへ出力。
  - `ParticleEnergy`, `ftagNulabel`, `Interaction_*`, `Secondary*`, truthベクタ群を同時出力。
- FAIR側反映:
  - `AHCALSimHit`（`cellID/Edep/Nmip/HitTime/TimeOfArrival`）へ直接マッピング。
  - `SimData`へ`injected_*`, interaction, secondary系を段階投入。

### 5) Pythia8 decayer
- `src/Pythia8Decayer.cc`, `src/Pythia8DecayerPhysics.cc`
  - 外部decayerとしてPythia8をGeant4 decay processへ接続。
  - Tau/B/D系など対象PDGに対し、既存decay tableの置換ロジックあり。
- FAIR側反映:
  - Phase 4でadapter/service層として移植。
  - configでON/OFF制御し、OFF時はGeant4標準decayへフォールバック。
