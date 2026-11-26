#!/usr/bin/env bash
# Clone tvami/QuirkPheno at the pinned commit, apply the charged-quirk decay patch and build.
# Usage: QP_WORK=/path/to/workdir ./apply_and_build.sh    (needs Geant4, Pythia8, HepMC3, ROOT-free build env active)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${QP_WORK:-$PWD/qd_work}"
PIN=d91cb27
mkdir -p "$WORK"
[ -d "$WORK/QuirkPheno" ] || git clone https://github.com/tvami/QuirkPheno.git "$WORK/QuirkPheno"
cd "$WORK/QuirkPheno"
git checkout -q "$PIN"
git apply --check "$HERE/patches/0001-quirk-decay.patch" && git apply "$HERE/patches/0001-quirk-decay.patch"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
cp -r "$HERE/macros" .
echo "built: $WORK/QuirkPheno/build/quirksim/quirksim ; run from $WORK/QuirkPheno with: -m macros/ident.mac"
