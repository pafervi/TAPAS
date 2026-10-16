#!/usr/bin/env bash
# Smoke tests for TAPAS.
# Covers: package import, CLI help, both examples in --dry-run with the
# expected unique counts, and one real end-to-end run that verifies the
# written files on disk. Fast (<30 s total), no network, no external data.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$REPO_ROOT/tapas.py"
EX01="$REPO_ROOT/examples/01_general_enumeration"
EX02="$REPO_ROOT/examples/02_refinement"

EXPECTED_01=190
EXPECTED_02=160

fail=0
tmpdir=$(mktemp -d)
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

# ---------------------------------------------------------------------------
# T1: package imports and exposes the programmatic API
# ---------------------------------------------------------------------------
echo "[T1] package import and public API"
out=$(cd "$REPO_ROOT" && python -c "
import tapas
for name in ('run_general','run_refine','canonical_label','symmetry_permutations',
             'write_config','parse_config_info','validate_and_center',
             'shells_for_vacancy'):
    assert hasattr(tapas, name), f'missing: {name}'
assert tapas.__version__
print('ok', tapas.__version__)
" 2>&1)
if [[ "$out" == ok* ]]; then
    echo "  OK  $out"
else
    echo "  FAIL  $out"
    fail=1
fi

# ---------------------------------------------------------------------------
# T2: CLI help renders without error
# ---------------------------------------------------------------------------
echo "[T2] CLI --help"
if python "$SCRIPT" --help >/dev/null 2>&1; then
    echo "  OK"
else
    echo "  FAIL  --help exited non-zero"
    fail=1
fi

# ---------------------------------------------------------------------------
# T3: example 01 (general) dry-run count
# ---------------------------------------------------------------------------
echo "[T3] example 01 general enumeration (dry-run)"
got=$(cd "$EX01" && python "$SCRIPT" POSCAR_CeO2_110_demo \
    --nvac 2 --mode all --vac-layers surf,subsurf,inner --red-shell both \
    --host-coord 3 4 --r-nn 2.6 --r-nnn 4.8 \
    --min-vac-dist 3.5 --max-bonds-lost 1 \
    --max-red-combos 15 --seed 42 --symprec 0.01 \
    --out "$tmpdir/t3" --dry-run 2>&1 \
    | awk '/Unique configurations:/ {print $3}')

if [[ "$got" == "$EXPECTED_01" ]]; then
    echo "  OK  unique=$got (expected $EXPECTED_01)"
else
    echo "  FAIL  unique=$got (expected $EXPECTED_01)"
    fail=1
fi

# ---------------------------------------------------------------------------
# T4: example 02 (refinement + --balanced + --dedup-local) dry-run count
# ---------------------------------------------------------------------------
echo "[T4] example 02 refinement (dry-run, --balanced --dedup-local)"
got=$(cd "$EX02" && python "$SCRIPT" POSCAR_CeO2_110_demo \
    --refine input_cfg/ --balanced --dedup-local --symprec 0.01 \
    --out "$tmpdir/t4" --dry-run 2>&1 \
    | awk '/Unique configurations:/ {print $3}')

if [[ "$got" == "$EXPECTED_02" ]]; then
    echo "  OK  unique=$got (expected $EXPECTED_02)"
else
    echo "  FAIL  unique=$got (expected $EXPECTED_02)"
    fail=1
fi

# ---------------------------------------------------------------------------
# T5: real end-to-end run writes the expected files on disk
# ---------------------------------------------------------------------------
echo "[T5] example 01 real run writes all cfg folders + summary.csv"
out_dir="$tmpdir/t5"
(cd "$EX01" && python "$SCRIPT" POSCAR_CeO2_110_demo \
    --nvac 2 --mode all --vac-layers surf,subsurf,inner --red-shell both \
    --host-coord 3 4 --r-nn 2.6 --r-nnn 4.8 \
    --min-vac-dist 3.5 --max-bonds-lost 1 \
    --max-red-combos 15 --seed 42 --symprec 0.01 \
    --out "$out_dir" > /dev/null 2>&1)

n_cfg=$(ls -d "$out_dir"/cfg_* 2>/dev/null | wc -l)
n_poscar=$(find "$out_dir" -name POSCAR | wc -l)
n_info=$(find "$out_dir" -name config_info.txt | wc -l)
has_summary=0; [[ -f "$out_dir/summary.csv" ]] && has_summary=1

if [[ "$n_cfg" == "$EXPECTED_01" && "$n_poscar" == "$EXPECTED_01" \
   && "$n_info" == "$EXPECTED_01" && "$has_summary" == 1 ]]; then
    echo "  OK  cfg=$n_cfg POSCAR=$n_poscar config_info=$n_info summary.csv=yes"
else
    echo "  FAIL  cfg=$n_cfg POSCAR=$n_poscar config_info=$n_info summary.csv=$has_summary"
    fail=1
fi

# ---------------------------------------------------------------------------
# T6: config_info.txt from the real run is in English and parseable
# ---------------------------------------------------------------------------
echo "[T6] config_info.txt is English and parseable"
cfg_info="$out_dir/cfg_00000/config_info.txt"
if ! grep -q "VACANCIES (indices in the ORIGINAL" "$cfg_info"; then
    echo "  FAIL  missing English VACANCIES header"
    fail=1
elif grep -qE "VACANTES|composicion|generado|POOLS por vacante" "$cfg_info"; then
    echo "  FAIL  Spanish strings leaked into config_info.txt"
    fail=1
else
    # Also verify parse_config_info reads it back correctly
    out=$(cd "$REPO_ROOT" && python -c "
from tapas import parse_config_info
sp, n, vac, pools = parse_config_info('$cfg_info', 'O')
assert n == 2
assert len(vac) == 2
assert all(v in pools for v in vac)
print('ok')
" 2>&1)
    if [[ "$out" == "ok" ]]; then
        echo "  OK  English + parse_config_info round-trip"
    else
        echo "  FAIL  parse_config_info error: $out"
        fail=1
    fi
fi

# ---------------------------------------------------------------------------
echo ""
if [[ "$fail" -eq 0 ]]; then
    echo "[test] All 6 tests passed."
else
    echo "[test] FAILED. See output above."
    exit 1
fi
