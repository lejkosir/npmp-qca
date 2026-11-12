# Reruns every experiment used in the report, from scratch.
#   powershell -ExecutionPolicy Bypass -File reproduce.ps1
# Results land in runs\<timestamp>_all_seed<N>\ (report.md, summary.csv, segments.csv).

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

# Main experiment (seed 1) rebuilds the simulator first; seeds 2 and 3 check
# that the verdicts do not depend on the engine's random cell order.
python sweep.py all --seed 1 --name seed1 --rebuild
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
foreach ($seed in 2, 3) {
    python sweep.py all --seed $seed --name "seed$seed"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
