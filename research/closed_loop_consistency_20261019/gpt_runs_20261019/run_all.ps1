$ErrorActionPreference = 'Stop'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONIOENCODING='utf-8'
$runRoot = $PSScriptRoot
function Invoke-RunPython([string[]]$RunArgs) {
    & python @RunArgs
    if ($LASTEXITCODE -ne 0) { throw "Python failed: $($RunArgs -join ' ')" }
}
Invoke-RunPython @("$runRoot/experiment_a/run_a.py",'--workers','3')
Invoke-RunPython @("$runRoot/experiment_a/summarize_a.py")
Invoke-RunPython @("$runRoot/experiment_a/paired_analysis.py")
Invoke-RunPython @("$runRoot/comparators/run_sml.py",'0','1')
Invoke-RunPython @("$runRoot/comparators/run_indirect.py",'0','1')
Invoke-RunPython @("$runRoot/comparators/run_indirect_ivgn.py",'0','1')
Invoke-RunPython @("$runRoot/comparators/summarize.py")
Invoke-RunPython @("$runRoot/experiment_b/run_b.py",'--reps','20','--workers','2')
Invoke-RunPython @("$runRoot/experiment_b/summarize_b.py")
foreach ($zPos in @('0','-0.5','0.5','1')) {
    Invoke-RunPython @("$runRoot/experiment_c/run_c.py",$zPos)
    Invoke-RunPython @("$runRoot/experiment_c/run_c_sml.py",$zPos)
}
Invoke-RunPython @("$runRoot/experiment_c/summarize_c.py")
Invoke-RunPython @("$runRoot/plot_results.py")
Invoke-RunPython @("$runRoot/inventory.py",'--verify')
Invoke-RunPython @("$runRoot/build_report.py")
