$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$compiler = Join-Path $PSScriptRoot 'tmp/latex-runtime/tectonic.exe'
$source = Join-Path $PSScriptRoot '论文大纲_相对通道动态辨识与刚柔解耦自校准.tex'
if (-not (Test-Path -LiteralPath $compiler)) { throw 'Portable Tectonic compiler is missing.' }
& $compiler -X compile --outdir $PSScriptRoot --outfmt pdf --keep-logs --untrusted $source
if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed: $LASTEXITCODE" }
