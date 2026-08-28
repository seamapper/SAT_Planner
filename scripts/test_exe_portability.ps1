param([string]$ExeName = "SAT_Planner_v2026.37.exe")

$env:PATH = "C:\Windows\System32;C:\Windows"
Remove-Item Env:GDAL_DATA, Env:PROJ_LIB, Env:PROJ_DATA, Env:CONDA_PREFIX -ErrorAction SilentlyContinue

$exe = Join-Path $PSScriptRoot "..\dist\$ExeName" | Resolve-Path
$log = Join-Path $env:USERPROFILE "sat_planner_geospatial_error.log"
if (Test-Path $log) { Remove-Item $log }

Write-Output "Testing $ExeName"
$p = Start-Process -FilePath $exe -PassThru
Start-Sleep -Seconds 8
if (Get-Process -Id $p.Id -ErrorAction SilentlyContinue) {
    Write-Output "STILL_RUNNING"
    Stop-Process -Id $p.Id -Force
} else {
    Write-Output "EXITED"
}
if (Test-Path $log) {
    Write-Output "--- ERROR LOG ---"
    Get-Content $log
} else {
    Write-Output "No geospatial error log"
}
