param(
    [string]$SourceRoot = "Y:\long\10X_datasets\Xenium\Atera\WTA_Preview_FFPE_Breast_Cancer_outs",
    [string]$RemoteHost = "sscb-a100.scilifelab.se",
    [string]$RemoteDir = "/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$requiredFiles = @(
    "cell_feature_matrix.h5",
    "cells.csv.gz",
    "experiment.xenium",
    "metrics_summary.csv",
    "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
    "xenium_explorer_annotations.geojson"
)

if (-not (Test-Path -LiteralPath $SourceRoot)) {
    throw "SourceRoot does not exist: $SourceRoot"
}

& ssh $RemoteHost "mkdir -p '$RemoteDir'"
if ($LASTEXITCODE -ne 0) {
    throw "Failed to create remote directory $RemoteDir on $RemoteHost"
}

foreach ($fileName in $requiredFiles) {
    $sourcePath = Join-Path $SourceRoot $fileName
    if (-not (Test-Path -LiteralPath $sourcePath)) {
        throw "Required input file is missing: $sourcePath"
    }
    Write-Host "Syncing $fileName -> ${RemoteHost}:$RemoteDir/"
    & scp $sourcePath "${RemoteHost}:$RemoteDir/"
    if ($LASTEXITCODE -ne 0) {
        throw "scp failed for $sourcePath"
    }
}

& ssh $RemoteHost "ls -lh '$RemoteDir'"
