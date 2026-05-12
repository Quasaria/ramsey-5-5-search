# ============================================================================
# run_warmstarts.ps1 - Lance v5 en parallele depuis les top-K circulants.
# ============================================================================
#
# Usage :
#   .\run_warmstarts.ps1                       # defauts : top 8, 4h chacun
#   .\run_warmstarts.ps1 -NumStarts 21 -Parallel 8 -Hours 4
#
# Prerequis :
#   1. Avoir lance une fois : python ramsey_55_circulant.py enumerate --top_k 50
#      Ca cree le dossier circulants_top/ avec rank_01..rank_50 + _index.json.
#   2. Le script s'attend a trouver ramsey_55_search.py dans le dossier courant.
#
# Sortie :
#   warmstart_runs/
#     rank_01/
#       ramsey55_best.json
#       ramsey55_best.matrix.txt
#       log.txt
#     rank_02/
#       ...
#
# Chaque run est isole dans son sous-dossier (le script copie le .py dedans
# et lance python depuis la, ainsi les fichiers de sauvegarde ne s'ecrasent
# pas mutuellement). A la fin, on affiche un resume des K5 atteints.
# ============================================================================

param(
    [int]$NumStarts = 8,
    [int]$Parallel = 0,
    [double]$Hours = 4.0,
    [string]$CirculantDir = "circulants_top",
    [string]$OutputDir = "warmstart_runs",
    [string]$ScriptName = "ramsey_55_search.py",
    [double]$TStart = 8.0,
    [double]$WeightK4 = 10.0,
    [double]$VertexMoveProb = 0.10,
    [double]$TwoFlipProb = 0.40,
    [int]$TwoFlipLookahead = 6,
    [int]$RestartAfter = 100000,
    [int]$PerturbSize = 30
)

$ErrorActionPreference = "Stop"

# ----------------------------------------------------------------------------
# Verifications
# ----------------------------------------------------------------------------
$indexPath = Join-Path $CirculantDir "_index.json"
if (-not (Test-Path $indexPath)) {
    Write-Host "ERREUR : index introuvable : $indexPath" -ForegroundColor Red
    Write-Host "Lance d'abord : python ramsey_55_circulant.py enumerate --top_k 50"
    exit 1
}
if (-not (Test-Path $ScriptName)) {
    Write-Host "ERREUR : $ScriptName introuvable dans le dossier courant" -ForegroundColor Red
    exit 1
}

$index = Get-Content $indexPath -Raw | ConvertFrom-Json
$entries = $index.entries | Sort-Object rank | Select-Object -First $NumStarts
if ($entries.Count -lt $NumStarts) {
    Write-Host "Note : $($entries.Count) circulants disponibles (demande: $NumStarts)" -ForegroundColor Yellow
    $NumStarts = $entries.Count
}

if ($Parallel -le 0) { $Parallel = $NumStarts }
$timeSeconds = [int]($Hours * 3600)

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

Write-Host ""
Write-Host "===== WARMSTART PARALLEL =====" -ForegroundColor Cyan
Write-Host "  Circulants : top $NumStarts depuis $CirculantDir"
Write-Host "  Parallele  : $Parallel runs simultanes max"
Write-Host "  Duree/run  : $Hours h ($timeSeconds s)"
Write-Host "  Total CPU  : $($NumStarts * $Hours) h-CPU"
Write-Host "  Sortie     : $OutputDir/"
Write-Host ""

# ----------------------------------------------------------------------------
# Lancement
# ----------------------------------------------------------------------------
$jobs = @()
$startedRanks = @{}

foreach ($entry in $entries) {
    $rank = $entry.rank
    $runDir = Join-Path $OutputDir ("rank_{0:D2}" -f $rank)
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null

    Copy-Item $ScriptName -Destination $runDir -Force
    $srcCirculant = Join-Path $CirculantDir $entry.filename
    $dstCirculant = Join-Path $runDir "warmstart_input.json"
    Copy-Item $srcCirculant -Destination $dstCirculant -Force

    while ((Get-Job -State Running).Count -ge $Parallel) {
        Start-Sleep -Seconds 5
    }

    Write-Host "-> Lancement rank $rank (K5_init=$($entry.k5_total), |S|=$($entry.S_size))" -ForegroundColor Green

    $job = Start-Job -Name ("rank_{0:D2}" -f $rank) -ScriptBlock {
        param($workdir, $time, $tStart, $w4, $vmp, $tfp, $tfl, $ra, $ps)
        Set-Location $workdir
        $logPath = "log.txt"
        & python .\ramsey_55_search.py `
            --start_from warmstart_input.json `
            --time $time `
            --t_start $tStart `
            --weight_k4 $w4 `
            --vertex_move_prob $vmp `
            --two_flip_prob $tfp `
            --two_flip_lookahead $tfl `
            --restart_after $ra `
            --perturb_size $ps `
            *>&1 | Out-File -FilePath $logPath -Encoding utf8
    } -ArgumentList `
        (Resolve-Path $runDir).Path, $timeSeconds, $TStart, $WeightK4, `
        $VertexMoveProb, $TwoFlipProb, $TwoFlipLookahead, $RestartAfter, $PerturbSize

    $jobs += $job
    $startedRanks[$job.Id] = $rank
}

Write-Host ""
Write-Host "Tous les jobs sont lances ($($jobs.Count) total)." -ForegroundColor Cyan
Write-Host "Pour suivre : Get-Job  |  Receive-Job <id>  |  Get-Content $OutputDir\rank_XX\log.txt -Tail 5"
Write-Host ""

# ----------------------------------------------------------------------------
# Attente avec pulse
# ----------------------------------------------------------------------------
Write-Host "Surveillance en cours (Ctrl+C pour arreter le suivi - les jobs continuent)..."
$pulseEvery = 60
while ((Get-Job -State Running).Count -gt 0) {
    Start-Sleep -Seconds $pulseEvery
    $running = Get-Job -State Running
    $done = Get-Job -State Completed
    $failed = Get-Job -State Failed
    $now = Get-Date -Format "HH:mm:ss"
    Write-Host "[$now] running=$($running.Count) done=$($done.Count) failed=$($failed.Count)"
}

Write-Host ""
Write-Host "===== TOUS LES JOBS TERMINES =====" -ForegroundColor Cyan

# ----------------------------------------------------------------------------
# Resume : extraire le best_K5 de chaque run
# ----------------------------------------------------------------------------
Write-Host ""
Write-Host "Resume des resultats :"
Write-Host "  rank   K5_init   K5_final  status"
Write-Host "  ----   -------   --------  ------"

$results = @()
foreach ($entry in $entries) {
    $rank = $entry.rank
    $runDir = Join-Path $OutputDir ("rank_{0:D2}" -f $rank)
    $bestPath = Join-Path $runDir "ramsey55_best.json"
    $logPath = Join-Path $runDir "log.txt"

    $k5Final = "?"
    $status = "no output"
    if (Test-Path $bestPath) {
        try {
            $bestData = Get-Content $bestPath -Raw | ConvertFrom-Json
            if ($bestData.info.k5 -ne $null) {
                $k5Final = $bestData.info.k5
                $status = "ok"
            }
        } catch {
            $status = "parse error"
        }
    }

    $results += [PSCustomObject]@{
        Rank = $rank
        K5Init = $entry.k5_total
        K5Final = $k5Final
        Status = $status
        RunDir = $runDir
    }

    $line = ("  {0,4}   {1,7}   {2,8}  {3}" -f $rank, $entry.k5_total, $k5Final, $status)
    if ($k5Final -is [int] -and $k5Final -lt $entry.k5_total) {
        Write-Host $line -ForegroundColor Green
    } else {
        Write-Host $line
    }
}

$valid = $results | Where-Object { $_.K5Final -is [int] }
if ($valid.Count -gt 0) {
    $best = $valid | Sort-Object K5Final | Select-Object -First 1
    Write-Host ""
    Write-Host ("MEILLEUR GLOBAL : rank {0} avec K5={1}" -f $best.Rank, $best.K5Final) -ForegroundColor Cyan
    $bestJsonPath = Join-Path $best.RunDir "ramsey55_best.json"
    Write-Host ("  -> " + $bestJsonPath)
    Write-Host ""
    Write-Host "Pour verifier ce resultat :"
    Write-Host ("  python ramsey_55_circulant.py verify " + $bestJsonPath)
}

Get-Job | Remove-Job -Force
