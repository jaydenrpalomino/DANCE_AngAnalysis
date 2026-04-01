#!/bin/bash
####
# DANCE Analysis symlink setup (run from DANCE_Analysis/)
# Jayden R. Palomino — 2025
####

set -euo pipefail

# -----------------------
# Helpers
# -----------------------
ensure_dir() {
  local d="$1"
  if [ ! -d "$d" ]; then
    echo "📁 Creating directory: $d"
    mkdir -p "$d"
  fi
}

safe_link() {
  local target="$1"
  local linkname="$2"

  # ensure parent of link exists
  ensure_dir "$(dirname "$linkname")"

  # if target missing, warn and skip without deleting anything
  if [ ! -e "$target" ]; then
    echo "⚠️  Target missing: $target — skipping link: $linkname"
    return 0
  fi

  # refuse to replace a real directory
  if [ -d "$linkname" ] && [ ! -L "$linkname" ]; then
    echo "⚠️  $linkname exists as a real directory — refusing to replace it. Skipping."
    return 0
  fi

  # remove existing file or symlink
  if [ -L "$linkname" ] || [ -e "$linkname" ]; then
    echo "🔗 Removing existing: $linkname"
    rm -f "$linkname"
  fi

  ln -s "$target" "$linkname"
  echo "✅ Linked $linkname → $target"
}

require_here() {
  # minimal sanity check that we’re in DANCE_Analysis/
  if [ ! -d "./Config" ]; then
    echo "❌ This script must be run from the DANCE_Analysis/ directory (no ./Config found)."
    exit 1
  fi
}

# hard requirement on DANCE_Alpha_Calibrator being present
require_alpha_calibrator() {
  if [ ! -d "../DANCE_Alpha_Calibrator" ] || [ ! -d "../DANCE_Alpha_Calibrator/DANCE_Alpha_Database" ]; then
    echo "❌ Missing dependency: ../DANCE_Alpha_Calibrator (and/or its DANCE_Alpha_Database)."
    echo "   This repository is required before linking."
    echo "   Get it here: https://github.com/lansce-nuclear-physics/DANCE_Alpha_Calibrator"
    echo
    echo "   From the parent directory of DANCE_Analysis, run:"
    echo "     cd .."
    echo "     git clone https://github.com/lansce-nuclear-physics/DANCE_Alpha_Calibrator.git"
    echo
    echo "   Then re-run this script from DANCE_Analysis/."
    exit 1
  fi
}

# -----------------------
# Start
# -----------------------
echo "=== DANCE Analysis symlink setup (cwd: $PWD) ==="
require_here
require_alpha_calibrator

# Ensure current-dir Gates exists
ensure_dir "./Gates"

# Ensure common external mount points exist
ensure_dir /data2/lansce/jr2514
ensure_dir /data2/lansce/jr2514/raw_data

# ---- Top-level links (current dir) ----
safe_link /data2/lansce/jr2514/stage0_bin/test        stage0_bin
safe_link /data2/lansce/jr2514/stage0_root/test       stage0_root
safe_link /data2/lansce/jr2514/stage1_root/test       stage1_root
safe_link /data2/lansce/jr2514/stage0_simulated  stage0_simulated
safe_link /data2/lansce/jr2514/raw_data          raw_data

# raw data
safe_link /data/lansce/7862-A/raw_data /data2/lansce/jr2514/raw_data

# cfg files
safe_link cfg_files/stage1_Co60_AngularAnalysis_wf.cfg             stage1.cfg
safe_link cfg_files/stage0_caen2018_wf.cfg       stage0.cfg

# ---- Config/ ----
pushd Config >/dev/null
safe_link DanceMap_23Nov15.txt  DanceMap.txt
safe_link TMatrix_57Fe.txt      TMatrix.txt
popd >/dev/null

# ---- ../DANCE_Alpha_Calibrator/DANCE_Alpha_Database ----
ensure_dir "../DANCE_Alpha_Calibrator/DANCE_Alpha_Database"
pushd ../DANCE_Alpha_Calibrator/DANCE_Alpha_Database >/dev/null
safe_link DANCE_Alpha_Database_107980_107982.root  DANCE_Alpha_Database_2018.root
safe_link DANCE_Alpha_Database_2018.root           DANCE_Alpha_Database.root
popd >/dev/null

# ---- ../Calibrations ----
DANCE_ANALYSIS_DIR="$PWD"
ensure_dir "../Calibrations"
pushd ../Calibrations >/dev/null
safe_link "$DANCE_ANALYSIS_DIR/calib_ideal_2018.dat" calib_ideal.dat
popd >/dev/null

# ---- ./Gates ----
pushd ./Gates >/dev/null
safe_link Alpha_107982.dat     Alpha.dat
safe_link Gamma_109430.dat     Gamma.dat
safe_link Retrigger_071918.dat Retrigger.dat
safe_link Pileup_2018.dat      Pileup.dat
popd >/dev/null

echo "🎉 Done. All links set up."