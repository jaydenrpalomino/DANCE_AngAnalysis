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

  # ensure parent of link exists (should be current dir or subdir)
  ensure_dir "$(dirname "$linkname")"

  # remove any existing file/symlink (don’t recurse)
  if [ -L "$linkname" ] || [ -e "$linkname" ]; then
    echo "🔗 Removing existing: $linkname"
    unlink "$linkname" || rm -f "$linkname"
  fi

  # if target missing, create placeholder (so the link is set up regardless)
  if [ ! -e "$target" ]; then
    echo "⚠️  Target missing: $target — creating placeholder."
    if [[ "$target" == */ ]]; then
      mkdir -p "$target"
    else
      ensure_dir "$(dirname "$target")"
      : > "$target"
    fi
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

# NEW: hard requirement on DANCE_Alpha_Calibrator being present
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
require_alpha_calibrator   # <-- added: stop early if missing

# Ensure current-dir Gates exists (your specific ask)
ensure_dir "./Gates"

# Ensure common external mount points exist (placeholders if missing)
ensure_dir /mnt/hygelac-data/30/dance/jr2514/ang
ensure_dir /mnt/hygelac-data/24/dance/caen2018

# ---- Top-level links (current dir) ----
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_bin        stage0_bin
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_root       stage0_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage1_root       stage1_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_simulated  stage0_simulated
safe_link /mnt/hygelac-data/24/dance/caen2018                     raw_data

safe_link cfg_files/stage1_La139.cfg      stage1.cfg
safe_link cfg_files/stage0_caen2018.cfg   stage0.cfg

# ---- Config/ ----
pushd Config >/dev/null
safe_link DanceMap_au2019.txt  DanceMap.txt
safe_link TMatrix_2019.txt     TMatrix.txt
popd >/dev/null

# ---- ../DANCE_Alpha_Calibrator/DANCE_Alpha_Database ----
pushd ../DANCE_Alpha_Calibrator/DANCE_Alpha_Database >/dev/null
safe_link DANCE_Alpha_Database_113391_113393.root  DANCE_Alpha_Database.root
popd >/dev/null

# ---- ../Calibrations ----
pushd ../Calibrations >/dev/null
safe_link ../DANCE_Analysis/calib_ideal_2019.dat  calib_ideal.dat
popd >/dev/null

# ---- ./Gates (ensure local Gates only) ----
pushd ./Gates >/dev/null
safe_link Alpha_au2019.dat     Alpha.dat
safe_link Gamma_au2019.dat     Gamma.dat
safe_link Retrigger_au2019.dat Retrigger.dat
safe_link Pileup_au2019.dat    Pileup.dat
popd >/dev/null

echo "🎉 Done. All links set up."
