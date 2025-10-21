#!/bin/bash
####
# Automatic symbolic link setup for DANCE Analysis
# Jayden R. Palomino — 2025
#
# This script safely creates or replaces all symbolic links and
# makes any missing directories or config files needed.
####

set -e  # Exit on critical error
set -u  # Treat unset vars as error

# -----------------------
# Helper functions
# -----------------------

# Create directory if missing
ensure_dir() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        echo "📁 Creating directory: $dir"
        mkdir -p "$dir"
    fi
}

# Safe link creation — replaces existing links or files, ensures directory exists
safe_link() {
    local target=$1
    local linkname=$2

    # Ensure parent directory for the link exists
    ensure_dir "$(dirname "$linkname")"

    # Unlink existing if needed
    if [ -L "$linkname" ] || [ -e "$linkname" ]; then
        echo "🔗 Removing existing: $linkname"
        unlink "$linkname" || rm -f "$linkname"
    fi

    # Ensure target’s parent dir exists if it’s a local relative path
    if [[ "$target" != /* ]]; then
        ensure_dir "$(dirname "$target")"
    fi

    # Warn if target missing (but still make link)
    if [ ! -e "$target" ]; then
        echo "⚠️  Target missing ($target), creating empty placeholder."
        if [[ "$target" == */ ]]; then
            mkdir -p "$target"
        else
            ensure_dir "$(dirname "$target")"
            touch "$target"
        fi
    fi

    ln -s "$target" "$linkname"
    echo "✅ Linked $linkname → $target"
}

# -----------------------
# Core setup
# -----------------------

echo "=== Setting up DANCE Analysis symbolic links ==="

# Main working directory (adjust if needed)
ensure_dir /mnt/hygelac-data/30/dance/jr2514/ang
ensure_dir /mnt/hygelac-data/24/dance/caen2018

safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_bin stage0_bin
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_root stage0_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage1_root stage1_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_simulated stage0_simulated
safe_link /mnt/hygelac-data/24/dance/caen2018 raw_data

safe_link cfg_files/stage1_La139.cfg stage1.cfg
safe_link cfg_files/stage0_caen2018.cfg stage0.cfg

# -----------------------
# Config directory setup
# -----------------------
ensure_dir Config
cd Config
safe_link DanceMap_au2019.txt DanceMap.txt
safe_link TMatrix_2019.txt TMatrix.txt
cd ..

# -----------------------
# Alpha Calibrator setup
# -----------------------
ensure_dir DANCE_Alpha_Calibrator/DANCE_Alpha_Database
cd DANCE_Alpha_Calibrator/DANCE_Alpha_Database
safe_link DANCE_Alpha_Database_113391_113393.root DANCE_Alpha_Database.root
cd ../../

# -----------------------
# Calibrations setup
# -----------------------
ensure_dir Calibrations
cd Calibrations
safe_link calib_ideal_2019.dat calib_ideal.dat
cd ..

# -----------------------
# Gates setup
# -----------------------
ensure_dir DANCE_Analysis/Gates
cd DANCE_Analysis/Gates
safe_link Alpha_au2019.dat Alpha.dat
safe_link Gamma_au2019.dat Gamma.dat
safe_link Retrigger_au2019.dat Retrigger.dat
safe_link Pileup_au2019.dat Pileup.dat
cd ../../

echo "🎉 DANCE symbolic link setup complete!"
