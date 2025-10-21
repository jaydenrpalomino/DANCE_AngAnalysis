#!/bin/bash
#### 
# For first-time setup of symbolic links
# Jayden R. Palomino for Data in 2025
####

# Define a helper function for linking safely
safe_link() {
    local target=$1
    local linkname=$2

    # If symlink or file exists, unlink/remove it
    if [ -L "$linkname" ] || [ -e "$linkname" ]; then
        echo "Unlinking existing: $linkname"
        unlink "$linkname"
    fi

    # If target directory doesn't exist, make it
    if [[ "$target" == */* ]]; then
        local parent_dir
        parent_dir=$(dirname "$target")
        if [ ! -d "$parent_dir" ]; then
            echo "Creating directory: $parent_dir"
            mkdir -p "$parent_dir"
        fi
    fi

    # Create symlink if target exists
    if [ -e "$target" ]; then
        ln -s "$target" "$linkname"
        echo "Linked $linkname → $target"
    else
        echo "⚠️ Warning: target not found: $target"
    fi
}

# --- Main symbolic link setup ---

safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_bin stage0_bin
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_root stage0_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage1_root stage1_root
safe_link /mnt/hygelac-data/30/dance/jr2514/ang/stage0_simulated stage0_simulated
safe_link cfg_files/stage1_La139.cfg stage1.cfg
safe_link cfg_files/stage0_caen2018.cfg stage0.cfg
safe_link /mnt/hygelac-data/24/dance/caen2018 raw_data

# --- Config directory ---
cd Config || exit
safe_link DanceMap_au2019.txt DanceMap.txt
safe_link TMatrix_2019.txt TMatrix.txt

# --- Alpha Calibrator ---
cd ../../DANCE_Alpha_Calibrator/DANCE_Alpha_Database || exit
safe_link DANCE_Alpha_Database_113391_113393.root DANCE_Alpha_Database.root

# --- Calibrations ---
cd ../../Calibrations || exit
safe_link calib_ideal_2019.dat calib_ideal.dat

# --- Gates ---
cd ../DANCE_Analysis/Gates || exit
safe_link Alpha_au2019.dat Alpha.dat
safe_link Gamma_au2019.dat Gamma.dat
safe_link Retrigger_au2019.dat Retrigger.dat
safe_link Pileup_au2019.dat Pileup.dat

echo "✅ All symbolic links configured."
