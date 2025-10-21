####
#
# For first time setup of symbolic links
# Jayden R Palomino for Data in 2025
#
####

#for changing symbolic links of file paths in DANCE Analysis
unlink stage0_bin
ln -s /mnt/hygelac-data/30/dance/jr2514/ang/stage0_bin stage0_bin
unlink stage0_root
ln -s /mnt/hygelac-data/30/dance/jr2514/ang/stage0_root stage0_root
unlink stage1_root
ln -s /mnt/hygelac-data/30/dance/jr2514/ang/stage1_root stage1_root
unlink stage0_simulated
ln -s /mnt/hygelac-data/30/dance/jr2514/ang/stage0_simulated stage0_simulated
unlink stage1.cfg
ln -s cfg_files/stage1_La139.cfg stage1.cfg
unlink stage0.cfg
ln -s cfg_files/stage0_caen2018.cfg stage0.cfg
unlink stage1_simulated.cfg
ln -s cfg_files/stage1_simulated.cfg stage1_simulated.cfg
unlink raw_data
ln -s mnt/hygelac-data/24/dance/caen2018 raw_data


cd Config
unlink DanceMap.txt
ln -s DanceMap_au2019.txt DanceMap.txt
unlink TMatrix.txt
ln -s TMatrix_2019.txt TMatrix.txt

#need to link directories in Alpha Calibrator
cd ../../DANCE_Alpha_Calibrator/DANCE_Alpha_Database
unlink DANCE_Alpha_Database.root
ln -s DANCE_Alpha_Database_113391_113393.root DANCE_Alpha_Database.root

#for calibrations and gates
#Calibrations is its own folder on the same level as DANCE_Analysis
cd ../../Calibrations
unlink calib_ideal.dat
ln -s calib_ideal_2019.dat calib_ideal.dat

cd ../DANCE_Analysis/Gates
unlink Alpha.dat
ln -s Alpha_au2019.dat Alpha.dat
#ln -s Alpha_112775.dat Alpha.dat
#the following is a test gate drawn from run 113391
#ln -s Alpha_109430.dat Alpha.dat 
#ln -s Alpha_107982.dat Alpha.dat


unlink Gamma.dat
ln -s Gamma_au2019.dat Gamma.dat

unlink Retrigger.dat
ln -s Retrigger_au2019.dat Retrigger.dat
#ln -s Retrigger_113210.dat Retrigger.dat

unlink Pileup.dat
ln -s Pileup_au2019.dat Pileup.dat


