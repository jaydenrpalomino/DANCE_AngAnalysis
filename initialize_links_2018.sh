####
#
# For first time setup of symbolic links
# Jayden R Palomino for Data in 2018
#
####

#for changing symbolic links of file paths in DANCE Analysis
unlink stage0_bin
ln -s /data2/lansce/jr2514/stage0_bin stage0_bin
unlink stage0_root
ln -s /data2/lansce/jr2514/stage0_root stage0_root
unlink stage1_bin
ln -s /data2/lansce/jr2514/stage1_bin stage1_bin
unlink stage1_root
ln -s /data2/lansce/jr2514/stage1_root stage1_root
unlink stage0_simulated
ln -s /data2/lansce/jr2514/stage0_simulated stage0_simulated
unlink stage0_bin_automated
ln -s /data2/lansce/jr2514/stage0_bin_automated stage0_bin_automated
unlink raw_data
ln -s /data2/lansce/jr2514/raw_data raw_data

cd Config
unlink DanceMap.txt
ln -s DanceMap_23Nov15.txt DanceMap.txt
unlink TMatrix.txt
ln -s TMatrix_57Fe.txt TMatrix.txt

#need to link directories in Alpha Calibrator
cd ../../DANCE_Alpha_Calibrator

unlink ParamOutput
ln -s /data2/lansce/jr2514/Alpha_Calibrator_ParamOutput ParamOutput
unlink RootOutput
ln -s /data2/lansce/jr2514/Alpha_Calibrator_RootOutput RootOutput
unlink quadratic_terms.txt
ln -s quadratic_terms_2018.txt quadratic_terms.txt

cd DANCE_Alpha_Database
unlink DANCE_Alpha_Database.root
ln -s DANCE_Alpha_Database_2018.root DANCE_Alpha_Database.root

#for calibrations and gates
#Calibrations is its own folder on the same level as DANCE_Analysis
cd ../../Calibrations
unlink calib_ideal.dat
ln -s calib_ideal_2018.dat calib_ideal.dat

cd ../DANCE_Analysis_Ang/Gates
unlink Alpha.dat
#ln -s Alpha_116470.dat Alpha.dat
#ln -s Alpha_112775.dat Alpha.dat
#the following is a test gate drawn from run 113391
#ln -s Alpha_109430.dat Alpha.dat 
ln -s Alpha_107982.dat Alpha.dat


unlink Gamma.dat
ln -s Gamma_109430.dat Gamma.dat

unlink Retrigger.dat
ln -s Retrigger_071918.dat Retrigger.dat
#ln -s Retrigger_113210.dat Retrigger.dat
unlink Pileup.dat
ln -s Pileup2018.dat Pileup.dat

cd ..
#Need to have a stage1.cfg file in here too
unlink stage1.cfg
ln -s cfg_files/stage1_Cd114.cfg stage1.cfg


cd /data2/lansce/jr2514
unlink raw_data
ln -s /data/lansce/7862-A/raw_data raw_data


