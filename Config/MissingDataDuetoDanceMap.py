import ROOT
import os
no_data_runnums = []

for runnum in range(107923,108922): 
    pathtofile = f'/home/jr2514/DANCE/DANCE_Analysis/stage0_root/Stage0_Histograms_Run_{runnum}_500ns_CW_0ns_CBT_0ns_DEBT.root'
    if not os.path.isfile(pathtofile):
        print(f"File does not exist: {pathtofile}")
    else:  
        file = ROOT.TFile(pathtofile)
        hID = file.Get("ID")

        bin200 = hID.FindBin(200)
        bin201 = hID.FindBin(201)

        # Get the bin content
        content200 = hID.GetBinContent(bin200)
        content201 = hID.GetBinContent(bin201)

        if content200 != 0 and content201 != 0:
            no_data_runnums.append(runnum)

with open('runs_bothdata_map_11_15_201', 'w') as file:
    for element in no_data_runnums:
        file.write(f'{element}\n')
    
