#!/usr/bin/bash                                                                                                                       
cd /storage/cylin/bin/pipeline/
source activate py27_anaconda
python2.7 /storage/cylin/bin/pipeline/dynamicEnhancer.py -g hg19 -d /storage/cylin/grail/projects/Test_aquas/HS_Kronos_CHIP_data_table.txt -n 22rv1_Prostate_Treat_AR,22rv1_Prostate_Control_AR -r /storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/rose/22rv1_Prostate_Treat_AR_ROSE,/storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/rose/22rv1_Prostate_Control_AR_ROSE -o /storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/dynamic/Treat_+_control_AR -m -a

#python2.7 /storage/cylin/bin/pipeline/dynamicEnhancer.py -g hg19 -d /storage/cylin/grail/projects/Test_aquas/HS_Kronos_CHIP_data_table.txt -n 22rv1_Prostate_Treat_FOXA1,22rv1_Prostate_Control_FOXA1 -r /storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/rose/22rv1_Prostate_Treat_FOXA1_ROSE,/storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/rose/22rv1_Prostate_Control_FOXA1_ROSE -o /storage/cylin/grail/projects/Test_aquas/Kronos_prostate_cancer/dynamic/Treat_+_control_FOXA1 -m -a
