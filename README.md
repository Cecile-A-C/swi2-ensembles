# swi2-ensembles

Scripts and files for the ensemble-based history matching and multiobjective optimization under uncertainty described in 
"Cécile Coulon, Jeremy T. White, Alexandre Pryet, Laura Gatel & Jean-Michel Lemieux, A multi-objective, ensemble-based approach for pumping optimization in an island aquifer considering parameter, observation and climate uncertainty"

## File organization

## 1) ies

1.1) *pestpp_setup_ies_v4.py* used the files in **input_setup_files** to create the **send_to_ies_v4** folder and its contents
1.2) PESTPP-IES was run on a cluster using the files in **send_to_ies_v4**
1.3) the results of the IES run were saved in **results_v4**
1.4) *analyze_ies_outputs.py* analyzed the results from the IES run and plotted figures 6, 7 of Coulon et al.

## 2) postproc

2.1) *postproc_v3.py* generated the sea level and recharge ensembles using the output files of 1) and *SWB2_recharge_scenarios.csv* and plotted figures 3, 8 of Coulon et al.

## 3) swp

3.1) *pstpp_setup_swp_v2.py* used the output files of 1) and 2) to create the **climatechange** and **noclimatechange** folders in **send_to_swp_v2**
3.2) PESTPP-SWP was run on a cluster using the files in **send_to_swp_v2**
3.3) the results of the SWP runs were saved in **results_v2**
3.4) *compare_swp_runs.py* compared the results from the SWP runs and plotted figure 9 of Coulon et al.

## 4) mou

4.1) *pstpp_setup_mou_v3.py* used the output files of 3) to create the **climatechange** and **noclimatechange** folders in **send_to_mou_v4**
4.2) PESTPP-MOU was run on a cluster using the files in **send_to_mou_v4**
4.3) the results of the MOU runs were saved in **results_v4**
4.4) *compare_mou_opt.py* compared the results from the MOU runs and plotted figures 10, 11 of Coulon et al.


## Typical folders used to run PEST++

All the folders used for PEST++ runs had a similar structure: (**send_to_ies_v4**, **send_to_swp_v2**, **send_to_mou_v4**)

1) **model_ouputs**: folder containing the model output files
2) **model_params**: folder containing the model input files
3) **pest_files**: folder containing the PEST instruction and template files
4) *preproc_IDM_20m_qgis_v6.bin*: binary file containing spatial inputs
5) *model_run_steady.py* or *model_run_opt.py*: script for the model run
6) *functions_model_run.py*: script containing functions used by *model_run.py*
7) *pestpp_ies.pst*, *pestpp_swp.pst* or *pestpp_mou.pst*: PEST control file

With additional files depending on the PEST++ program used.
