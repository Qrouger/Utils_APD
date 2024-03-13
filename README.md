Generate MSA-plot :

 	-Use APD_plot_prediction_quality.py
		python APD_plot_prediction_quality.py

	-Necessary file :
                If you have done mono-oligomer protéins : (Coverage, PAE and Plddt)
			$Pickle_name.pkl		
   			result_model_ptm_pred_{num}.pkl
   			ranking_debug.json
      
		If you have done only multi-oligomer protéins : (Coverage)
  			$Pickle_name.pkl		
			result_model_{num}_multimer_v3_pred_0.pkl 

more information about the script : https://blog.biostrand.ai/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2



Generate Distogram :

	-Use APD_plot_distogram.py
		python APD_plot_distogram.py $Path --pickle_name $Pickle_name --plot
	
	
	-Necessary file :
		$Pickle_name.pkl
		result_*.pkl		(file need to be gunzip)
more information about the script : https://github.com/clami66/dgram2dmap

Generate table of residue distance :

	-Use APD_calculate_residue_dist.py :
		python APD_calculate_residue_dist.py
  
  	-Necessary file :
		ranked_*.pdb

Allows to have the fundamental information present on the Distogram. A table of atoms distances.
