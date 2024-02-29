Generate MSA-plot :
  
	-Use APD_plot_prediction_quality.py
		python APD_plot_prediction_quality.py

	-Necessary file :
		$Pickle_name.pkl
		result_model_multimer_v2_pred_{num}.pkl for multimer
		result_model_ptm_pred_{num}.pkl
		ranking_debug.json
more information about the script : https://blog.biostrand.ai/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2
This script can be used only after done a mono-holigomer of the protéins


Generate Distogram :

	-Use APD_plot_distogram.py
		python plot_distogram.py $Path --pickle_name $Pickle_name --plot
	
	
	-Necessary file :
		$Pickle_name.pkl
		result_*.pkl		(file need to be gunzip)
more information about the script : https://github.com/clami66/dgram2dmap

Generate table of residue distance :

	-Use APD_calculate_residue_dist.py :
		python APD_calculate_residue_dist.py
  
  	-Necessary file :
		ranked_*.pdb
