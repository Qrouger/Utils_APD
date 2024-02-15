Generate MSA-plot :
  
	-Use APD_plot_prediction_quality.py
		python APD_plot_prediction_quality.py

	-necessary file :
		$Pickle_name.pkl
		result_model_multimer_v2_pred_{num}.pkl for multimer
		result_model_ptm_pred_{num}.pkl
		ranking_debug.json

Generate Distogram :

	-Use plot_distogram.py
		python plot_distogram.py $Path --pickle_name $Pickle_name --plot
	
	
	-necessary file :
		$Pickle_name.pkl
		result_*.pkl		(file need to be gunzip)
