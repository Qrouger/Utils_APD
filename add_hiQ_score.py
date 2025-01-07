    cmd4 = f"singularity exec --no-home --bind result_homo_oligo:/mnt /home/emmanuel/Downloads/alpha-analysis_jax_0.4.sif run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
    os.system(cmd4)
    with open("predictions_with_good_interpae.csv", "r") as file1 :
        reader = csv.DictReader(file1)
        all_lines = "jobs,pi_score,iptm_ptm,hiQ_score\n"
        all_homo = dict()
        for row in reader :
            job = row['jobs']
            if 'homo' in job and row['pi_score'] != 'No interface detected' :
                if job not in all_homo.keys() :
                    all_homo[job] = (row['pi_score'],1,row)
                else :
                    sum_pi_score = float(all_homo[job][0]) + float(row['pi_score'])
                    sum_int = all_homo[job][1] + 1
                    all_homo[job] = (sum_pi_score,sum_int,row)
        for key in all_homo.keys() :
            row = all_homo[key][2]
            hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm is always same for each homo of same protein
            line = key+","+str(all_homo[key][0])+","+row['iptm_ptm']+","+str(hiQ_score)+"\n"
            all_lines = all_lines + line
    with open("predictions_with_good_interpae.csv", "w") as file2 :
        file2.write(all_lines)
