for name in 3gbn 2fjg aayl49 aayl49_ML aayl51 1mlc 1n8z 4fqi 1mhp
do
  python ../models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order HLA
  python ../models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order AHL
  python ../models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order LAH
done