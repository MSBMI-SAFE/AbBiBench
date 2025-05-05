for name in 2fjg 1mlc 1n8z aayl49_ML aayl51 aayl49 3gbn 1mhp 4fqi
do
  python models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order HLA
  python models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order AHL
  python models/ProtGPT2/get_model_log_likelihood.py --name $name --chain_order LAH
done