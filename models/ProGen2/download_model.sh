model=progen2-base

wget -P checkpoints/${model} https://storage.googleapis.com/sfr-progen-research/checkpoints/${model}.tar.gz
tar -xvf checkpoints/${model}/${model}.tar.gz -C checkpoints/${model}/