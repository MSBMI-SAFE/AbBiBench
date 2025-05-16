for model in progen2-small progen2-medium progen2-base progen2-large
do
    for chain_order in HLA AHL LAH
    do
        for name in 3gbn 2fjg aayl49 aayl49_ML aayl51 1mlc 1n8z 4fqi 1mhp
        do
            python ../models/ProGen2/get_model_log_likelihood.py --name $name --model $model --chain_order $chain_order --with_structure
        done
    done
done