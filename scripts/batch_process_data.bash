#!/bin/bash

mapfile -t files < <(ls)
for file in "${files[@]}"; do
    if [[ "$file" == *.py ]]; then
        echo "* $file"
    fi
done

# ns-train spirulae --data data/adr_food/251211_chao_3 --pipeline.model.apply-loss-for-mask True --pipeline.model.randomize_background False --pipeline.model.mcmc_cap_max 100000 --pipeline.model.sh-degree 3 --viewer.quit_on_train_completion True
