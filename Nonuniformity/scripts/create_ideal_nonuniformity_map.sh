#!/bin/bash
energies=(1 4 8)
particle=${1}
energy=${energies[${2}]}

case ${particle} in
"electron")
python3 ${TAO_CALIB_PATH}/Nonuniformity/utils/create_ideal_map.py \
--output ${TAO_CALIB_PATH}/Nonuniformity/data/map/ideal_nonuniformity_${particle}_${energy}MeV.pkl \
--input_dir ${TAO_CALIB_PATH}/change_data/nonuniformity/electron \
--energy ${energy} \
--particle electron
;;
"gamma")
python3 ${TAO_CALIB_PATH}/Nonuniformity/utils/create_ideal_map.py \
--output ${TAO_CALIB_PATH}/Nonuniformity/data/map/ideal_nonuniformity_${particle}_${energy}MeV.pkl \
--input_dir ${TAO_CALIB_PATH}/change_data/nonuniformity/gammas \
--energy ${energy} \
--particle gamma
;;
"positron")
python3 ${TAO_CALIB_PATH}/Nonuniformity/utils/create_ideal_map.py \
--output ${TAO_CALIB_PATH}/Nonuniformity/data/map/ideal_nonuniformity_${particle}_${energy}MeV.pkl \
--input_dir ${TAO_CALIB_PATH}/change_data/nonuniformity/positron \
--energy ${energy} \
--particle positron
;;
esac

