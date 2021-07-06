#!/bin/bash

mode=${1}
radius_cut=150
case ${mode} in
"normal")
# most normal case
python3 generate_gamma_map.py --event_num 50000
;;
"radius_cut")
# most normal case
python3 generate_gamma_map.py --fit_can_file ./data/fit_gamma_radius_cut.root \
--gamma_nonuniformity_file ./data/gamma_calibed_radius_cut_nonuniformity.root \
--gamma_nonuniformity_map ./data/map/gamma_calibed_radius_cut_nonuniformity.csv \
--radius_cut ${radius_cut} --event_num 50000
;;
"dead_sipm")
# open the switch of dead sipm
python3 generate_gamma_map.py --fit_can_file ./data/fit_gamma_sipmdead.root \
--gamma_nonuniformity_file ./data/gamma_calibed_sipmdead_nonuniformity.root \
--gamma_nonuniformity_map ./data/map/gamma_calibed_sipmdead_nonuniformity.csv \
--event_num 50000 --open_sipm_dead
;;
"dead_sipm_radius_cut")
# open the switch of dead sipm
python3 generate_gamma_map.py --fit_can_file ./data/fit_gamma_sipmdead_radius_cut.root \
--gamma_nonuniformity_file ./data/gamma_calibed_sipmdead_radius_cut_nonuniformity.root \
--gamma_nonuniformity_map ./data/map/gamma_calibed_sipmdead_radius_cut_nonuniformity.csv \
--event_num 50000 --open_sipm_dead --radius_cut ${radius_cut}
;;
esac
