#!/bin/bash
top_path=/dybfs/users/xuhangkun/SimTAO/offline
source ${top_path}/SourceDesign/setup.sh
SourceName=${1}
mode=${2}
energy=${3}
version=${4}
#python3 utils/TaoUsefulDistribution.py -path ${top_path}/change_data/neutron_design/2Weight_Enclosure_Ref_0.95 -file ${SourceName}.root -source_name ${SourceName} -mode ${mode} -output ${top_path}/SourceDesign/data/2weight/${SourceName}.root -energy ${energy} -radiu 400
python3 utils/TaoUsefulDistribution.py -path ${top_path}/change_data/neutron_design/Nake -file ${SourceName}_v${version}.root -source_name ${SourceName} -mode ${mode} -output ${top_path}/SourceDesign/data/nake/${SourceName}_v${version}.root -energy ${energy} -radiu 400
