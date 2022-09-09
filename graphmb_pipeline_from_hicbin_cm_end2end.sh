#!/bin/bash


source ~/miniconda3/etc/profile.d/conda.sh


while getopts ":l:d:f:g:c:m:" opt; do
    case $opt in
        l)
            HICBIN_OUTPUT_LOCATION=$OPTARG;;
        g)
            GRAPHMB_OUT_PREFIX=$OPTARG;;
        d) 
            DEPTHS_FILE=$OPTARG;;
        f)
            FASTA=$OPTARG;;
        c)
            CHECKM_OUT_PREFIX=$OPTARG;;
        m)
            MESSAGE=$OPTARG;;
        \?)
            echo "Invalid parameter!"
            exit;;
    esac
done

NORM_CM="${HICBIN_OUTPUT_LOCATION}/contact_map_normalised.tsv"
RAW_CM="${HICBIN_OUTPUT_LOCATION}/contact_map_raw.tsv"

raw_cm="${HICBIN_OUTPUT_LOCATION}/contact_map_raw.npz"
norm_cm="${HICBIN_OUTPUT_LOCATION}/contact_map_normalized.npz"
names="${HICBIN_OUTPUT_LOCATION}/sequences_info.pkl"


conda activate ml

echo $HICBIN_OUTPUT_LOCATION


# echo "Processing raw contact map"
# python ~/tools/GNN_plus_HiC/create_contact_map_from_hicbin.py \
# -n $names -m $norm_cm -o $RAW_CM
# echo "DONE"


# echo "Preparing contact graph for raw contact map"

# hicbin_graphmb_input_raw="${HICBIN_OUTPUT_LOCATION}/graphmb_input_raw/"
# python ~/tools/GNN_plus_HiC/preprocess_files.py --graphmb -c $RAW_CM \
#                                                 -f $FASTA \
#                                                 -o $hicbin_graphmb_input_raw \
#                                                 --force 




graphmb_out_raw="${GRAPHMB_OUT_PREFIX}_raw/"
graphmb_out_norm="${GRAPHMB_OUT_PREFIX}_normalised/"

# echo "Running GRAPHMB on raw contact map"

# graphmb --outdir $graphmb_out_raw \
#         --depth $DEPTHS_FILE \
#         --assembly_name $FASTA \
#         --graph_file "${hicbin_graphmb_input_raw}/assembly_graph.gfa" \
#         --numcores 16 \
#         --minbin 500000 \
#         --mincontig 2000



# echo "Running CHECKM on raw GRAPHMB"


# conda activate checkm

# checkm lineage_wf --reduced_tree --pplacer_threads 16 --threads 16 -x fa \
#                 "${graphmb_out_raw}_bins/" \
#                 "${CHECKM_OUT_PREFIX}/graphmb_${MESSSAGE}_raw/"

conda activate ml

echo "Processing normalised contact map"
python ~/tools/GNN_plus_HiC/create_contact_map_from_hicbin.py \
-n $names -m $norm_cm -o $NORM_CM
echo "DONE"


echo "Preparing contact graph for normalised contact map"

hicbin_graphmb_input_norm="${HICBIN_OUTPUT_LOCATION}/graphmb_input_norm/"
python ~/tools/GNN_plus_HiC/preprocess_files.py --graphmb -c $NORM_CM \
                                                -f $FASTA \
                                                -o $hicbin_graphmb_input_norm \
                                                --force

echo "Running GRAPHMB on normalized contact map"

graphmb --outdir $graphmb_out_norm \
        --depth $DEPTHS_FILE \
        --assembly_name $FASTA \
        --graph_file "${hicbin_graphmb_input_norm}/assembly_graph.gfa" \
        --numcores 32 \
        --minbin 500000 \
        --mincontig 2000

conda activate checkm
echo "Running CHECKM on normalized GRAPHMB"
checkm lineage_wf --reduced_tree --pplacer_threads 16 --threads 16 -x fa \
                "${graphmb_out_norm}_bins/" \
                "${CHECKM_OUT_PREFIX}/graphmb_${MESSSAGE}_normalised/"
