
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ml

CORES="32"

while getopts ":c:p:s:" opt; do
    case $opt in
        c)
            CONFIG_DIR=$OPTARG;;
        p)
            CORES=$OPTARG;;

        s)
            SNAKEFILE=$OPTARG;;

        \?)
            echo "Invalid parameter!"
            exit;;
    esac
done





# for config in "${CONFIG_DIR}/*";
#     do
#         echo $config

# done

echo "Running pipelines with configs from dir: ${CONFIG_DIR}"

cd $CONFIG_DIR
for config in *;
    do

        configfile="${CONFIG_DIR}/${config}"
        cd ..
        echo "Using snakemake pipeline with config ${configfile}"
        snakemake --configfile $configfile --cores $CORES --jobs $CORES --snakefile $SNAKEFILE
        cd $CONFIG_DIR 

        echo "DONE! :)";
done