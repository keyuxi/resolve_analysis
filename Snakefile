import os
# from scripts.util import *

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_20220116.yaml"
#######################################

# --- Define Global Variables --- #

DATADIR = config["dataDir"]
EXPDIR = os.path.normpath(DATADIR + "/../") + "/"
SLIDES = config["slides"]

# hardcoded tile numbers
TILES = expand("{letter}{number}", letter=["A", "B", "C", "D"], number=["1", "2"])
# master metadata dictionary for all tiles
slide_dict = {slide: 
    {"dir": os.path.join(DATADIR, slide),
     "results": expand(os.path.join(DATADIR, slide, "%s_{tile}_results.txt"%slide), tile=TILES),
     "dapi": expand(os.path.join(DATADIR, slide, "%s_{tile}_DAPI.tiff"%slide), tile=TILES)} 
    for slide in SLIDES}

print("\nslide_dict: \n", slide_dict, "\n\n")


# --- Define Required Output --- #

rule all:
    input:
        # expand(DATADIR + "segmentation_{slide}/{tile}/segmentation_counts.tsv", slide=slide_dict.keys(), tile=TILES)
        DATADIR + "cellpose_32801-slideC/cellpose_mask_A1.tiff"

# --- Rules --- #

rule run_cellpose:
    input:
        os.path.join(DATADIR, "{slide}", "%s_{tile}_DAPI.tiff"%slide)
    output:
        DATADIR + "cellpose_{slide}/cellpose_mask_{tile}.tiff"
    params:
        diameter = 65,
        mask_threshold = 0.2
        cluster_memory = "32G"
    conda:
        "envs/cellpose.yml"
    shell:
        """
        python scripts/cellpose.py \
            --diameter {params.diameter} --mask_threshold {mask_threshold} \
            -i {input} -o {output}
        """

rule run_baysor:
    input:
        molecules = DATADIR + "{slide}/{slide}_{tile}_results.txt",
        dapi_seg = DATADIR + "cellpose_{slide}/cellpose_mask_{tile}.tif"
    output:
        DATADIR + "segmentation_{slide}/{tile}/segmentation_counts.tsv"
    threads:
        1
    params:
        scale = "45.0",
        scale_std = "50%",
        exclude_genes = "'Adra2a'",
        prior_segmentation_confidence = "0.2",
        outdir = DATADIR + "segmentation_{slide}/{tile}",
        cluster_memory = "16G",
        cluster_time = "02:00:00"
    shell:
        """
        scripts/add_header.sh {input.molecules}
        {config[baysor]} run \
            -s {params.scale} --scale-std={params.scale_std} -c config/config.toml \
            -p --save-polygons=geojson --exclude-genes={params.exclude_genes} --force-2d \
            --prior-segmentation-confidence {params.prior_segmentation_confidence} \
            -o {params.outdir} {input.molecules} \
            {input.dapi_seg}
        """