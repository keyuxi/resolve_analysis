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

# print("\nslide_dict: \n", slide_dict, "\n\n")


# --- Define Required Output --- #

rule all:
    input:
        expand(DATADIR + "segmentation_with_nuclei_prior_{slide}/{tile}/segmentation_counts.tsv", slide=slide_dict.keys(), tile=TILES)
        # DATADIR + "cellpose_32801-slideC/cellpose_mask_A2.tiff"

# --- Rules --- #

rule run_cellpose:
    input:
        os.path.join(DATADIR, "{slide}", "{slide}_{tile}_DAPI.tiff")
    output:
        DATADIR + "cellpose_{slide}/cellpose_mask_{tile}.tiff"
    threads:
        1
    params:
        diameter = "65",
        mask_threshold = "-0.1",
        cluster_memory = "32G"
    conda:
        "envs/cellpose_export.yml"
        # Had to manually install pytorch then export
        # because it runs out of memory on sherlock using pip
        # and need to tell it's cpuonly
    shell:
        """
        python resolve/run_cellpose.py \
            --diameter {params.diameter} --mask_threshold {params.mask_threshold} \
            -i {input} -o {output}
        """

rule run_baysor:
    input:
        molecules = DATADIR + "{slide}/{slide}_{tile}_results.txt",
        dapi_seg = DATADIR + "cellpose_{slide}/cellpose_mask_{tile}.tiff"
    output:
        DATADIR + "segmentation_with_nuclei_prior_{slide}/{tile}/segmentation_counts.tsv"
    threads:
        1
    params:
        scale = "45.0",
        scale_std = "50%",
        exclude_genes = "'Adra2a'",
        prior_segmentation_confidence = "0.5",
        outdir = DATADIR + "segmentation_with_nuclei_prior_{slide}/{tile}",
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