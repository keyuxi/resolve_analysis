import os
# from scripts.util import *

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_20230303.yaml"
#######################################

# --- Define Global Variables --- #

DATADIR = config["dataDir"]
EXPDIR = os.path.normpath(DATADIR + "/../") + "/"
SLIDES = config["slides"]

# hardcoded tile numbers
TILES = expand("{letter}{number}-1", letter=["A", "B", "C", "D"], number=["1", "2"])
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
        # expand(DATADIR + "segmentation_with_strong_nuclei_prior_{slide}/{tile}/segmentation_counts.tsv", slide=slide_dict.keys(), tile=TILES)
        # DATADIR + "segmentation_32801-slideC/C2/segmentation_counts.tsv"
        expand("/home/users/kyx/oak/data/spatial/segmentation_P22344_slideF4_20230221/{tile}/segmentation_counts.tsv", tile=TILES)
        #"/home/users/kyx/oak/data/spatial/segmentation_P22344_slideF4_20230221/A2-1/segmentation_counts.tsv"

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
        cluster_memory = "32G",
        singularity_image = "/home/users/kyx/oak/software/singularity/cellpose_2.1.1_cv1.sif"
    # singularity:
        # "/home/users/kyx/oak/software/singularity/cellpose_2.1.1_cv1.sif"
    # conda:
        # "envs/cellpose_export.yml"
        # Had to manually install pytorch then export
        # because it runs out of memory on sherlock using pip
        # and need to tell it's cpuonly
    shell:
        """
        singularity exec {params.singularity_image} \
            python resolve/run_cellpose.py \
             --diameter {params.diameter} --mask_threshold {params.mask_threshold} \
             -i {input} -o {output}
        """

rule run_baysor:
    input:
        molecules = DATADIR + "{slide}/{slide}_{tile}_results.txt",
        dapi_seg = DATADIR + "cellpose_{slide}/cellpose_mask_{tile}.tiff"
    output:
        DATADIR + "segmentation_{slide}/{tile}/segmentation_counts.tsv"
    threads:
        1
    params:
        scale = "45.0",
        scale_std = "50%",
        exclude_genes = "'Adra2a'",
        prior_segmentation_confidence = "0.9",
        outdir_mount = "/data/segmentation_{slide}/{tile}",
        molecules_mount = "/data/{slide}/{slide}_{tile}_results.txt",
        dapi_seg_mount = "/data/cellpose_{slide}/cellpose_mask_{tile}.tiff",
        datadir = DATADIR,
        cluster_memory = "16G",
        cluster_time = "02:00:00"
    shell:
        """
        scripts/add_header.sh {input.molecules}
        export SINGULARITY_BIND="{params.datadir}:/data"
        singularity exec --cleanenv {config[baysor]} /Baysor/bin/Baysor run\
            -s {params.scale} --scale-std={params.scale_std} -c config/config.toml \
            -p --save-polygons=geojson \
            --force-2d \
            --prior-segmentation-confidence {params.prior_segmentation_confidence} \
            -o {params.outdir_mount} {params.molecules_mount} \
            {params.dapi_seg_mount}
        """