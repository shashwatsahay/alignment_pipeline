print(
    """
    \tAlignment pipeline for UMI based sequencing reads
    \tAuthor: Shashwat Sahay
    \tEmail: shashwat.sahay@charite.de
    \tVersion: 0.2.0

    """
)


def get_data_time():
    now = datetime.now()

    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    return "========= " + dt_string + " ========= "


print(get_data_time(), "Starting run")

#######################################################################
################# Sanity check on the input data ######################
#######################################################################

if "sample" not in config:
    raise ValueError("Sample name not provided")

if "pid" not in config:
    raise ValueError("Patient ID not provided")

if "metadata" not in config:
    raise ValueError("Metadata file not provided")
else:
    sep = ","
    if config["metadata"].split(".")[-1] in ["tsv", "txt"]:
        sep = "\t"
    elif config["metadata"].split(".")[-1] not in ["csv", "tsv", "txt"]:
        raise ValueError(
            "Can't resolve the separator for the metadata file from extension"
        )
    metadata = pd.read_csv(config["metadata"], dtype=str, sep=sep)
    metadata = metadata[
        (metadata["SAMPLE_TYPE"] == config["sample"])
        & (metadata["PATIENT_ID"] == config["pid"])
    ]
    if metadata.shape[0] == 0:
        raise ValueError(
            "No metadata found for the given sample: %s and patient ID: %s"
            % (config["sample"], config["pid"])
        )

##########################################################
###### Setting Sequencing type and related parameters ####
##########################################################

if "SeqType" not in config:
    raise ValueError("Sequencing type not provided")
else:
    seq_type = config["SeqType"]
    if seq_type not in ["WGS", "WES", "Panel"]:
        raise ValueError("Invalid sequencing type provided")
    if seq_type in ["WES", "Panel"]:
        if "target_regions" not in config:
            raise ValueError("Target regions not provided")
        else:
            target_regions = config["target_regions"]

        if "chrom_sizes" not in config:
            raise ValueError("Chromosome sizes not provided")
        else:
            chrom_sizes = config["chrom_sizes"]
        if "dict_genome" not in config:
            raise ValueError("Genome not provided")
        else:
            dict_genome = config["dict_genome"]

        if "bait_regions" not in config:
            print(
                get_data_time(),
                "Bait regions not provided will be calculated from target regions",
            )
        else:
            bait_regions = config["bait_regions"]


##########################################################
################ Setting Library kit #####################
##########################################################

if "library_prep_kit" in config:
    library_prep_kit = config["library_prep_kit"]
    if library_prep_kit == "":
        library_prep_kit = "Unknown"
else:
    library_prep_kit = "Unknown"


##########################################################
################## Setting Genome ########################
##########################################################

if "genome" not in config:
    raise ValueError("Genome not provided")
else:
    genome = config["genome"]

if "dbsnp" not in config:
    raise ValueError("dbSNP file not provided")
else:
    dbsnp = config["dbsnp"]


##########################################################
############# Setting working directory ##################
##########################################################

if "work_dir" not in config:
    raise ValueError("Work directory not provided")
else:
    wrkdir = Path(config["work_dir"])

if "log_dir" not in config:
    logdir = wrkdir / "logs"
else:
    logdir = Path(config["log_dir"])

if "scratch_dir" not in config:
    scratch_dir = tempfile.gettempdir()
else:
    scratch_dir = config[
        "scratch_dir"
    ]  # snakemake doesnt accept path in resources has to be a string
    if not os.path.exists(scratch_dir):
        print(get_data_time(), "Scratch directory does not exist")
        os.makedirs(scratch_dir)

print(get_data_time(), "Setting working directory to %s" % wrkdir)
print(get_data_time(), "Setting log directory to %s" % logdir)
print(get_data_time(), "Setting temp directory to %s" % scratch_dir)


#########################################################################
############# Setting variables related to trimming #####################
#########################################################################

if "trim_adapters" not in config or (not config["trim_adapters"]):
    print(get_data_time(), "Adapter trimming has been turned off")
    config["trim_adapters"] = False
    adapter_seq_r1 = None
    adapter_seq_r2 = None
else:
    if config["trim_adapters"]:
        if "adapter_forward_read" not in config:
            raise ValueError("Adapter sequence for R1 not provided")
        else:
            adapter_seq_r1 = config["adapter_forward_read"]
        if "adapter_reverse_read" not in config:
            raise ValueError("Adapter sequence for R3 not provided")
        else:
            adapter_seq_r2 = config["adapter_reverse_read"]


if "max_coverage" not in config:
    print(get_data_time(), "Setting default value for max_coverage to 10000")
    max_coverage = 10000
else:
    max_coverage = config["max_coverage"]


# Create a set of valid combination for run lane and read ids

LANE = metadata["LANE_NO"].unique().tolist()
RUN_ID = metadata["RUN_ID"].unique().tolist()


def filter_combinator(combinator, allow_list):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # print(wc_comb)
            if frozenset(wc_comb) in allow_list:
                yield wc_comb

    # print("hi")
    return filtered_combinator


allow_list = set()
for run_id in RUN_ID:
    for lane in metadata[metadata.RUN_ID == run_id]["LANE_NO"].unique().tolist():
        allow_list.add(
            frozenset(
                {
                    ("run_id", run_id),
                    ("sample", config["pid"] + "_" + config["sample"]),
                    ("lane", lane),
                }
            )
        )


filtered_product = filter_combinator(product, allow_list)

allow_list_2 = set()
for read in ["R1", "R2"]:
    for run_id in RUN_ID:
        for lane in metadata[metadata.RUN_ID == run_id]["LANE_NO"].unique().tolist():
            allow_list_2.add(
                frozenset(
                    {
                        ("run_id", run_id),
                        ("sample", config["pid"] + "_" + config["sample"]),
                        ("lane", lane),
                        ("read", read),
                    }
                )
            )

filtered_product_2 = filter_combinator(product, allow_list_2)

#######################################################
######## extra params #################################
### Since each tool can allow for mutliple parmas #####
### we can provide them directly in the config file ###
#######################################################

if "extra_cutadapt_params" in config:
    extra_cutadapt_params = config["extra_cutadapt_params"]
else:
    extra_cutadapt_params = ""


# print(metadata)
# print("blah")
# print(filtered_product)
# print("blah")
# print(allow_list)


def getFastq(wildcards, metadata, wrkdir, trim_adapters):
    if trim_adapters:
        fastq_r1 = [
            wrkdir
            / "fastq"
            / wildcards.run_id
            / "cutadapt"
            / (wildcards.sample + "_R1_" + wildcards.lane + "_trim.fastq.gz"),
            wrkdir
            / "fastq"
            / wildcards.run_id
            / "cutadapt"
            / (wildcards.sample + "_R2_" + wildcards.lane + "_trim.fastq.gz"),
        ]
    else:
        fastq_r1 = (
            metadata[
                (metadata["RUN_ID"] == wildcards.run_id)
                & (metadata["LANE_NO"] == wildcards.lane)
                & (metadata["SAMPLE_NAME"] == wildcards.sample)
            ]
            .sort_values("READ")["FASTQ_FILE"]
            .values
        )
    if len(fastq_r1) != 2:
        print(fastq_r1)
        raise ValueError("Expected to two file R1 and R2")

    return fastq_r1
