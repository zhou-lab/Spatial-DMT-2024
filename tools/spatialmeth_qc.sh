#!/usr/bin/env bash
################################################################################
##
## Quality Control script for BISCUIT output
##
## Output from this script can be fed into MultiQC to produce a nice HTML output
## showing the different BISCUIT QC metrics
##
## Notes:
##   1.) biscuit, samtools, bedtools, and awk all must be in PATH for
##       script to work
##
## Created by:
##   Wanding Zhou
##
## Creation date:
##   May 2019
##
## Update notes:
##   Dec 2019 -
##     - Clean up code to make more readable
##     - Catch empty files, alert user, and remove files
##   Mar 2020 -
##     - Reworking of QC script
##         -- Removed repeats of functions and code
##         -- Removed python references for using BASH throughout
##         -- All sorts now have parallel option
##         -- Consistent tmp file naming scheme
##         -- <unset> is now <unused>
##         -- Assorted other changes
##     - Bug fix in Coefficient of Variation calculation
##   Apr 2020 -
##     - Refactoring QC script to reduce time spent running
##   May 2020 -
##     - Changing flags as necessary for updates to subcommand flags
##     - Adding PATH check for GNU parallel
##   Jun 2020 -
##     - Adding PATH check for GNU awk
##   Jul 2021 -
##     - Create qc subcommand and integrate into QC.sh
##   Jan 2023 -
##     - Fix bug so new versions of GNU parallel work
##     - Remove GNU parallel dependency
##
################################################################################

set -euo pipefail
# Improves sort/awk performance
export LC_ALL=C

# Check for biscuit, samtools, bedtools, awk in PATH
function check_path {
  if [[ `which biscuit 2>&1 > /dev/null` ]]; then
      >&2 echo "biscuit does not exist in PATH"
      exit 1
  else
      >&2 echo "Using biscuit found at: `which biscuit`"
  fi
  if [[ `which samtools 2>&1 > /dev/null` ]]; then
      >&2 echo "samtools does not exist in PATH"
      exit 1
  else
      >&2 echo "Using samtools found at: `which samtools`"
  fi
  if [[ `which bedtools 2>&1 > /dev/null` ]]; then
      >&2 echo "bedtools does not exist in PATH"
      exit 1
  else
      >&2 echo "Using bedtools found at: `which bedtools`"
  fi
  if [[ `which awk 2>&1 > /dev/null` ]]; then
      >&2 echo "awk does not exist in PATH"
      exit 1
  else
      if awk --version | grep -q GNU; then
          >&2 echo "Using GNU awk found at: `which awk`"
      else
          >&2 echo "It doesn't appear you are using GNU awk"
          >&2 echo "Try adding GNU awk at the front of PATH"
          exit 1
      fi
  fi
}

# Check that certain variables have been set and files exist
function check_variables {
    VARS="
    BISCUIT_CPGS
    BISCUIT_TOPGC
    BISCUIT_BOTGC
    in_vcf
    "

    for var in ${VARS}; do
        if [[ ${!var} != "<unused>" && ! -f ${!var} ]]; then
            >&2 echo "${var}: ${!var} does not exist"
            exit 1
        fi
    done
}

# Wait for jobs, exit if they're unsuccessful
# Via: https://stackoverflow.com/questions/1131484/wait-for-bash-background-jobs-in-script-to-be-finished/
function wait_for_jobs {
    while true; do
	wait -n || {
	    code="$?"
	    ([[ $code = "127" ]] && exit 0 || exit "$code")
	    break
	}
    done;
}

# Workhorse function for processing BISCUIT QC
function biscuitQC {
    # Simple check for necessary command line tools
    check_path

    # Check variables and their associated files exist
    #check_variables

    # Create outdir if it does not exist
    if [ ! -d ${outdir} ]; then
        mkdir -p ${outdir}
    fi

    >&2 echo -e "Starting BISCUIT QC at `date`\n"

    # MAPQ, Insert Size, Duplicate Rates, Strand Info, and Read-averaged conversion
    if [[ "${single_end}" == true ]]; then
        biscuit qc -s ${genome} ${in_bam} ${outdir}/${sample}
    else
        biscuit qc ${genome} ${in_bam} ${outdir}/${sample}
    fi

    ## NOTE: covdist + cv_table have been deliberately dropped from this script.
    ## The original implementation derived them via `paste <(BISCUIT_CPGS) <(yame
    ## subset | yame unpack)`, which silently zips the reference CpG BED against
    ## the per-sample/per-cell .cg row-by-row. When the contig set in BISCUIT_CPGS
    ## doesn't exactly match the contig set the .cg was built over (alt contigs,
    ## chrUn, chrM), the rows shift and depths are reported for the wrong CpGs.
    ## Sample-level coverage/depth are now computed from `yame summary` in the
    ## qc_curated_metrics rule -- no row-order assumption, contig-agnostic.

    if [[ -f ${in_vcf} ]]; then
        echo "BISCUITqc Conversion Rate by Base Average Table" \
            > ${outdir}/${sample}_totalBaseConversionRate.txt
        biscuit vcf2bed -e -t c ${in_vcf} | \
        awk '{ beta_sum[$6] += $8; beta_cnt[$6] += 1 } END {
            print "CA\tCC\tCG\tCT"
            if ( beta_cnt["CA"] < 20 ) {
                ca_frac = -1
            } else {
                ca_frac = beta_sum["CA"] / beta_cnt["CA"]
            }
            if ( beta_cnt["CC"] < 20 ) {
                cc_frac = -1
            } else {
                cc_frac = beta_sum["CC"] / beta_cnt["CC"]
            }
            if ( beta_cnt["CG"] < 20 ) {
                cg_frac = -1
            } else {
                cg_frac = beta_sum["CG"] / beta_cnt["CG"]
            }
            if ( beta_cnt["CT"] < 20 ) {
                ct_frac = -1
            } else {
                ct_frac = beta_sum["CT"] / beta_cnt["CT"]
            }
            print ca_frac"\t"cc_frac"\t"cg_frac"\t"ct_frac
        }' >> ${outdir}/${sample}_totalBaseConversionRate.txt
    else
        >&2 echo -ne "Input VCF not supplied. "
        >&2 echo -ne "${sample}_totalBaseConversionRate.txt will not be generated."
    fi

    ###################################
    ## Remove temporary files
    ###################################
    if [[ ! "${keep_tmp}" == true ]]; then
        rm -f ${outdir}/${sample}*.tmp.*
    fi

    >&2 echo -e "\nFinished BISCUIT QC at `date`"
}

# Print helpful usage information
function usage {
    >&2 echo -e "\nUsage: QC.sh [-h,--help] [-s,--single-end] [-v,--vcf] [-o,--outdir] [-k,--keep-tmp-files] assets_directory genome sample_name in_bam\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tassets_directory    : Path to assets directory"
    >&2 echo -e "\tgenome              : Path to reference FASTA file used in alignment"
    >&2 echo -e "\tsample_name         : Prefix of output files"
    >&2 echo -e "\tinput_bam           : Aligned BAM from BISCUIT\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help           : Print help message and exit"
    >&2 echo -e "\t-s,--single-end     : Input BAM is from single end data [DEFAULT: Assumes paired-end]"
    >&2 echo -e "\t-v,--vcf            : Path to VCF output from BISCUIT [DEFAULT: <unused>]"
    >&2 echo -e "\t-o,--outdir         : Output directory [DEFAULT: BISCUITqc]"
    >&2 echo -e "\t-k,--keep-tmp-files : Flag to keep temporary files for debugging [DEFAULT: Delete files]\n"
}

# Initialize default values for optional inputs
in_vcf="<unused>"
outdir="BISCUITqc"
keep_tmp=false
single_end=false

# Process command line arguments
OPTS=$(getopt \
    --options hsv:o:k \
    --long help,single-end,vcf:,outdir:,keep-bed-files \
    --name "$(basename "$0")" \
    -- "$@"
)
eval set -- ${OPTS}

while true; do
    case "$1" in
        -h|--help )
            usage
            exit 0
            ;;
        -s|--single-end )
            single_end=true
            shift
            ;;
        -v|--vcf )
            in_vcf="$2"
            shift 2
            ;;
        -o|--outdir )
            outdir="$2"
            shift 2
            ;;
        -k|--keep-tmp-files )
            keep_tmp=true
            shift
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure there are the correct number of inputs
if [[ $# -ne 4 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
assets=$1
genome=$2
sample=$3
in_bam=$4

# Do some checks on the given inputs
if [[ ! -d "$assets" ]]; then
    >&2 echo "Assets directory missing: $assets"
    exit 1
fi

if [[ ! -f "${genome}.fai" ]]; then
    >&2 echo "Cannot locate fai-indexed reference: ${genome}.fai"
    >&2 echo "Please provide a viable path to the reference genome FASTA file."
    exit 1
fi

if [[ ! -f "${in_bam}" ]]; then
    >&2 echo "Cannot locate aligned BAM: ${in_bam}"
    >&2 echo "Please provide an existing aligned BAM."
    exit 1
fi

# Set variables for supplementary BED files
BISCUIT_CPGS="${assets}/cpg.bed.gz"
BISCUIT_TOPGC="${assets}/windows100bp.gc_content.top10p.bed.gz"
BISCUIT_BOTGC="${assets}/windows100bp.gc_content.bot10p.bed.gz"

>&2 echo "## Running BISCUIT QC script with following configuration ##"
>&2 echo "=============="
>&2 echo "Sample Name        : ${sample}"
>&2 echo "Input BAM          : ${in_bam}"
>&2 echo "Input VCF          : ${in_vcf}"
>&2 echo "Output Directory   : ${outdir}"
>&2 echo "Assets Directory   : ${assets}"
>&2 echo "Reference          : ${genome}"
>&2 echo "Keep *.tmp.* files : ${keep_tmp}"
>&2 echo "CG file            : ${cg_file}"
>&2 echo "Barcode            : ${barcode}"
>&2 echo "Single-end data    : ${single_end}"
>&2 echo "CPGS               : ${BISCUIT_CPGS}"
>&2 echo "TOPGC              : ${BISCUIT_TOPGC}"
>&2 echo "BOTGC              : ${BISCUIT_BOTGC}"
>&2 echo "=============="

biscuitQC
