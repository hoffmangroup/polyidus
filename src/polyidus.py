from argparse import ArgumentParser
from datetime import datetime
import os
import psutil
from utils import polyidusEngine


def report_memory(process_name):
    cur_time = str(datetime.now())
    process = psutil.Process(os.getpid())
    membytes = process.memory_info().rss
    if membytes > 1e9:
        cur_mem = "{} GBs".format(membytes / 1e9)
    else:
        cur_mem = "{} MBs".format(membytes / 1e6)
    print(
        "Memory usage at {} after {} is {}".format(
            cur_time, process_name, cur_mem))


def check_indices(hostindex, viralindex):
    index_dict = {"host": hostindex,
                  "viral": viralindex}
    for indexname, indexprefix in index_dict.items():
        indexfolder = os.path.dirname(indexprefix)
        indexstr = os.path.basename(indexprefix)
        if not os.path.exists:
            raise ValueError(
                "{} for {} index doesn't exist".format(
                    indexfolder, indexname))
        indexfiles = os.listdir(indexfolder)
        indexfiles = [each for each in indexfiles
                      if indexstr in each]
        if len(indexfiles) == 0:
            raise ValueError(
                "Index folder for {} was empty".format(
                    indexname))


def is_tool(name):
    from shutil import which
    return which(name) is not None


def check_other_programs():
    programs = ["samtools", "bedtools"]
    for program in programs:
        if not is_tool(program):
            raise ValueError(
                "{} doesn't exist in PATH".format(program))


def check_aligner(aligner):
    aligner_exists = is_tool(aligner)
    if not aligner_exists:
        raise ValueError(
            "Aligner {} doesn't exist in PATH".format(aligner))


if __name__ == "__main__":
    parser = ArgumentParser(
        '''
        Identify viral integration sites in FASTQ files.
        Requires a host reference genome index,
        a viral reference genome index, and the FASTQ file(s)
        of the experiment to identify the exact
        integration sites.
        Make sure you have generated genome indices for
        both the host and virus genomes.
        You can use bwa instead of bowtie2 as well.
        Make sure the aligner (bowtie2 or bwa), samtools,
        and bedtools also exist in $PATH.''',
        epilog='''
        Citation: Karimzadeh M., Arlidge C., Rostami A.,
        Lupien M., Bratman S., and Hoffman M. M.
        Integration of human papillomavirus into the human
        genome promotes cancer by modifying local
        chromatin and transcription.''')
    parser.add_argument(
        "hostindex",
        help="Path to host index prefix (either bwa or bowtie2)")
    parser.add_argument(
        "viralindex",
        help="Path to viral index prefix (either bwa or bowtie2)")
    parser.add_argument(
        "--fastq",
        nargs="*",
        help="Path to fastq files (if more than 1, assumes second "
        "is the second pair.")
    parser.add_argument(
        "--outdir",
        default=os.path.join(os.getcwd(), "polyidusOutput"),
        help="Path to output folder")
    parser.add_argument(
        "--aligner",
        choices=["bwa", "bowtie2"],
        default="bowtie2",
        help="Choose from bwa or bowtie2 (default)")
    args = parser.parse_args()
    report_memory("initialization")
    check_indices(args.hostindex, args.viralindex)
    check_aligner(args.aligner)
    check_other_programs()
    os.makedirs(args.outdir, exist_ok=True)
    polyidusObj = polyidusEngine(
        args.hostindex, args.viralindex,
        args.fastq, args.outdir, args.aligner)
    polyidusObj.align_files()
    report_memory("aligning virus and host fastq files")
    polyidusObj.find_approximate_integrations()
    report_memory("initial investigation of BAM files")
    polyidusObj.find_exact_integrations()
    report_memory("saving exact integration sites")
