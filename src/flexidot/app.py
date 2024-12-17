from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import sys
import time

from matplotlib import rc as mplrc
from matplotlib import rcParams
import pylab as P

from flexidot._version import __version__
from flexidot.plotting import selfdotplot, pairdotplot, polydotplot
from flexidot.utils.logs import init_logging
from flexidot.utils.file_handling import read_gff_color_config
from flexidot.utils.utils import time_track

# Matplotlib settings
P.switch_backend(
    "agg"
)  # bugfix for _tkinter.TclError on CentOs 7 servers, see Github Issue #5

# Font settings
mplrc("pdf", fonttype=42, compression=0)

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = [
    "Helvetica",
    "Verdana",
    "Tahoma",
    "DejaVu Sans",
    "Droid Sans Mono",
    "Sans",
    "Liberation",
    "Ubuntu",
    "Arial",
]


def parse_args():
    parser = ArgumentParser(
        prog="flexidot",
        description="FlexiDot: Flexible dotplot generation tool",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # Input/output arguments
    parser.add_argument(
        "-i",
        "--infiles",
        required=True,
        nargs="+",
        help="Input fasta files (fasta file name or space-separated file list.)",
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        default="flexidot_output",
        help="File prefix to be added to the generated filenames.",
    )

    # Collage arguments
    parser.add_argument(
        "-c",
        "--collage",
        action="store_true",
        default=False,
        help="Combine multiple dotplots in a collage.",
    )
    parser.add_argument(
        "--n_col",
        type=int,
        default=4,
        help="Number of columns per page (if collage is ON.",
    )
    parser.add_argument(
        "--n_row",
        type=int,
        default=5,
        help="Number of rows per page (if collage is ON).",
    )

    # File type
    parser.add_argument(
        "-f",
        "--filetype",
        choices=["png", "pdf", "svg"],
        default="png",
        help="Output file format: png, pdf, svg",
    )

    # Sorting
    parser.add_argument(
        "-s",
        "--sort",
        action="store_true",
        default=False,
        help="Sort sequences alphabetically by name.",
    )

    # Calculation parameters
    parser.add_argument(
        "-k",
        "--wordsize",
        type=int,
        default=10,
        help="Wordsize (kmer length) for dotplot comparison.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        action="append",
        choices=["0", "1", "2"],
        help="Mode of FlexiDot dotplotting. 0 = self [default], 1 = paired, 2 = poly (matrix with all-against-all dotplots). Call -m multiple times to run multiple modes.",
    )
    parser.add_argument(
        "-t",
        "--type_seq",
        choices=["aa", "nuc"],
        default="nuc",
        help="Biological sequence type: aa (amino acid) or nuc (nucleotide).",
    )
    parser.add_argument(
        "-w",
        "--wobble_conversion",
        action="store_true",
        default=False,
        help="Ambiguity handling for relaxed matching.",
    )
    parser.add_argument(
        "-S",
        "--substitution_count",
        type=int,
        default=0,
        help="Number of substitutions allowed per window.",
    )
    parser.add_argument(
        "-r",
        "--norev",
        action="store_true",
        default=False,
        help="Do not calculate reverse complementary matches (only for nucleotide sequences.)",
    )
    parser.add_argument(
        "-O",
        "--only_vs_first_seq",
        action="store_true",
        default=False,
        help="Limit pairwise comparisons to the 1st sequence only (if plotting mode=1 paired.)",
    )

    # Graphic formatting
    parser.add_argument("-A", "--line_width", type=float, default=1, help="Line width")

    parser.add_argument("-B", "--line_col_for", default="black", help="Line color")

    parser.add_argument(
        "-C", "--line_col_rev", default="green", help="Reverse line color"
    )

    parser.add_argument(
        "-D",
        "--x_label_pos",
        choices=["top", "bottom"],
        default="top",
        help="Position of the X-label. Default: 'top'",
    )

    parser.add_argument("-E", "--label_size", type=int, default=10, help="Font size")

    parser.add_argument(
        "-F",
        "--spacing",
        type=float,
        default=0.04,
        help="Spacing between dotplots (if plotting mode=2 polyplot).",
    )

    parser.add_argument(
        "-L",
        "--length_scaling",
        action="store_true",
        default=False,
        help="Scale plot size for pairwise comparison.",
    )

    parser.add_argument(
        "-M",
        "--mirror_y_axis",
        action="store_true",
        default=False,
        help="Flip y-axis (bottom-to-top or top-to-bottom)",
    )

    parser.add_argument("-P", "--plot_size", type=int, default=10, help="Plot size")

    parser.add_argument(
        "-R",
        "--representation",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Region of plot to display. Only if plotting mode is 2: polyplot\n\
                                    0 = full [default]\n\
                                    1 = upper\n\
                                    2 = lower",
    )

    parser.add_argument(
        "-T",
        "--title_length",
        type=int,
        default=50,
        help="Limit title length for comparisons. Default: 50 characters",
    )

    # GFF shading
    parser.add_argument(
        "-g",
        "--gff",
        nargs="+",
        default=None,
        help="GFF3 files for markup in self-dotplots. Provide a space-delimited list of GFF files.",
    )
    parser.add_argument(
        "-G",
        "--gff_color_config",
        default=None,
        help="Config file for custom GFF shading.",
    )

    # longest common subsequence (LCS) shading
    parser.add_argument(
        "-x",
        "--lcs_shading",
        action="store_true",
        default=False,
        help="Shade subdotplot based on longest common subsequence (LCS).",
    )

    parser.add_argument(
        "-X",
        "--lcs_shading_num",
        type=int,
        default=5,
        help="Number of shading intervals.",
    )

    parser.add_argument(
        "-y",
        "--lcs_shading_ref",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Reference for LCS shading.\n\
                        0 = maximal LCS length [default]\n\
                        1 = maximally possible length (length of shorter sequence in pairwise comparison) \n\
                        2 = given interval sizes - DNA [default 100 bp] or proteins [default 10 aa] - see -Y",
    )

    parser.add_argument(
        "-Y",
        "--lcs_shading_interval_len",
        type=int,
        default=50,
        help="Length of intervals for LCS shading (only if --lcs_shading_ref=2) [default for nucleotides = 50; default for amino acids = 10]",
    )

    parser.add_argument(
        "-z",
        "--lcs_shading_ori",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Shade subdotplots based on LCS\n\
                            0 = forward [default]\n\
                            1 = reverse, or\n\
                            2 = both strands (forward shading above diagonal, reverse shading on diagonal and below; if using --user_matrix_file, best LCS is used below diagonal)",
    )

    # Custom matrix shading
    parser.add_argument(
        "-u", "--user_matrix_file", help="Matrix file for shading above diagonal."
    )
    parser.add_argument(
        "-U",
        "--user_matrix_print",
        action="store_true",
        default=False,
        help="Display matrix entries above diagonal.",
    )

    # Other options
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'INFO'",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    parser.add_argument("--logfile", default=None, help="Name of log file")

    return parser.parse_args()


###############################
#        Function Call        #
###############################


def main():
    # parse command line arguments
    args = parse_args()

    # Set up logging
    init_logging(loglevel=args.loglevel, logfile=args.logfile)

    # If loglevel is DEBUG, print arguments
    if args.loglevel == "DEBUG":
        logging.debug(args)

    # Log command line arguments
    logging.info("cmd: {0}".format(" ".join(sys.argv)))

    # TODO: Check if files in seq list exist
    # args.infiles
    # TODO: Check if GFF files exist
    # TODO: Check valid kmer length

    # Set up variables
    alphabetic_sorting = args.sort
    convert_wobbles = args.wobble_conversion
    filetype = args.filetype
    gff = args.gff
    gff_color_config_file = args.gff_color_config
    input_user_matrix_file = args.user_matrix_file
    label_size = args.label_size
    lcs_shading_interval_len = args.lcs_shading_interval_len
    lcs_shading_num = args.lcs_shading_num
    lcs_shading_ori = args.lcs_shading_ori
    lcs_shading_ref = args.lcs_shading_ref
    lcs_shading = args.lcs_shading
    length_scaling = args.length_scaling
    max_N_percentage = 49
    mirror_y_axis = args.mirror_y_axis
    modes = args.mode
    multi = args.collage
    ncols = args.n_col
    nrows = args.n_row
    only_vs_first_seq = args.only_vs_first_seq
    plot_size = args.plot_size
    prefix = args.output_prefix
    norevcomp = args.norev
    seq_list = args.infiles
    spacing = args.spacing
    substitution_count = args.substitution_count
    title_clip_pos = "B"  # TODO: Was processed out of title_length in old args
    title_length = args.title_length
    user_matrix_print = args.user_matrix_print
    wordsize = args.wordsize
    line_col_rev = args.line_col_rev
    line_col_for = args.line_col_for

    # TODO: Unaccounted for arguments
    # commandline
    line_width = args.line_width

    # Set True if nucleotide sequence
    if args.type_seq == "nuc":
        type_nuc = True
    elif args.type_seq == "aa":
        type_nuc = False

    # Set x label position
    if args.x_label_pos == "top":
        x_label_pos_top = True
    elif args.x_label_pos == "bottom":
        x_label_pos_top = False

    # This doesn't seem to be used
    # if convert_wobbles:
    #    if type_nuc:
    #        ambiq_res = "N"
    #    else:
    #        ambiq_res = "X"

    # read gff color config file if provided
    # TODO: Figure out if this is necessary
    if args.gff:
        if gff_color_config_file:
            text = "\n%s\n\nReading GFF color configuration file\n%s\n\n=> %s\n" % (
                50 * "=",
                28 * "-",
                gff_color_config_file,
            )
            logging.info(text)
        gff_feat_colors = read_gff_color_config(gff_color_config_file)
    else:
        gff_feat_colors = {}
        if gff_color_config_file:
            logging.warning(
                f"Provide GFF annotation files to use configuration file: {gff_color_config_file}"
            )

    # If color is set to white, reverse complementary matches are skipped
    # TODO: Figure out if "white" actually stops reverse complement calculations
    if norevcomp:  # if norev is set
        line_col_rev = "white"  # reverse matches not calculated
    elif not type_nuc:
        logging.warning("Reverse complement deactivated for proteins.")
        line_col_rev = "white"  # reverse matches not calculated

    # TODO: Log mode types by name
    mode_text = []
    for item in modes:
        mode_text.append(str(item))

    logging.info("Requested plotting modes: %s" % (", ".join(mode_text)))

    # create dotplots
    ##########################################
    # Init empty list for image file names
    list_of_png_names = list()

    # self dotplots
    t1 = time.time()
    if "0" in modes:
        logging.info("Calling selfdotplot")
        list_of_png_names = selfdotplot(
            seq_list,
            wordsize,
            alphabetic_sorting=alphabetic_sorting,
            convert_wobbles=convert_wobbles,
            filetype=args.filetype,
            gff_color_dict=gff_feat_colors,
            gff_files=gff,
            label_size=label_size,
            line_col_rev=args.line_col_rev,
            line_col_for=args.line_col_for,
            max_N_percentage=max_N_percentage,
            mirror_y_axis=mirror_y_axis,
            multi=multi,
            ncols=ncols,
            nrows=nrows,
            plot_size=plot_size,
            prefix=prefix,
            substitution_count=substitution_count,
            title_clip_pos=title_clip_pos,
            title_length=title_length,
            type_nuc=type_nuc,
            line_width=line_width,
        )
        t1 = time_track(t1)
        if list_of_png_names:
            logging.info("-> Image file(s): %s\n" % ", ".join(list_of_png_names))
        else:
            logging.info("No image files were created!\n")
        print(50 * "=")

    # paired dotplots
    if "1" in modes:
        if multi:
            logging.info("Calling pairdotplot with collage")
            list_of_png_names = pairdotplot(
                seq_list,
                wordsize,
                alphabetic_sorting=alphabetic_sorting,
                convert_wobbles=convert_wobbles,
                filetype=filetype,
                label_size=label_size,
                length_scaling=length_scaling,
                line_col_rev=args.line_col_rev,
                line_col_for=args.line_col_for,
                max_N_percentage=max_N_percentage,
                mirror_y_axis=mirror_y_axis,
                multi=multi,
                ncols=ncols,
                nrows=nrows,
                only_vs_first_seq=only_vs_first_seq,
                plot_size=plot_size,
                prefix=prefix,
                substitution_count=substitution_count,
                title_clip_pos=title_clip_pos,
                title_length=title_length,
                type_nuc=type_nuc,
                x_label_pos_top=x_label_pos_top,
                line_width=line_width,
            )
            t1 = time_track(t1)
        else:
            if not length_scaling:
                logging.info(
                    "Pairwise dotplot with individual output files scaled by sequence length automatically."
                )

            logging.info("Calling pairdotplot")
            list_of_png_names = pairdotplot(
                seq_list,
                wordsize,
                alphabetic_sorting=alphabetic_sorting,
                convert_wobbles=convert_wobbles,
                filetype=filetype,
                label_size=label_size,
                length_scaling=True,
                line_col_rev=args.line_col_rev,
                line_col_for=args.line_col_for,
                max_N_percentage=max_N_percentage,
                mirror_y_axis=mirror_y_axis,
                multi=multi,
                ncols=ncols,
                nrows=nrows,
                only_vs_first_seq=only_vs_first_seq,
                plot_size=plot_size,
                prefix=prefix,
                substitution_count=substitution_count,
                title_clip_pos=title_clip_pos,
                title_length=title_length,
                type_nuc=type_nuc,
                line_width=line_width,
            )
            t1 = time_track(t1)
        if list_of_png_names:
            logging.info("-> Image file(s): %s\n" % ", ".join(list_of_png_names))
        else:
            logging.info("No image files were created!\n")

        print(50 * "=")

    # all-against-all dotplot
    if "2" in modes:
        list_of_png_names = polydotplot(
            seq_list,
            wordsize,
            alphabetic_sorting=alphabetic_sorting,
            convert_wobbles=convert_wobbles,
            filetype=filetype,
            gff_color_dict=gff_feat_colors,
            gff_files=gff,
            input_user_matrix_file=input_user_matrix_file,
            label_size=label_size,
            line_col_rev=args.line_col_rev,
            line_col_for=args.line_col_for,
            lcs_shading_interval_len=lcs_shading_interval_len,
            lcs_shading_num=lcs_shading_num,
            lcs_shading_ori=lcs_shading_ori,
            lcs_shading_ref=lcs_shading_ref,
            lcs_shading=lcs_shading,
            max_N_percentage=max_N_percentage,
            mirror_y_axis=mirror_y_axis,
            plot_size=plot_size,
            prefix=prefix,
            representation=args.representation,
            spacing=spacing,
            substitution_count=substitution_count,
            title_clip_pos=title_clip_pos,
            title_length=title_length,
            type_nuc=type_nuc,
            user_matrix_print=user_matrix_print,
            line_width=line_width,
        )
        t1 = time_track(t1)
        if list_of_png_names:
            logging.info("-> Image file(s): %s\n" % ", ".join(list_of_png_names))
        else:
            logging.info("No image files were created!")

        print(50 * "=")

    logging.info("Thank you for using FlexiDot.")


######################
# FlexiDot Execution #
######################
if __name__ == "__main__":
    main()
