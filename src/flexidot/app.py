import logging
import sys

from matplotlib import rc as mplrc
from matplotlib import rcParams
import pylab as P

from flexidot._version import __version__
from flexidot.plotting import selfdotplot, pairdotplot, polydotplot
from flexidot.utils.args import parse_args
from flexidot.utils.checks import print_summary, check_kmer_length
from flexidot.utils.file_handling import read_gff_color_config
from flexidot.utils.logs import init_logging

# Matplotlib settings

# Switch to non-interactive backend to avoid _tkinter.TclError on CentOs 7 servers see Github Issue #5
P.switch_backend("agg") 

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

###############################
#        Function Call        #
###############################

def main():
    # Parse command line arguments
    args = parse_args()

    # Set up logging
    init_logging(loglevel=args.loglevel, logfile=args.logfile)
    
    # Print summary of arguments
    print_summary(args)

    # Log version and command line arguments if debug is enabled
    logging.debug("FlexiDot version: %s" % __version__)
    logging.debug("  ${0}\n\n".format(" ".join(sys.argv)))
    logging.debug("Command line arguments: %s" % args)
    
    # Check valid kmer length
    check_kmer_length(args.wordsize)

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
    max_N_percentage = args.max_n
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

    # Read gff color config file if provided
    if args.gff:
        if gff_color_config_file:
            logging.info(f"Reading GFF color configuration file: {gff_color_config_file}")
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
    mode_names = {"0": "self", "1": "paired", "2": "poly"}
    for item in modes:
        mode_text.append("\n" + str(item) + ": " + mode_names[item])

    logging.info("\nRequested plotting modes: %s" % (", ".join(mode_text)))

    # Create dotplots
    ##########################################
    
    # Init empty list for image file names
    list_of_png_names = list()

    # self dotplots
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
            line_col_rev=line_col_rev,
            line_col_for=line_col_for,
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
        #t1 = time_track(t1)
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
                line_col_rev=line_col_rev,
                line_col_for=line_col_for,
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
            #t1 = time_track(t1)
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
                line_col_rev=line_col_rev,
                line_col_for=line_col_for,
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
            #t1 = time_track(t1)
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
            line_col_rev=line_col_rev,
            line_col_for=line_col_for,
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
        #t1 = time_track(t1)
        if list_of_png_names:
            logging.info("\n-> Image file(s): %s\n" % ",\n".join(list_of_png_names))
        else:
            logging.info("No image files were created!")

    logging.info("Finished. Thank you for using FlexiDot.")


######################
# FlexiDot Execution #
######################
if __name__ == "__main__":
    main()
