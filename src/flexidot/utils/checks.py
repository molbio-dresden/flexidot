import logging
import os
import sys

from flexidot._version import __version__


def check_kmer_length(kmer_length):
    """
    Logs a warning if the kmer length is less than 10.

    Parameters:
    kmer_length (int): The kmer length to check.
    """
    if kmer_length < 10:
        logging.warning(
            'Kmer length is less than 10. This may result in less accurate dotplots.'
        )


def check_input(filename):
    """
    Check if input files exist
    """
    # Check if the input file exists
    if not os.path.isfile(filename):
        logging.error('Input file does not exist. Quitting.')
        raise FileNotFoundError(f"Input sequence file '{filename}' does not exist.")
    sys.stderr.write(f'\033[92m  Input file found: {filename}\033[0m\n')


def check_output_dir(outdir):
    """
    Check if the output directory exists
    """
    # Check if the output directory exists
    if not os.path.isdir(outdir):
        logging.error('Output directory does not exist. Quitting.')
        raise FileNotFoundError(f"Output directory '{outdir}' does not exist.")
    sys.stderr.write(f'\033[92m  Output directory found: {outdir}\033[0m\n\n')


def print_summary(args):
    """
    Print a summary of the selected options to stderr.
    Verifies the existence of any provided filenames.
    """
    # Define color for headings
    heading_color = '\033[1;34m'
    reset_color = '\033[0m'

    # Log command line arguments
    sys.stderr.write(f'\n{heading_color}CMD line input{reset_color}\n')
    sys.stderr.write('  ${0}\n\n'.format(' '.join(sys.argv)))

    # Input/output arguments
    sys.stderr.write(f'{heading_color}Input/Output Arguments:{reset_color}\n')
    sys.stderr.write(f'  Input files [-i]: {args.infiles}\n')
    for infile in args.infiles:
        check_input(infile)
    sys.stderr.write(f'  Output prefix [-o]: {args.output_prefix}\n')
    sys.stderr.write(f'  Output directory [--outdir]: {args.outdir}\n')
    check_output_dir(args.outdir)

    # Collage arguments
    sys.stderr.write(f'{heading_color}Collage Arguments:{reset_color}\n')
    sys.stderr.write(f'  Collage [-c]: {args.collage}\n')
    sys.stderr.write(f'  Number of columns per page [--n_col]: {args.n_col}\n')
    sys.stderr.write(f'  Number of rows per page [--n_row]: {args.n_row}\n\n')

    # File type
    sys.stderr.write(f'{heading_color}File Type:{reset_color}\n')
    sys.stderr.write(f'  Output file format [-f]: {args.filetype}\n\n')

    # Sorting
    sys.stderr.write(f'{heading_color}Sorting:{reset_color}\n')
    sys.stderr.write(f'  Sort sequences alphabetically [-s]: {args.sort}\n\n')

    # Calculation parameters
    sys.stderr.write(f'{heading_color}Calculation Parameters:{reset_color}\n')
    sys.stderr.write(f'  Wordsize (kmer length) [-k]: {args.wordsize}\n')
    mode_descriptions = {
        '0': 'self [default]',
        '1': 'paired',
        '2': 'poly (matrix with all-against-all dotplots)',
    }
    modes = [mode_descriptions[mode] for mode in args.mode]
    sys.stderr.write(f'  Modes [-m]: {", ".join(modes)}\n')
    sys.stderr.write(f'  Sequence type [-t]: {args.type_seq}\n')
    sys.stderr.write(f'  Wobble conversion [-w]: {args.wobble_conversion}\n')
    sys.stderr.write(f"  Maximum 'N' percentage [-n_max]: {args.max_n}%\n")
    sys.stderr.write(f'  Substitution count [-S]: {args.substitution_count}\n')
    sys.stderr.write(f'  No reverse complementary matches [-r]: {args.norev}\n')
    sys.stderr.write(f'  Only vs first sequence [-O]: {args.only_vs_first_seq}\n\n')

    # Graphic formatting
    sys.stderr.write(f'{heading_color}Graphic Formatting:{reset_color}\n')
    sys.stderr.write(f'  Line width [-A]: {args.line_width}\n')
    sys.stderr.write(f'  Line color (forward) [-B]: {args.line_col_for}\n')
    sys.stderr.write(f'  Line color (reverse) [-C]: {args.line_col_rev}\n')
    sys.stderr.write(f'  X-label position [-D]: {args.x_label_pos}\n')
    sys.stderr.write(f'  Font size [-E]: {args.label_size}\n')
    sys.stderr.write(f'  Spacing between dotplots [-F]: {args.spacing}\n')
    sys.stderr.write(f'  Length scaling [-L]: {args.length_scaling}\n')
    sys.stderr.write(f'  Mirror y-axis [-M]: {args.mirror_y_axis}\n')
    sys.stderr.write(f'  Plot size [-P]: {args.plot_size}\n')
    representation_descriptions = {0: 'full [default]', 1: 'upper', 2: 'lower'}
    sys.stderr.write(
        f'  Representation [-R]: {representation_descriptions[args.representation]}\n'
    )
    sys.stderr.write(f'  Title length [-T]: {args.title_length}\n\n')

    # GFF shading
    sys.stderr.write(f'{heading_color}GFF Shading:{reset_color}\n')
    sys.stderr.write(f'  GFF files [-g]: {args.gff}\n')
    if args.gff:
        for gff_file in args.gff:
            check_input(gff_file)
    if args.gff_color_config:
        sys.stderr.write(f'  GFF color config [-G]: {args.gff_color_config}\n')
        check_input(args.gff_color_config)
        sys.stderr.write('\n')

    # Longest Common Subsequence (LCS) shading
    sys.stderr.write(
        f'\n{heading_color}Longest Common Subsequence (LCS) Shading:{reset_color}\n'
    )
    sys.stderr.write(f'  LCS shading [-x]: {args.lcs_shading}\n')
    sys.stderr.write(f'  LCS shading intervals [-X]: {args.lcs_shading_num}\n')
    lcs_shading_ref_descriptions = {
        0: 'maximal LCS length [default]',
        1: 'maximally possible length',
        2: 'given interval sizes',
    }
    sys.stderr.write(
        f'  LCS shading reference [-y]: {lcs_shading_ref_descriptions[args.lcs_shading_ref]}\n'
    )
    sys.stderr.write(
        f'  LCS shading interval length [-Y]: {args.lcs_shading_interval_len}\n'
    )
    lcs_shading_ori_descriptions = {
        0: 'forward [default]',
        1: 'reverse',
        2: 'both strands',
    }
    sys.stderr.write(
        f'  LCS shading orientation [-z]: {lcs_shading_ori_descriptions[args.lcs_shading_ori]}\n\n'
    )

    # Custom matrix shading
    sys.stderr.write(f'{heading_color}Custom Matrix Shading:{reset_color}\n')
    sys.stderr.write(f'  User matrix file [-u]: {args.user_matrix_file}\n')
    if args.user_matrix_file:
        check_input(args.user_matrix_file)
    sys.stderr.write(f'  User matrix print [-U]: {args.user_matrix_print}\n\n')

    # Other options
    sys.stderr.write(f'{heading_color}Other Options:{reset_color}\n')
    sys.stderr.write(f'  Logging level [--loglevel]: {args.loglevel}\n')
    sys.stderr.write(f'  Log file [--logfile]: {args.logfile}\n')
    sys.stderr.write(f'  Version [-v]: {__version__}\n\n')
