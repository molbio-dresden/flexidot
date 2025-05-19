# Old arg handling from flexidot.py
def usage():
    """
    usage and help
    """

    print("""\n\n FLEXIDOT
    -------------------------------------------------------------------

    Version:
    1.06

    Citation: 
    Kathrin M. Seibt, Thomas Schmidt, Tony Heitkam (2018) 
    "FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence investigation"
    Bioinformatics 34 (20), 3575–3577, doi: 10.1093/bioinformatics/bty395

 
    General usage: 
    $ python flexidot.py -a [ARGUMENTS]
    $ python flexidot.py -i <fasta_file_name> [ARGUMENTS]

    
    ARGUMENTS
    -------------------------------------------------------------------

    
    INPUT/OUTPUT OPTIONS... required are [-a] OR [-i]

    -a, --auto_fas                 Imports all fasta files from current directory (*.fasta, *.fas, *.fa, *.fna)
                                    -i is not needed, if -a is activated
                                    [inactive by default]

    -i, --in_file                  Input fasta file (fasta file name or comma-separated file list)
                                   > Provide multiple files: Recall -i or provide comma-separated file names

    -o, --output_file_prefix       File prefix to be added to the generated filenames [default = NONE]

    -c, --collage_output           Multiple dotplots are combined in a collage 
                                    Y or 1 = ON [default]
                                    N or 0 = OFF 

    -m, --m_col                    Number of columns per page [default = 4] (only if --collage_output is ON)

    -n, --n_row                    Number of rows per page [default = 5] (only if --collage_output is ON)

    -f, --filetype                 Output file format
                                    0 = PNG [default]
                                    1 = PDF
                                    2 = SVG 

    -s, --alphabetic_sorting       Sort sequences alphabetically according to titles 
                                    Y or 1 = ON
                                    N or 0 = OFF [default]

 
    CALCULATION PARAMETERS...

    -k, --wordsize                 Wordsize (kmer length) for dotplot comparison [default = 10]

    -p, --plotting_mode            Mode of FlexiDot dotplotting
                                    0 = self [default]
                                    1 = paired
                                    2 = poly (matrix with all-against-all dotplots)
                                   > Run multiple plotting modes: Recall -p or provide comma-separated numbers

    -t, --type_nuc                 Type of residue is nucleotide
                                    Y or 1 = nucleotide [default]
                                    N or 0 = amino acid

    -w, --wobble_conversion        Ambiguity handling for relaxed matching
                                    Y or 1 = ON
                                    N or 0 = OFF [default]

    -S, --substitution_count       Number of substitutions (mismatches) allowed per window for relaxed matching 
                                    [default = 0]

    -r, --rc_option                Find reverse complementary matches (only if type_nuc=y)
                                    Y or 1 = ON  [default]
                                    N or 0 = OFF 

    -O, --only_vs_first_seq        Limit pairwise comparisons to match all sequences to 1st sequence only 
                                   (only if --plotting_mode=1) 
                                    Y or 1 = ON 
                                    N or 0 = OFF [default]
 
    GRAPHIC FORMATTING...

    -A, --line_width               Line width [default = 1]

    -B, --line_col_for             Line color [default = black]

    -C, --line_col_rev             Reverse line color [default = green]

    -D, --x_label_pos              Position of the X-label
                                    Y or 1 = top [default]
                                    N or 0 = bottom

    -E, --label_size               Font size [default = 10]

    -F, --spacing                  Spacing between all-against-all dotplots (only if --plotting_mode=2)
                                    [default = 0.04]

    -L, --length_scaling           Scale plot size for pairwise comparison (only if --plotting_mode=1)  
                                    Y or 1 = Scaling ON (axes scaled according to sequence length)
                                    N or 0 = Scaling OFF (squared plots) [default]

    -M, --mirror_y_axis            Flip y-axis bottom to top (cartesian coordinate system)
                                    Y or 1 = y-axis bottom to top 
                                    N or 0 = y-axis top to bottom [default]

    -P, --plot_size                Plotsize [default = 10]

    -R, --representation           Region of plot to display (only if --plotting_mode=2)
                                    0 = full [default]
                                    1 = upper
                                    2 = lower

    -T, --title_length             Limit title length for dotplot comparisons
                                    [default = 20]
                                   Position of selection can be specified by appending a letter (e.g. -T 20E)
                                    B = beginning [default]
                                    E = end 


    GFF SHADING (for -p/--plotting_mode=0 only)...

    -g, --input_gff_files          GFF3 file used for markup in self-dotplots
                                    (provide multiple files: Recall -g or provide comma-separated file names)

    -G, --gff_color_config_file    Tab-delimited config file for custom gff shading
                                    column 1: feature type
                                    column 2: color
                                    column 3: alpha
                                    column 4: zoom factor (for small regions)


    LCS SHADING OPTIONS (for -p/--plotting_mode=2 only)...

    -x, --lcs_shading              Shade subdotplot based on the length of the longest common substring (LCS)
                                    Y or 1 = ON
                                    N or 0 = OFF [default]

    -X, --lcs_shading_num          Number of shading intervals (hues) for LCS (-x) and user matrix shading (-u) 
                                    [default = 5]

    -y, --lcs_shading_ref          Reference for LCS shading
                                    0 = maximal LCS length [default]
                                    1 = maximally possible length (length of shorter sequence in pairwise comparison)
                                    2 = given interval sizes - DNA [default 100 bp] or proteins [default 10 aa] - see -Y

    -Y, --lcs_shading_interval_len Length of intervals for LCS shading (only if --lcs_shading_ref=2) 
                                    [default for nucleotides = 50; default for amino acids = 10] 

    -z, --lcs_shading_ori          Shade subdotplots according to LCS on
                                    0 = forward [default],
                                    1 = reverse, or
                                    2 = both strands (forward shading above diagonal, reverse shading on diagonal and below;
                                                      if using --input_user_matrix_file, best LCS is used below diagonal)


    CUSTOM USER MATRIX SHADING OPTIONS (for -p/--plotting_mode=2 only)...

    -u, --input_user_matrix_file   Shading above diagonal according to values in matrix file specified by the user
                                    (tab-delimited or comma-separated matrix with sequence name in column 1 and numbers in columns 2-n 
                                     e.g. identity matrix from multiple sequence alignment - strings are ignored)

    -U, --user_matrix_print        Display provided matrix entries in the fields above diagonal of all-against-all dotplot
                                    Y or 1 = ON
                                    N or 0 = OFF [default]


    OTHERS...

    -h, --help                     Help screen

    -v, --verbose                  Verbose




    """)


def check_input(argv, trial_mode=False):
    """
    commandline argument parsing
    """

    global log_txt, aa_bp_unit

    # helpers for argument parsing
    ######################################

    arguments = [
        "-a",
        "--auto_fas",
        "a",
        "auto_fas",
        "-i",
        "--input_fasta",
        "i:",
        "input_fasta=",
        "-o",
        "--output_file_prefix",
        "o:",
        "output_file_prefix=",
        "-c",
        "--collage_output",
        "c:",
        "collage_output=",
        "-m",
        "--m_col",
        "m:",
        "m_col=",
        "-n",
        "--n_row",
        "n:",
        "n_row=",
        "-f",
        "--filetype",
        "f:",
        "filetype=",
        "-t",
        "--type_nuc",
        "t:",
        "type_nuc=",
        "-g",
        "--input_gff_files",
        "g:",
        "input_gff_files",
        "-G",
        "--gff_color_config_file",
        "G:",
        "gff_color_config_file",
        "-k",
        "--wordsize",
        "k:",
        "wordsize=",
        "-p",
        "--plotting_mode",
        "p:",
        "plotting_mode=",
        "-w",
        "--wobble_conversion",
        "w:",
        "wobble_conversion=",
        "-S",
        "--substitution_count",
        "S:",
        "substitution_count=",
        "-r",
        "--rc_option",
        "r:",
        "rc_option=",
        "-O",
        "--only_vs_first_seq",
        "O:",
        "only_vs_first_seq=",
        "-s",
        "--alphabetic_sorting",
        "s:",
        "alphabetic_sorting=",
        "-x",
        "--lcs_shading",
        "x:",
        "lcs_shading=",
        "-X",
        "--lcs_shading_num",
        "X:",
        "lcs_shading_num=",
        "-y",
        "--lcs_shading_ref",
        "y:",
        "lcs_shading_ref=",
        "-Y",
        "--lcs_shading_interval_len",
        "Y:",
        "lcs_shading_interval_len=",
        "-z",
        "--lcs_shading_ori",
        "z:",
        "lcs_shading_ori=",
        "-u",
        "--input_user_matrix_file",
        "u:",
        "input_user_matrix_file=",
        "-U",
        "--user_matrix_print",
        "U:",
        "user_matrix_print=",
        "-P",
        "--plot_size",
        "P:",
        "plot_size=",
        "-A",
        "--line_width",
        "A:",
        "line_width=",
        "-B",
        "--line_col_for",
        "B:",
        "line_col_for=",
        "-C",
        "--line_col_rev",
        "C:",
        "line_col_rev=",
        "-D",
        "--x_label_pos",
        "D:",
        "x_label_pos=",
        "-E",
        "--label_size",
        "E:",
        "label_size=",
        "-F",
        "--spacing",
        "F:",
        "spacing=",
        "-L",
        "--length_scaling",
        "L:",
        "length_scaling=",
        "-M",
        "--mirror_y_axis",
        "M:",
        "mirror_y_axis=",
        "-R",
        "--representation",
        "R:",
        "representation=",
        "-T",
        "--title_length",
        "T:",
        "title_length=",
        "-h",
        "--help",
        "h",
        "help",
        "-v",
        "--verbose",
        "v",
        "verbose",
    ]

    arguments_sysargv = tuple(arguments[0::4] + arguments[1::4])
    arguments_opts = "".join(arguments[2::4])
    arguments_args = arguments[3::4]

    # setting defaults
    ######################################

    auto_fas = False  # 0
    input_fasta = []
    output_file_prefix = None
    collage_output = True  # 1
    m_col = 4
    n_row = 5
    filetype = 0
    type_nuc = True
    input_gff_files = []
    gff_color_config_file = ""

    wordsize = 10
    plotting_modes = [0]
    wobble_conversion = False  # 0
    substitution_count = 0
    rc_option = True  # 1
    alphabetic_sorting = False  # 0
    only_vs_first_seq = False  # 0

    lcs_shading = False  # 0
    lcs_shading_num = 4
    lcs_shading_ref = 0
    lcs_shading_interval_len = (
        50  # interval default changes to "10" for amino acids [type_nuc = n]
    )
    lcs_shading_ori = 0

    input_user_matrix_file = ""
    user_matrix_print = False

    plot_size = 10
    line_width = 1
    line_col_for = "black"
    line_col_rev = "#009243"
    x_label_pos = True  # 0
    label_size = 10
    spacing = 0.04
    length_scaling = False  # 0
    title_length = 20  # float("Inf")
    title_clip_pos = "B"  # B (begin), E (end)
    max_N_percentage = 49  # fixed value, no user input
    mirror_y_axis = False
    representation = 0

    aa_bp_unit = "bp"

    verbose = False  # 0

    filetype_dict = {0: "png", 1: "pdf", 2: "svg"}
    lcs_shading_ref_dict = {
        0: "maximal LCS length",
        1: "maximally possible length",
        2: "given interval sizes",
    }
    plotting_mode_dict = {0: "self", 1: "paired", 2: "all-against-all"}
    lcs_shading_ori_dict = {0: "forward", 1: "reverse complement", 2: "both"}
    representation_dict = {0: "full", 1: "upper", 2: "lower"}

    # return default parameters for testing purposes
    if trial_mode:
        print("ATTENTION: YOU ARE IN THE TRIAL MODE!!!\n\n")

        commandline = "trial_mode\n"

        parameters = [
            commandline,
            auto_fas,
            input_fasta,
            output_file_prefix,
            collage_output,
            m_col,
            n_row,
            filetype_dict[filetype],
            type_nuc,
            input_gff_files,
            gff_color_config_file,
            wordsize,
            plotting_modes,
            wobble_conversion,
            substitution_count,
            rc_option,
            alphabetic_sorting,
            only_vs_first_seq,
            lcs_shading,
            lcs_shading_num,
            lcs_shading_ref,
            lcs_shading_interval_len,
            lcs_shading_ori,
            input_user_matrix_file,
            user_matrix_print,
            plot_size,
            line_width,
            line_col_for,
            line_col_rev,
            x_label_pos,
            label_size,
            spacing,
            length_scaling,
            title_length,
            title_clip_pos,
            max_N_percentage,
            mirror_y_axis,
            representation,
            verbose,
        ]
        return parameters

    # read arguments
    ######################################

    commandline = ""
    for arg in sys.argv:
        commandline += arg + " "

    log_txt = "\n...reading input arguments..."
    print(log_txt)

    if len(sys.argv) < 2:
        print("\nERROR: More arguments are needed. Exit...")
        log_txt += "\nERROR: More arguments are needed. Exit..."
        usage()
        sys.exit()

    elif sys.argv[1] not in arguments_sysargv:
        print(
            "\nINPUT ERROR: Input argument %s unknown. Please check the help screen."
            % sys.argv[1]
        )
        log_txt += (
            "\nINPUT ERROR: Input argument %s unknown. Please check the help screen."
            % sys.argv[1]
        )
        # usage()
        sys.exit()

    try:
        opts, args = getopt.getopt(sys.argv[1:], arguments_opts, arguments_args)

    except getopt.GetoptError:
        print(
            "\nINPUT ERROR (getopt): Input argument %s unknown. Please check the help screen."
            % sys.argv[1:]
        )
        log_txt += (
            "\nINPUT ERROR (getopt): Input argument %s unknown. Please check the help screen."
            % sys.argv[1:]
        )
        # usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("...fetch help screen")
            log_txt += "\n...fetch help screen"
            usage(), sys.exit()

        if opt in ("-v", "--verbose"):
            print("...verbose output")
            log_txt += "\n...verbose output"
            verbose = True

        elif opt in ("-i", "--input_fasta"):
            if "," in arg:
                arg_list = arg.split(",")
                for temp_file in arg_list:
                    if not os.path.exists(str(temp_file)):
                        message = "\nERROR: fasta_file '%s' was not found!" % str(
                            temp_file
                        )
                        sys.exit(message)
                    else:
                        input_fasta.append(str(temp_file))
                        print("fasta file #%i: %s" % (len(input_fasta), str(temp_file)))
                        log_txt += "\nfasta file #%i: %s" % (
                            len(input_fasta),
                            str(temp_file),
                        )
            else:
                if not os.path.exists(str(arg)):
                    message = "\nERROR: fasta_file '%s' was not found!" % str(arg)
                    log_txt += message
                    sys.exit(message)
                else:
                    input_fasta.append(str(arg))
                    print("fasta file #%i: %s" % (len(input_fasta), str(arg)))
                    log_txt += "\nfasta file #%i: %s" % (len(input_fasta), str(arg))

        elif opt in ("-a", "--auto_fas"):
            auto_fas = True

        # multiple gff files: reads them into a list
        elif opt in ("-g", "--input_gff_files"):
            # append gff file only if existing
            if "," in arg:
                arg_list = arg.split(",")
                for temp_file in arg_list:
                    if not os.path.exists(str(temp_file)):
                        message = "\nERROR: gff_file '%s' was not found!" % str(
                            temp_file
                        )
                        print(message)
                        log_txt += message
                        print("  -->Running FlexiDot without this gff file!")
                        log_txt += "\n  -->Running FlexiDot without this gff file!"
                    else:
                        print(
                            "GFF file #%i: %s" % (len(input_gff_files), str(temp_file))
                        )
                        log_txt += "\nGFF file #%i: %s" % (
                            len(input_gff_files),
                            str(temp_file),
                        )
                        input_gff_files.append(str(temp_file))
            else:
                if not os.path.exists(str(arg)):
                    message = "\nERROR: gff_file '%s' was not found!" % str(arg)
                    print(message)
                    log_txt += message
                    print("  -->Running FlexiDot without this gff file!")
                    log_txt += "\n  -->Running FlexiDot without this gff file!"
                else:
                    input_gff_files.append(str(arg))
                    print("GFF file #%i: %s" % (len(input_gff_files), str(arg)))
                    log_txt += "\nGFF file #%i: %s" % (len(input_gff_files), str(arg))

        elif opt in ("-G", "--gff_color_config_file"):
            if not os.path.exists(str(arg)):
                message = "\nERROR: gff_color_config_file '%s' was not found!" % str(
                    arg
                )
                print(
                    message
                    + "\n  -->Running FlexiDot with default gff coloring specification!"
                )
                log_txt += (
                    message
                    + "\n  -->Running FlexiDot with default gff coloring specification!"
                )
            else:
                gff_color_config_file = str(arg)

        elif opt in ("-u", "--input_user_matrix_file"):
            if not os.path.exists(str(arg)):
                message = "\nERROR: input_user_matrix_file '%s' was not found!" % str(
                    arg
                )
                print(
                    message
                    + "\n  -->Running FlexiDot without input_user_matrix_file %s!" % arg
                )
                log_txt += (
                    message + "\n  -->Running FlexiDot withdefault matrix shading file!"
                )
            else:
                input_user_matrix_file = str(arg)

        elif opt in ("-U", "--user_matrix_print"):
            user_matrix_print = check_bools(str(arg), default=user_matrix_print)

        elif opt in ("-o", "--output_file_prefix"):
            output_file_prefix = arg

        elif opt in ("-c", "--collage_output"):
            collage_output = check_bools(str(arg), default=collage_output)

        elif opt in ("-m", "--m_col"):
            try:
                m_col = int(arg)
            except:
                print("m_col - invalid argument - using default value")
                log_txt += "\nm_col - invalid argument - using default value"

        elif opt in ("-n", "--n_row"):
            try:
                n_row = int(arg)
            except:
                print("n_row - invalid argument - using default value")
                log_txt += "\nn_row - invalid argument - using default value"

        elif opt in ("-f", "--filetype"):
            if 0 <= int(arg) <= 2:
                filetype = int(arg)
            else:
                print(
                    "\nERROR: Please provide valid filetype argument. %s is out of range. It will be set to -f 0 [default]."
                    % (filetype)
                )
                log_txt += (
                    "\nERROR: Please provide valid filetype argument. %s is out of range. It will be set to -f 0 [default]."
                    % (filetype)
                )

        elif opt in ("-t", "--type_nuc"):
            type_nuc = check_bools(str(arg), default=type_nuc)

            if type_nuc == False:
                # interval default changed for amino acids
                lcs_shading_interval_len = 10
                aa_bp_unit = "aa"

        elif opt in ("-k", "--wordsize"):
            try:
                wordsize = int(arg)
            except:
                print("wordsize - invalid argument - using default value")
                log_txt += "\nwordsize - invalid argument - using default value"

        elif opt in ("-p", "--plotting_mode"):
            if "," in arg:
                temp_modes = arg.split(",")
                for item in temp_modes:
                    if item in ["0", "1", "2"]:
                        plotting_modes.append(int(item))
            elif arg in ["0", "1", "2"]:
                plotting_modes = [int(arg)]
            else:
                print(
                    "Please provide valid plotting_modes argument - e.g. 1 or 0,1,2 - using default [0]"
                )
                log_txt += "\nPlease provide valid plotting_modes argument - e.g. 1 or 0,1,2 - using default [0]"

        elif opt in ("-w", "--wobble_conversion"):
            wobble_conversion = check_bools(str(arg), default=wobble_conversion)

        elif opt in ("-S", "--substitution_count"):
            try:
                substitution_count = int(arg)
            except:
                print("substitution_count - invalid argument - using default value")
                log_txt += (
                    "\nsubstitution_count - invalid argument - using default value"
                )

        elif opt in ("-r", "--rc_option"):
            rc_option = check_bools(str(arg), default=rc_option)

        elif opt in ("-s", "--alphabetic_sorting"):
            alphabetic_sorting = check_bools(str(arg), default=alphabetic_sorting)

        elif opt in ("-O", "--only_vs_first_seq"):
            only_vs_first_seq = check_bools(str(arg), default=only_vs_first_seq)

        elif opt in ("-x", "--lcs_shading"):
            lcs_shading = check_bools(str(arg), default=lcs_shading)

        elif opt in ("-X", "--lcs_shading_num"):
            try:
                lcs_shading_num = int(arg) - 1
            except:
                print("lcs_shading_num - invalid argument - using default value")
                log_txt += "\nlcs_shading_num - invalid argument - using default value"

        elif opt in ("-y", "--lcs_shading_ref"):
            try:
                if 0 <= int(arg) <= 2:
                    lcs_shading_ref = int(arg)
                else:
                    print(
                        "\nERROR: lcs_shading_ref %s out of range. It will be set to -y 0 [default]."
                        % (lcs_shading_ref)
                    )
                    log_txt += (
                        "\nERROR: lcs_shading_ref %s out of range. It will be set to -y 0 [default]."
                        % (lcs_shading_ref)
                    )
            except:
                print("lcs_shading_ref - invalid argument - using default value")
                log_txt += "\nlcs_shading_ref - invalid argument - using default value"

        elif opt in ("-Y", "--lcs_shading_interval_len"):
            try:
                lcs_shading_interval_len = int(arg)
            except:
                print(
                    "lcs_shading_interval_len - invalid argument - using default value"
                )
                log_txt += "\nlcs_shading_interval_len - invalid argument - using default value"

        elif opt in ("-z", "--lcs_shading_ori"):
            if 0 <= int(arg) <= 2:
                lcs_shading_ori = int(arg)
            else:
                print(
                    "\nERROR: Please provide valid lcs_shading_ori argument. %s is out of range. It will be set to -z 0 [default]."
                    % (lcs_shading_ori)
                )
                log_txt += (
                    "\nERROR: Please provide valid lcs_shading_ori argument. %s is out of range. It will be set to -z 0 [default]."
                    % (lcs_shading_ori)
                )

        elif opt in ("-P", "--plot_size"):
            try:
                plot_size = float(arg)
            except:
                print("plot_size - invalid argument - using default value")
                log_txt += "\nplot_size - invalid argument - using default value"

        elif opt in ("-A", "--line_width"):
            try:
                line_width = float(arg)
            except:
                print("line_width - invalid argument - using default value")
                log_txt += "\nline_width - invalid argument - using default value"

        elif opt in ("-B", "--line_col_for"):
            if mcolors.is_color_like(arg):
                line_col_for = arg
            else:
                print("line_col_for - invalid argument - using default value")
                log_txt += "\nline_col_for - invalid argument - using default value"

        elif opt in ("-C", "--line_col_rev"):
            if mcolors.is_color_like(arg):
                line_col_rev = arg
            else:
                print("line_col_rev - invalid argument - using default value")
                log_txt += "\nline_col_rev - invalid argument - using default value"

        elif opt in ("-D", "--x_label_pos"):
            x_label_pos = check_bools(str(arg), default=x_label_pos)

        elif opt in ("-E", "--label_size"):
            try:
                label_size = float(arg)
            except:
                print("label_size - invalid argument - using default value")
                log_txt += "\nlabel_size - invalid argument - using default value"

        elif opt in ("-F", "--spacing"):
            try:
                spacing = float(arg)
            except:
                print("spacing - invalid argument - using default value")
                log_txt += "\nspacing - invalid argument - using default value"

        elif opt in ("-L", "--length_scaling"):
            length_scaling = check_bools(str(arg), default=length_scaling)

        elif opt in ("-M", "--mirror_y_axis"):
            mirror_y_axis = check_bools(str(arg), default=mirror_y_axis)

        elif opt in ("-R", "--representation"):
            if 0 <= int(arg) <= 2:
                representation = int(arg)
            else:
                print(
                    "\nERROR: Please provide valid representation argument. %s is out of range. It will be set to -R 0 [default]."
                    % (representation)
                )
                log_txt += (
                    "\nERROR: Please provide valid representation argument. %s is out of range. It will be set to -R 0 [default]."
                    % (representation)
                )

        elif opt in ("-T", "--title_length"):
            try:
                title_length = int(arg)
            except:
                try:
                    title_length = int(str(arg)[:-1])
                    if arg[-1].upper() in ["B", "E"]:  # B (beginning), E (end)
                        title_clip_pos = arg[-1].upper()
                    else:
                        print(
                            "title_length position information invalid - using default value"
                        )
                        log_txt += "\ntitle_length position information invalid - using default value"
                except:
                    print("title_length - invalid argument - using default value")
                    log_txt += "\ntitle_length - invalid argument - using default value"

    # start logging file
    logprint(commandline, start=True, printing=False, prefix=output_file_prefix)
    logprint(log_txt, start=False, printing=False)
    
    


    # print chosen arguments
    ######################################

    text = "\n%s\n" % (70 * "-")
    text += "\n" + "INPUT/OUTPUT OPTIONS...\n"
    text += (
        "\n"
        + "Input fasta file:                                  "
        + ", ".join(input_fasta)
    )
    text += "\n" + "Automatic fasta collection from current directory: " + str(auto_fas)
    text += (
        "\n"
        + "Collage output:                                    "
        + str(collage_output)
    )
    text += "\n" + "Number of columns per page:                        " + str(m_col)
    text += "\n" + "Number of rows per page:                           " + str(n_row)
    text += (
        "\n"
        + "File format:                                       "
        + filetype_dict[filetype]
    )
    text += "\n" + "Residue type is nucleotide:                        " + str(type_nuc)

    text += "\n" + "\n\nCALCULATION PARAMETERS...\n"
    text += "\n" + "Wordsize:                                          " + str(wordsize)
    text += (
        "\n"
        + "Sustitution count:                                 "
        + str(substitution_count)
    )
    text += (
        "\n"
        + "Plotting mode:                                     "
        + str(plotting_modes).replace("[", "").replace("]", "")
        + "\n"
        + 51 * " "
    )
    for item in plotting_modes:
        text += plotting_mode_dict[item] + " "
    text += (
        "\n"
        + "Ambiguity handling:                                "
        + str(wobble_conversion)
    )
    text += (
        "\n" + "Reverse complement scanning:                       " + str(rc_option)
    )
    text += (
        "\n"
        + "Alphabetic sorting:                                "
        + str(alphabetic_sorting)
    )

    if 1 in plotting_modes:
        text += (
            "\n"
            + "Only matching sequences to first entry:            "
            + str(only_vs_first_seq)
        )

    if 0 in plotting_modes and input_gff_files != []:
        text += (
            "\n"
            + "Input gff files:                                   "
            + ", ".join(input_gff_files)
        )
        if gff_color_config_file != "":
            text += (
                "\n"
                + "GFF color config file:                             "
                + gff_color_config_file
            )
    text += (
        "\n"
        + "Prefix for output files:                           "
        + str(output_file_prefix)
    )

    if 2 in plotting_modes:
        text += (
            "\n" + "\n\nLCS SHADING OPTIONS (plotting_mode 'all-against-all' only)...\n"
        )
        text += (
            "\n"
            + "LCS shading:                                       "
            + str(lcs_shading)
        )
        text += (
            "\n"
            + "LCS shading interval number:                       "
            + str(lcs_shading_num + 1)
        )
        text += (
            "\n"
            + "LCS shading reference:                             "
            + lcs_shading_ref_dict[lcs_shading_ref]
        )
        if lcs_shading_ref == 2:
            text += (
                "\n"
                + "LCS shading interval size [%s]:                    " % (aa_bp_unit)
                + str(lcs_shading_interval_len)
            )
        text += (
            "\n"
            + "LCS shading orientation:                           "
            + lcs_shading_ori_dict[lcs_shading_ori]
        )
        if input_user_matrix_file != "":
            text += (
                "\n"
                + "Custom user shading matrix file:                   "
                + input_user_matrix_file
            )
            text += (
                "\n"
                + "Print user matrix values (instead of dotplot):     "
                + str(user_matrix_print)
            )
        text += (
            "\n"
            + "Displayed plot region:                             "
            + representation_dict[representation]
        )

    text += "\n" + "\n\nGRAPHIC FORMATTING...\n"
    text += (
        "\n" + "Plot size:                                         " + str(plot_size)
    )
    text += (
        "\n" + "Line width:                                        " + str(line_width)
    )
    text += "\n" + "Line color:                                        " + line_col_for
    text += "\n" + "Reverse line color:                                " + line_col_rev
    text += (
        "\n" + "X label position:                                  " + str(x_label_pos)
    )
    text += (
        "\n" + "Label size:                                        " + str(label_size)
    )
    text += "\n" + "Spacing:                                           " + str(spacing)
    if mirror_y_axis:
        text += (
            "\n"
            + "Y-axis mirrored (bottom to top)                   "
            + str(mirror_y_axis)
        )
    if title_clip_pos == "E":
        text += (
            "\n"
            + "Title length (limit number of characters):         "
            + "last"
            + str(title_length)
            + "characters"
        )
    else:
        text += (
            "\n"
            + "Title length (limit number of characters):         "
            + "first"
            + str(title_length)
            + "characters"
        )
    text += (
        "\n"
        + "Length scaling:                                    "
        + str(length_scaling)
    )
    text += "\n%s\n" % (70 * "-")
    logprint(text)

    # collect settings
    parameters = [
        commandline,
        auto_fas,
        input_fasta,
        output_file_prefix,
        collage_output,
        m_col,
        n_row,
        filetype_dict[filetype],
        type_nuc,
        input_gff_files,
        gff_color_config_file,
        wordsize,
        plotting_modes,
        wobble_conversion,
        substitution_count,
        rc_option,
        alphabetic_sorting,
        only_vs_first_seq,
        lcs_shading,
        lcs_shading_num,
        lcs_shading_ref,
        lcs_shading_interval_len,
        lcs_shading_ori,
        input_user_matrix_file,
        user_matrix_print,
        plot_size,
        line_width,
        line_col_for,
        line_col_rev,
        x_label_pos,
        label_size,
        spacing,
        length_scaling,
        title_length,
        title_clip_pos,
        max_N_percentage,
        mirror_y_axis,
        representation,
        verbose,
    ]

    return parameters