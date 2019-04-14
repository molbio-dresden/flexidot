#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

"""
FlexiDot Version 1.06

FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence investigation

Kathrin M. Seibt, Thomas Schmidt and Tony Heitkam 
Institute of Botany, TU Dresden, Dresden, 01277, Germany

Bioinformatics (2018) Vol. 34 (20), 3575–3577, doi 10.1093/bioinformatics/bty395
"""


###############################
#        Requirements         #
###############################

# import system modules
import os, glob
import time, datetime
import sys
import shutil, getopt
import unicodedata

def module_install_command(module_name, upgrade=False):
    """
    create installation commands for Python modules and print information 
    """
    if upgrade:
        load_command = "python -m pip install --upgrade %s" % module_name
    else:
        load_command = "python -m pip install %s" % module_name

    try:
        logprint("Installing Python module: %s\n\t%s\n" % (module_name, load_command))
    except:
        print "Installing Python module: %s\n\t%s\n" % (module_name, load_command)

    return load_command

def load_modules():
    """
    load Python modules, if possible - otherwise try to install them
    """
    # make module names global
    global cllct, gridspec, patches, rcParams, mplrc, P, Color, SeqIO, np, mcolors, rgb2hex, regex

    # matplotlib
    try:
        import matplotlib.collections as cllct
    except:
        command = module_install_command("matplotlib", upgrade=True)
        try:
            os.system(command) 
            print "\n"
            import matplotlib.collections as cllct
        except:
            print "Please install module matplotlib manually"
    import matplotlib.colors as mcolors
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as patches
    import pylab as P
    P.switch_backend('agg') # bugfix for _tkinter.TclError on CentOs 7 servers, see Github Issue #5

    # specify matplotlib font settings
    from matplotlib import rc as mplrc
    mplrc('pdf', fonttype=42, compression=0)
    from matplotlib import rcParams
    rcParams['font.family']     = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Verdana', 'Tahoma' , 'DejaVu Sans', 'Droid Sans Mono', 'Sans', 'Liberation', 'Ubuntu', 'Arial', ]

    # colour for color gradient palette
    try:
        from colour import Color
    except:
        command = module_install_command("colour")
        try:
            os.system(command) 
            print "\n"
            from colour import Color
        except:
            print "Please install module colour manually"

    # color converter
    try:
        from colormap import rgb2hex
    except:
        command = module_install_command("colormap")
        # additional module easydev.tools required by colormap
        command2 = module_install_command("easydev")
        try:
            os.system(command) 
            os.system(command2) 
            print "\n"
            from colormap import rgb2hex
        except:
            print "Please install module colormap manually"

    # biopython
    try:
        from Bio import SeqIO
    except:
        command = module_install_command("biopython")
        try:
            os.system(command) 
            print "\n"
            from Bio import SeqIO
        except:
            print "Please install module biopython manually"

    # numpy
    try:
        import numpy as np
    except:
        command = module_install_command("numpy")
        try:
            os.system(command) 
            print "\n"
            import numpy as np
        except:
            print "Please install module numpy manually"

    # regex for pattern matching
    try:
        import regex
    except:
        command = module_install_command("regex")
        try:
            os.system(command) 
            print "\n"
            import regex
        except:
            print "Please install module regex manually"



###############################
#        Usage & Input        #
###############################

def usage():
    """
    usage and help
    """

    print """\n\n FLEXIDOT
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




    """

def check_input(argv, trial_mode=False):
    """
    commandline argument parsing
    """

    global log_txt, aa_bp_unit

    # helpers for argument parsing
    ######################################

    arguments = ["-a", "--auto_fas",                 "a", "auto_fas",
                 "-i", "--input_fasta",              "i:", "input_fasta=",
                 "-o", "--output_file_prefix",       "o:", "output_file_prefix=",
                 "-c", "--collage_output",           "c:", "collage_output=",
                 "-m", "--m_col",                    "m:", "m_col=",
                 "-n", "--n_row",                    "n:", "n_row=",
                 "-f", "--filetype",                 "f:", "filetype=",
                 "-t", "--type_nuc",                 "t:", "type_nuc=",
                 "-g", "--input_gff_files",          "g:", "input_gff_files",
                 "-G", "--gff_color_config_file",    "G:", "gff_color_config_file",
                 "-k", "--wordsize",                 "k:", "wordsize=",
                 "-p", "--plotting_mode",            "p:", "plotting_mode=",
                 "-w", "--wobble_conversion",        "w:", "wobble_conversion=",
                 "-S", "--substitution_count",       "S:", "substitution_count=",
                 "-r", "--rc_option",                "r:", "rc_option=",
                 "-O", "--only_vs_first_seq",        "O:", "only_vs_first_seq=",
                 "-s", "--alphabetic_sorting",       "s:", "alphabetic_sorting=",
                 "-x", "--lcs_shading",              "x:", "lcs_shading=",
                 "-X", "--lcs_shading_num",          "X:", "lcs_shading_num=",
                 "-y", "--lcs_shading_ref",          "y:", "lcs_shading_ref=",
                 "-Y", "--lcs_shading_interval_len", "Y:", "lcs_shading_interval_len=",
                 "-z", "--lcs_shading_ori",          "z:", "lcs_shading_ori=",
                 "-u", "--input_user_matrix_file",   "u:", "input_user_matrix_file=",
                 "-U", "--user_matrix_print",        "U:", "user_matrix_print=",
                 "-P", "--plot_size",                "P:", "plot_size=",
                 "-A", "--line_width",               "A:", "line_width=",
                 "-B", "--line_col_for",             "B:", "line_col_for=",
                 "-C", "--line_col_rev",             "C:", "line_col_rev=",
                 "-D", "--x_label_pos",              "D:", "x_label_pos=",
                 "-E", "--label_size",               "E:", "label_size=",
                 "-F", "--spacing",                  "F:", "spacing=",
                 "-L", "--length_scaling",           "L:", "length_scaling=",
                 "-M", "--mirror_y_axis",            "M:", "mirror_y_axis=",
                 "-R", "--representation",           "R:", "representation=",
                 "-T", "--title_length",             "T:", "title_length=",
                 "-h", "--help",                     "h",  "help",
                 "-v", "--verbose",                  "v",  "verbose"]

    arguments_sysargv   = tuple(arguments[0::4] + arguments[1::4])
    arguments_opts      = "".join(arguments[2::4])
    arguments_args      = arguments[3::4]


    # setting defaults
    ######################################

    auto_fas                  = False # 0
    input_fasta               = [] 
    output_file_prefix        = None
    collage_output            = True  # 1
    m_col                     = 4
    n_row                     = 5
    filetype                  = 0
    type_nuc                  = True
    input_gff_files           = []
    gff_color_config_file     = ""

    wordsize                  = 10
    plotting_modes            = [0]
    wobble_conversion         = False # 0
    substitution_count        = 0
    rc_option                 = True  # 1
    alphabetic_sorting        = False # 0
    only_vs_first_seq         = False # 0

    lcs_shading               = False # 0
    lcs_shading_num           = 4
    lcs_shading_ref           = 0
    lcs_shading_interval_len  = 50   # interval default changes to "10" for amino acids [type_nuc = n]
    lcs_shading_ori           = 0

    input_user_matrix_file    = ""
    user_matrix_print         = False

    plot_size                 = 10
    line_width                = 1
    line_col_for              = "black"
    line_col_rev              = "#009243"
    x_label_pos               = True # 0
    label_size                = 10
    spacing                   = 0.04
    length_scaling            = False # 0
    title_length              = 20  # float("Inf")
    title_clip_pos            = "B" # B (begin), E (end)
    max_N_percentage          = 49  # fixed value, no user input 
    mirror_y_axis             = False
    representation            = 0
   
    aa_bp_unit                = "bp"

    verbose                   = False # 0

    filetype_dict             = {0: "png",                1: "pdf",                       2: "svg"}
    lcs_shading_ref_dict      = {0: "maximal LCS length", 1: "maximally possible length", 2: "given interval sizes"}
    plotting_mode_dict        = {0: "self",               1: "paired",                    2: "all-against-all"}
    lcs_shading_ori_dict      = {0: "forward",            1: "reverse complement",        2: "both"}
    representation_dict       = {0: "full",               1: "upper",                     2: "lower"}

    # return default parameters for testing purposes
    if trial_mode:
        print "ATTENTION: YOU ARE IN THE TRIAL MODE!!!\n\n"

        commandline = "trial_mode\n"

        parameters = [commandline, auto_fas, input_fasta, output_file_prefix, collage_output, m_col, n_row, filetype_dict[filetype], type_nuc, input_gff_files, gff_color_config_file, wordsize, plotting_modes, wobble_conversion, substitution_count, rc_option, alphabetic_sorting, only_vs_first_seq, lcs_shading, lcs_shading_num, lcs_shading_ref, lcs_shading_interval_len, lcs_shading_ori, input_user_matrix_file, user_matrix_print, plot_size, line_width, line_col_for, line_col_rev, x_label_pos, label_size, spacing, length_scaling, title_length, title_clip_pos, max_N_percentage, mirror_y_axis, representation, verbose]
        return parameters


    # read arguments
    ######################################

    commandline = ""
    for arg in sys.argv:
        commandline += arg + " "

    log_txt = "\n...reading input arguments..."
    print log_txt

    if len(sys.argv) < 2:
        print "\nERROR: More arguments are needed. Exit..."
        log_txt += "\nERROR: More arguments are needed. Exit..."
        usage()
        sys.exit()

    elif sys.argv[1] not in arguments_sysargv:
        print "\nINPUT ERROR: Input argument %s unknown. Please check the help screen." % sys.argv[1] 
        log_txt += "\nINPUT ERROR: Input argument %s unknown. Please check the help screen." % sys.argv[1] 
        # usage()
        sys.exit()

    try:
        opts, args = getopt.getopt(sys.argv[1:], arguments_opts, arguments_args)

    except getopt.GetoptError:
        print "\nINPUT ERROR (getopt): Input argument %s unknown. Please check the help screen." % sys.argv[1:]
        log_txt += "\nINPUT ERROR (getopt): Input argument %s unknown. Please check the help screen." % sys.argv[1:]
        # usage()
        sys.exit()

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print "...fetch help screen"
            log_txt += "\n...fetch help screen"
            usage(), sys.exit()

        if opt in ("-v", "--verbose"):
            print "...verbose output"
            log_txt += "\n...verbose output"
            verbose            = True

        elif opt in ("-i", "--input_fasta"):
            if "," in arg:
                arg_list = arg.split(",")
                for temp_file in arg_list: 
                    if not os.path.exists(str(temp_file)):
                        message = "\nERROR: fasta_file '%s' was not found!" % str(temp_file)
                        sys.exit(message)
                    else:
                        input_fasta.append(str(temp_file))
                        print "fasta file #%i: %s" % (len(input_fasta), str(temp_file))
                        log_txt += "\nfasta file #%i: %s" % (len(input_fasta), str(temp_file))
            else:
                if not os.path.exists(str(arg)):
                    message = "\nERROR: fasta_file '%s' was not found!" % str(arg)
                    log_txt += message
                    sys.exit(message)
                else:
                    input_fasta.append(str(arg))
                    print "fasta file #%i: %s" % (len(input_fasta), str(arg))
                    log_txt += "\nfasta file #%i: %s" % (len(input_fasta), str(arg))


        elif opt in ("-a", "--auto_fas"):
            auto_fas           = True


        # multiple gff files: reads them into a list
        elif opt in ("-g", "--input_gff_files"):

            # append gff file only if existing 
            if "," in arg:
                arg_list = arg.split(",")
                for temp_file in arg_list: 
                    if not os.path.exists(str(temp_file)):
                        message = "\nERROR: gff_file '%s' was not found!" % str(temp_file)
                        print message
                        log_txt += message
                        print "  -->Running FlexiDot without this gff file!"
                        log_txt += "\n  -->Running FlexiDot without this gff file!"
                    else:
                        print "GFF file #%i: %s" %(len(input_gff_files), str(temp_file))
                        log_txt += "\nGFF file #%i: %s" %(len(input_gff_files), str(temp_file))
                        input_gff_files.append(str(temp_file))
            else:
                if not os.path.exists(str(arg)):
                    message = "\nERROR: gff_file '%s' was not found!" % str(arg)
                    print message
                    log_txt += message
                    print "  -->Running FlexiDot without this gff file!"
                    log_txt += "\n  -->Running FlexiDot without this gff file!"
                else:
                    input_gff_files.append(str(arg))
                    print "GFF file #%i: %s" %(len(input_gff_files), str(arg))
                    log_txt += "\nGFF file #%i: %s" %(len(input_gff_files), str(arg))


        elif opt in ("-G", "--gff_color_config_file"):
            if not os.path.exists(str(arg)):
                message = "\nERROR: gff_color_config_file '%s' was not found!" % str(arg)
                print message + "\n  -->Running FlexiDot with default gff coloring specification!"
                log_txt += message + "\n  -->Running FlexiDot with default gff coloring specification!"
            else:
                gff_color_config_file = str(arg)


        elif opt in ("-u", "--input_user_matrix_file"):
            if not os.path.exists(str(arg)):
                message = "\nERROR: input_user_matrix_file '%s' was not found!" % str(arg)
                print message + "\n  -->Running FlexiDot without input_user_matrix_file %s!" % arg
                log_txt += message + "\n  -->Running FlexiDot withdefault matrix shading file!"
            else:
                input_user_matrix_file = str(arg)

        elif opt in ("-U", "--user_matrix_print"):
            user_matrix_print     = check_bools(str(arg), default=user_matrix_print)

        elif opt in ("-o", "--output_file_prefix"):
            output_file_prefix = arg

        elif opt in ("-c", "--collage_output"):
            collage_output     = check_bools(str(arg), default=collage_output)

        elif opt in ("-m", "--m_col"):
            try: m_col              = int(arg)
            except: 
                print "m_col - invalid argument - using default value"
                log_txt += "\nm_col - invalid argument - using default value"

        elif opt in ("-n", "--n_row"):
            try: n_row              = int(arg)
            except: 
                print "n_row - invalid argument - using default value"
                log_txt += "\nn_row - invalid argument - using default value"

        elif opt in ("-f", "--filetype"):
            if 0 <= int(arg) <= 2:
                filetype    = int(arg)
            else: 
                print "\nERROR: Please provide valid filetype argument. %s is out of range. It will be set to -f 0 [default]." %(filetype)
                log_txt += "\nERROR: Please provide valid filetype argument. %s is out of range. It will be set to -f 0 [default]." %(filetype)

        elif opt in ("-t", "--type_nuc"):
            type_nuc = check_bools(str(arg), default=type_nuc)

            if type_nuc == False:
                # interval default changed for amino acids
                lcs_shading_interval_len = 10
                aa_bp_unit = "aa"

        elif opt in ("-k", "--wordsize"):
            try: wordsize           = int(arg)
            except: 
                print "wordsize - invalid argument - using default value"
                log_txt += "\nwordsize - invalid argument - using default value"

        elif opt in ("-p", "--plotting_mode"):
            if "," in arg:
                temp_modes = arg.split(",")
                for item in temp_modes:
                    if item in ["0","1","2"]:
                        plotting_modes.append(int(item))
            elif arg in ["0","1","2"]:
                plotting_modes = [int(arg)]
            else:
                print "Please provide valid plotting_modes argument - e.g. 1 or 0,1,2 - using default [0]"
                log_txt += "\nPlease provide valid plotting_modes argument - e.g. 1 or 0,1,2 - using default [0]"

        elif opt in ("-w", "--wobble_conversion"):
            wobble_conversion  = check_bools(str(arg), default=wobble_conversion)

        elif opt in ("-S", "--substitution_count"):
            try: substitution_count         = int(arg)
            except: 
                print "substitution_count - invalid argument - using default value"
                log_txt += "\nsubstitution_count - invalid argument - using default value"

        elif opt in ("-r", "--rc_option"):
            rc_option          = check_bools(str(arg), default=rc_option)

        elif opt in ("-s", "--alphabetic_sorting"):
            alphabetic_sorting = check_bools(str(arg), default=alphabetic_sorting)

        elif opt in ("-O", "--only_vs_first_seq"):
            only_vs_first_seq = check_bools(str(arg), default=only_vs_first_seq)

        elif opt in ("-x", "--lcs_shading"):
            lcs_shading        = check_bools(str(arg), default=lcs_shading)

        elif opt in ("-X", "--lcs_shading_num"):
            try: lcs_shading_num    = int(arg) - 1
            except: 
                print "lcs_shading_num - invalid argument - using default value"
                log_txt += "\nlcs_shading_num - invalid argument - using default value"

        elif opt in ("-y", "--lcs_shading_ref"):
            try:
                if 0 <= int(arg) <= 2: 
                    lcs_shading_ref    = int(arg)
                else: 
                    print "\nERROR: lcs_shading_ref %s out of range. It will be set to -y 0 [default]." %(lcs_shading_ref)
                    log_txt += "\nERROR: lcs_shading_ref %s out of range. It will be set to -y 0 [default]." %(lcs_shading_ref)
            except: 
                print "lcs_shading_ref - invalid argument - using default value"
                log_txt += "\nlcs_shading_ref - invalid argument - using default value"

        elif opt in ("-Y", "--lcs_shading_interval_len"):
            try: lcs_shading_interval_len = int(arg)
            except: 
                print "lcs_shading_interval_len - invalid argument - using default value"
                log_txt += "\nlcs_shading_interval_len - invalid argument - using default value"

        elif opt in ("-z", "--lcs_shading_ori"):
            if 0 <= int(arg) <= 2:
                lcs_shading_ori    = int(arg)
            else: 
                print "\nERROR: Please provide valid lcs_shading_ori argument. %s is out of range. It will be set to -z 0 [default]." %(lcs_shading_ori)
                log_txt += "\nERROR: Please provide valid lcs_shading_ori argument. %s is out of range. It will be set to -z 0 [default]." %(lcs_shading_ori)

        elif opt in ("-P", "--plot_size"):
            try: plot_size         = float(arg)
            except: 
                print "plot_size - invalid argument - using default value"
                log_txt += "\nplot_size - invalid argument - using default value"


        elif opt in ("-A", "--line_width"):
            try: line_width         = float(arg)
            except: 
                print "line_width - invalid argument - using default value"
                log_txt += "\nline_width - invalid argument - using default value"

        elif opt in ("-B", "--line_col_for"):
            if mcolors.is_color_like(arg):
                line_col_for       = arg
            else: 
                print "line_col_for - invalid argument - using default value"
                log_txt += "\nline_col_for - invalid argument - using default value"

        elif opt in ("-C", "--line_col_rev"):
            if mcolors.is_color_like(arg):
                line_col_rev       = arg
            else: 
                print "line_col_rev - invalid argument - using default value"
                log_txt += "\nline_col_rev - invalid argument - using default value"

        elif opt in ("-D", "--x_label_pos"):
            x_label_pos        = check_bools(str(arg), default=x_label_pos)

        elif opt in ("-E", "--label_size"):
            try: label_size         = float(arg)
            except: 
                print "label_size - invalid argument - using default value"
                log_txt += "\nlabel_size - invalid argument - using default value"

        elif opt in ("-F", "--spacing"):
            try: spacing            = float(arg)
            except: 
                print "spacing - invalid argument - using default value"
                log_txt += "\nspacing - invalid argument - using default value"

        elif opt in ("-L", "--length_scaling"):
            length_scaling     = check_bools(str(arg), default=length_scaling)

        elif opt in ("-M", "--mirror_y_axis"):
            mirror_y_axis     = check_bools(str(arg), default=mirror_y_axis)

        elif opt in ("-R", "--representation"):
            if 0 <= int(arg) <= 2:
                representation    = int(arg)
            else: 
                print "\nERROR: Please provide valid representation argument. %s is out of range. It will be set to -R 0 [default]." %(representation)
                log_txt += "\nERROR: Please provide valid representation argument. %s is out of range. It will be set to -R 0 [default]." %(representation)

        elif opt in ("-T", "--title_length"):
            try: title_length     = int(arg)
            except:
                try: 
                    title_length = int(str(arg)[:-1])
                    if arg[-1].upper() in ["B", "E"]: # B (beginning), E (end)
                        title_clip_pos = arg[-1].upper() 
                    else:
                        print "title_length position information invalid - using default value"
                        log_txt += "\ntitle_length position information invalid - using default value"
                except:
                    print "title_length - invalid argument - using default value"
                    log_txt += "\ntitle_length - invalid argument - using default value"

    # start logging file
    logprint(commandline, start=True, printing=False, prefix=output_file_prefix)
    logprint(log_txt, start=False, printing=False)


    # print chosen arguments
    ######################################

    text = "\n%s\n" % (70 * "-")
    text += "\n" + "INPUT/OUTPUT OPTIONS...\n"
    text += "\n" + "Input fasta file:                                  " + ", ".join(input_fasta)
    text += "\n" + "Automatic fasta collection from current directory: " + str(auto_fas)
    text += "\n" + "Collage output:                                    " + str(collage_output)
    text += "\n" + "Number of columns per page:                        " + str(m_col)
    text += "\n" + "Number of rows per page:                           " + str(n_row)
    text += "\n" + "File format:                                       " + filetype_dict[filetype]
    text += "\n" + "Residue type is nucleotide:                        " + str(type_nuc)

    text += "\n" + "\n\nCALCULATION PARAMETERS...\n"
    text += "\n" + "Wordsize:                                          " + str(wordsize)
    text += "\n" + "Sustitution count:                                 " + str(substitution_count)
    text += "\n" + "Plotting mode:                                     " + str(plotting_modes).replace("[", "").replace("]", "") + "\n" + 51 * " "
    for item in plotting_modes:
        text += plotting_mode_dict[item] + " "
    text += "\n" + "Ambiguity handling:                                " + str(wobble_conversion)
    text += "\n" + "Reverse complement scanning:                       " + str(rc_option)
    text += "\n" + "Alphabetic sorting:                                " + str(alphabetic_sorting)
    
    if 1 in plotting_modes:
        text += "\n" + "Only matching sequences to first entry:            " + str(only_vs_first_seq)

    if 0 in plotting_modes and input_gff_files != []:
        text += "\n" + "Input gff files:                                   " + ", ".join(input_gff_files)
        if gff_color_config_file != "": 
            text += "\n" + "GFF color config file:                             " + gff_color_config_file
    text += "\n" + "Prefix for output files:                           " + str(output_file_prefix)

    if 2 in plotting_modes:
        text += "\n" + "\n\nLCS SHADING OPTIONS (plotting_mode 'all-against-all' only)...\n"
        text += "\n" + "LCS shading:                                       " + str(lcs_shading)
        text += "\n" + "LCS shading interval number:                       " + str(lcs_shading_num + 1)
        text += "\n" + "LCS shading reference:                             " + lcs_shading_ref_dict[lcs_shading_ref]
        if lcs_shading_ref == 2:
            text += "\n" + "LCS shading interval size [%s]:                    " % (aa_bp_unit) + str(lcs_shading_interval_len)
        text += "\n" + "LCS shading orientation:                           " + lcs_shading_ori_dict[lcs_shading_ori]
        if input_user_matrix_file != "":
            text += "\n" + "Custom user shading matrix file:                   " + input_user_matrix_file
            text += "\n" + "Print user matrix values (instead of dotplot):     " + str(user_matrix_print)
        text += "\n" + "Displayed plot region:                             " + representation_dict[representation]

    text += "\n" + "\n\nGRAPHIC FORMATTING...\n"
    text += "\n" + "Plot size:                                         " + str(plot_size)
    text += "\n" + "Line width:                                        " + str(line_width)
    text += "\n" + "Line color:                                        " + line_col_for
    text += "\n" + "Reverse line color:                                " + line_col_rev
    text += "\n" + "X label position:                                  " + str(x_label_pos)
    text += "\n" + "Label size:                                        " + str(label_size)
    text += "\n" + "Spacing:                                           " + str(spacing)
    if mirror_y_axis:
        text += "\n" + "Y-axis mirrored (bottom to top)                   " + str(mirror_y_axis)
    if title_clip_pos == "E": 
        text += "\n" + "Title length (limit number of characters):         " + "last" + str(title_length) + "characters"
    else:
        text += "\n" + "Title length (limit number of characters):         " + "first" + str(title_length) + "characters"
    text += "\n" + "Length scaling:                                    " + str(length_scaling)
    text += "\n%s\n" % (70 * "-")
    logprint(text)


    # collect settings
    parameters = [commandline, auto_fas, input_fasta, output_file_prefix, collage_output, m_col, n_row, filetype_dict[filetype], type_nuc, input_gff_files, gff_color_config_file, wordsize, plotting_modes, wobble_conversion, substitution_count, rc_option, alphabetic_sorting, only_vs_first_seq, lcs_shading, lcs_shading_num, lcs_shading_ref, lcs_shading_interval_len, lcs_shading_ori, input_user_matrix_file, user_matrix_print, plot_size, line_width, line_col_for, line_col_rev, x_label_pos, label_size, spacing, length_scaling, title_length, title_clip_pos, max_N_percentage, mirror_y_axis, representation, verbose]

    return parameters


###############################
#      Helper Functions       #
###############################

def alphabets(type_nuc=True):
    """
    provide ambiguity code for sequences
    """

    nucleotide_alphabet = ["A", "C", "G", "T"]

    nucleotide_alphabet_full = ["A", "C", "G", "T", "N", "B", "D", "H", 
                                "V", "Y", "R", "W", "S", "K", "M"]

    nucleotide_ambiguity_code = {"N": ["A", "C", "G", "T"], # any
                                 "B": ["C", "G", "T"], # not A
                                 "D": ["A", "G", "T"], # not C
                                 "H": ["A", "C", "T"], # not G
                                 "V": ["A", "C", "G"], # not T
                                 "Y": ["C", "T"], # pyrimidine
                                 "R": ["A", "G"], # purine
                                 "W": ["A", "T"], # weak
                                 "S": ["C", "G"], # strong
                                 "K": ["G", "T"], # keto
                                 "M": ["A", "C"]} # amino

    nucleotide_match_dict = {"N": "[ACGTNBDHVYRWSKM]", # any
                             "B": "[CGTNBDHVYRWSKM]", # not A
                             "D": "[AGTNBDHVYRWSKM]", # not C
                             "H": "[ACTNBDHVYRWSKM]", # not G
                             "V": "[ACGNBDHVYRWSKM]", # not T
                             "K": "[GTNBDHVYRWSK]", # keto - not A,C,M
                             "M": "[ACNBDHVYRWSM]", # amino - not G,T,K
                             "W": "[ATNBDHVYRWKM]", # weak - not C,G,S
                             "S": "[CGNBDHVYRSKM]", # strong - not A,G,W
                             "Y": "[CTNBDHVYWSKM]", # pyrimidine - not A,G,R
                             "R": "[AGNBDHVRWSKM]", # purine - not C,T,Y
                             "A": "[ANDHVRWM]",
                             "C": "[CNBHVYSM]",
                             "G": "[GNBDVRSK]",
                             "T": "[TNBDHYWK]"} 

    # nucleotide_match_dict = {"N": ".", # any
    #                          "B": "[^A]", # not A
    #                          "D": "[^C]", # not C
    #                          "H": "[^G]", # not G
    #                          "V": "[^T]", # not T
    #                          "K": "[^ACM]", # keto - not A,C,M
    #                          "M": "[^GTK]", # amino - not G,T,K
    #                          "W": "[^CGS]", # weak - not C,G,S
    #                          "S": "[^AGW]", # strong - not A,G,W
    #                          "Y": "[^AGR]", # pyrimidine - not A,G,R
    #                          "R": "[^CTY]", # purine - not C,T,Y
    #                          "A": "[ANDHVRWM]",
    #                          "C": "[CNBHVYSM]",
    #                          "G": "[GNBDVRSK]",
    #                          "T": "[TNBDHYWK]"} 

    aminoacid_alphabet = ["A", "R", "N", "D", "C", "E", "Q", "G", 
                          "H", "I", "L", "K", "M", "F", "P", "S", 
                          "T", "W", "Y", "V", "U", "O", "*"]

    aminoacid_alphabet_full = ["A", "R", "N", "D", "C", "E", "Q", "G", 
                               "H", "I", "L", "K", "M", "F", "P", "S", 
                               "T", "W", "Y", "V", "U", "O", "*", "J", 
                               "Z", "B", "X"]

    aminoacid_ambiguity_code  = {"J": ["I", "L"],
                                 "Z": ["Q", "E"], 
                                 "B": ["N", "D"], 
                                 "X": ["A", "R", "N", "D", "C", "E", "Q", "G", 
                                       "H", "I", "L", "K", "M", "F", "P", "S", 
                                       "T", "W", "Y", "V", "U", "O", "*"]} # any

    aminoacid_match_dict = {"J": "[ILJ]",
                            "Z": "[QEZ]", 
                            "B": "[NDB]", 
                            # "X": ".",
                            "X": "[ARNDCEQGHILKMFPSTWYVUO*XBZJ]",
                            "A": "[AX]",
                            "R": "[RX]",
                            "N": "[NXB]",
                            "D": "[DXB]",
                            "C": "[CX]",
                            "E": "[EXZ]",
                            "Q": "[QXZ]",
                            "G": "[GX]",
                            "H": "[HX]",
                            "I": "[IXJ]",
                            "L": "[LXJ]",
                            "K": "[KX]",
                            "M": "[MX]",
                            "F": "[FX]",
                            "P": "[PX]",
                            "S": "[SX]",
                            "T": "[TX]",
                            "W": "[WX]",
                            "Y": "[YX]",
                            "V": "[VX]",
                            "U": "[UX]",
                            "O": "[OX]",
                            "*": "[*X]"}

    aa_only = set(['E', 'F', 'I', 'J', 'L', 'O', 'Q', 'P', 'U', 'X', 'Z', '*'])
    # return nucleotide_alphabet, nucleotide_alphabet_full, nucleotide_ambiguity_code, aminoacid_alphabet, aminoacid_alphabet_full, aminoacid_ambiguity_code, aa_only

    if type_nuc:
        return nucleotide_alphabet, nucleotide_alphabet_full, nucleotide_ambiguity_code, nucleotide_match_dict
    else:
        return aminoacid_alphabet, aminoacid_alphabet_full, aminoacid_ambiguity_code, aminoacid_match_dict

def logprint(text, start=False, printing=True, prefix=""):
    """
    log output to log_file and optionally print
    """

    # define log file name and open file
    global log_file_name
    if start and trial_mode:
        log_file_name = "log_file.txt"
        if prefix != "" and prefix != None:
            if not prefix.endswith("-"):
                prefix = prefix + "-"
            log_file_name = prefix + log_file_name  
        log_file = open(log_file_name, 'w')
        log_file.write("Date: %s\n\n" % str(datetime.datetime.now()))
    elif start:
        date = datetime.date.today()
        time = str(datetime.datetime.now()).split(" ")[1].split(".")[0].replace(":", "-")
        log_file_name = "%s_%s_log_file.txt" % (date, time)
        if prefix != "" and prefix != None:
            if not prefix.endswith("-"):
                prefix = prefix + "-"
            log_file_name = prefix + log_file_name  
        log_file = open(log_file_name, 'w')
        log_file.write("Date: %s\n\n" % str(datetime.datetime.now()))
    else:
        log_file = open(log_file_name, 'a')

    # write log (and print)
    log_file.write(text + "\n")
    if printing:
        print text
    log_file.close()

def time_track(starting_time, show=True):
    """
    calculate time passed since last time measurement
    """
    now = time.time()
    delta = now - starting_time
    if show:
        text = "\n\t %s seconds\n" % str(delta)
        logprint(text, start=False, printing=True)
    return now

def calc_fig_ratio(ncols, nrows, plot_size, verbose=False):
    """
    calculate size ratio for given number of columns (ncols) and rows (nrows) 
    with plot_size as maximum width and length
    """
    ratio = ncols*1./nrows
    if verbose:
        text = " ".join([ncols, nrows, ratio])
        logprint(text, start=False, printing=True)
    if ncols >= nrows:
        figsize_x = plot_size
        figsize_y = plot_size / ratio
    else:
        figsize_x = plot_size * ratio
        figsize_y = plot_size
    return figsize_x, figsize_y

def shorten_name(seq_name, max_len=20, title_clip_pos="B"): #, delim="_"):
    """
    shorten sequence names (for diagram titles)
    """

    if len(seq_name) <= max_len:
        return seq_name

    # take last characters
    if title_clip_pos == "E":
        name = seq_name[len(seq_name)-max_len:]

    # take first characters
    else:
        name = seq_name[:max_len]

    """# keep first and last part if multiple parts separated by delimiter (e.g. species_prefix + sequence_id)
    if delim in seq_name:
        if seq_name.count(delim) >= 2:
            name = "%s..." % delim.join(seq_name.split(delim)[:1]) + seq_name.split(delim)[-1] # .replace("_000", "-")
        else:
            name = seq_name[:((max_len-2)//2)] + "..." + seq_name[((max_len-2)//2):]

        if len(name) > max_len:
            name = name[:((max_len-2)//2)] + "..." + name[((max_len-2)//2):]
    else:
        name = seq_name[:((max_len-2)//2)] + "..." + seq_name[((max_len-2)//2):]
    """

    return name

def unicode_name(name):
    """
    replace non-ascii characters in string (e.g. for use in matplotlib)
    """
    unicode_string = eval('u"%s"' % name)
    return unicodedata.normalize('NFKD', unicode_string).encode('ascii','ignore')

def check_bools(arg, update_log_txt = True, default=None):
    """
    converts commandline arguments into boolean
    """


    # convert valid arguments
    if str(arg).lower() == "y" or str(arg) == "1":  
        return True
    elif str(arg).lower() == "n" or str(arg) == "0":  
        return False

    # use default in case of invalid argument
    else:
        if update_log_txt:
            global log_txt
            log_txt += "using default for " + str(arg)
        else:
            try:
                logprint("using default for " + str(arg))
            except:
                print "using default for " + str(arg)
        return default

def create_color_list(number, color_map=None, logging=False, max_grey="#595959"):
    """
    create color list with given number of entries 
    grey by default, matplotlib color_map can be provided 
    """

    try:
        # create pylab colormap
        cmap = eval("P.cm." + color_map)
        # get descrete color list from pylab
        cmaplist = [cmap(i) for i in range(cmap.N)] # extract colors from map
        # determine positions for number of colors required
        steps = (len(cmaplist)-1)/(number)
        numbers = range(0, len(cmaplist), steps)

        # extract color and convert to hex code
        colors = []
        for idx in numbers[:-1]:
            rgb_color = cmaplist[idx]
            col = rgb2hex(rgb_color[0]*255, rgb_color[1]*255, rgb_color[2]*255)
            colors.append(col)

    # grey
    except:
        if not color_map == None:
            logprint("Invalid color_map (%s) provided! - Examples: jet, Blues, OrRd, bwr,..." % color_map)
            logprint("See https://matplotlib.org/users/colormaps.html\n")
        old_max_grey = "#373737"
        old_max_grey = "#444444"
        colors = list(Color("#FFFFFF").range_to(Color(max_grey), number)) # grey
        for idx in range(len(colors)): 
            colors[idx] = str(colors[idx]).replace("Color ", "") 
            if "#" in colors[idx] and len(colors[idx]) != 7:
                # print colors[idx]
                colors[idx] = colors[idx] + colors[idx][-(7-len(colors[idx])):]

    text = "%d Colors: %s" % (len(colors), ", ".join(colors))
    if logging: logprint(text, start=False, printing=True)

    if len(colors) < number:
        logprint("\nError in color range definition! %d colors missing\n" % (number - len(colors)))

    return colors


###############################
#        File Handling        #
###############################

def read_seq(input_fasta, verbose=False):
    """
    read fasta sequences from (all) file(s)
    """

    # check if file provided
    if input_fasta == [] or input_fasta == "":
        text = "Attention: No valid file names provided: >%s<" % input_fasta
        logprint(text, start=False, printing=True)
        return {}, []

    # combine sequence files, if required
    if type(input_fasta) == list:
        # concatenate fasta files
        if len(input_fasta) > 1:
            if verbose:
                print "concatenating fastas...",
                text = "concatenating fastas..."
            input_fasta_combi = concatenate_files(input_fasta)
            if verbose:
                print "done"
                text += "done"
                logprint(text, start=False, printing=False)
        else:
            input_fasta_combi = input_fasta[0]
    else:
        input_fasta_combi = input_fasta

    # read sequences
    if verbose:
        print "reading fasta...",
        text = "reading fasta...",
    try:
        seq_dict = SeqIO.index(input_fasta_combi, "fasta") 
    except ValueError:
        logprint("Error reading fasta sequences - please check input files, e.g. for duplicate names!")
        return {}, []
    except:
        logprint("Error reading fasta sequences - please check input files!")
        return {}, []

    if verbose:
        print "done"
        text += "done"
        logprint(text, start=False, printing=False)

    for seq in seq_dict:
        if "-" in seq_dict[seq].seq:
            # ungapped = seq_dict[seq].seq.ungap("-") # cannot be assigned back to sequence record
            text = "\nSequences degapped prior Analysis!!!"
            logprint(text, start=False, printing=True)
            return read_seq(degap_fasta(input_fasta), verbose=verbose)

    # get ordered sequence names
    sequences = []
    for item in SeqIO.parse(input_fasta_combi, "fasta"):
        sequences.append(item.id)
    return seq_dict, sequences

def read_gff_color_config(gff_color_config_file=""):
    """
    define coloring options for gff-based color shading of self-dotplots
    """

    # default aestetics for annotation shading (e.g. if no user config file is provided)
    # dictionary with feature_type as key and tuple(color, transparency, zoom) as value
    gff_feat_colors = {"orf":                     ("#b41a31",    0.2,  0),
                       "orf_rev":                 ("#ff773b",    0.3,  0),
                       "gene":                    ("#b41a31",    0.2,  0),
                       "cds":                     ("darkorange", 0.2,  0),
                       "exon":                    ("orange",     0.2,  0),
                       "intron":                  ("lightgrey",  0.2,  0),
                       "utr":                     ("lightblue",  0.2,  0),
                       "repeat_region":           ("green",      0.3,  0),
                       "repeat":                  ("green",      0.3,  0),
                       "tandem_repeat":           ("red",        0.3,  0),
                       "transposable_element":    ("blue",       0.3,  0),
                       "ltr_retrotransposon":     ("#cccccc",    0.5,  0),
                       "ltr-retro":               ("#cccccc",    0.5,  0),
                       "long_terminal_repeat":    ("#2dd0f0",    0.75, 2),
                       "ltr":                     ("#2dd0f0",    0.75, 2),
                       "pbs":                     ("purple",     0.75, 2),
                       "ppt":                     ("#17805a",    0.5,  2),
                       "target_site_duplication": ("red",        0.75, 2),
                       "misc_feature":            ("grey",       0.3,  0),
                       "misc_feat":               ("grey",       0.3,  0),
                       "misc":                    ("grey",       0.3,  0),
                       "others":                  ("grey",       0.5,  0)}
    if gff_color_config_file in ["", None] or not os.path.exists(str(gff_color_config_file)):
        return gff_feat_colors

    text = "Updating GFF color configuration with custom specifications\n"
    logprint(text, start=False, printing=True)

    # read custom gff_color_config_file 
    in_file = open(gff_color_config_file, 'rb')
    overwritten = set([])
    for line in in_file:
        if not line.startswith("#") and len(line.strip().split("\t")) >= 4:
            data = line.strip().split("\t")
            feat  = data[0].lower()
            color = data[1].lower()

            # check, if settings are valid
            if not mcolors.is_color_like(color):
                color = "grey"
                text = "Invalid color specified for %s: %s - default grey" % (data[0], data[1])
                logprint(text)
            try:
                alpha = float(data[2])
            except:
                alpha = 0.75
                text = "Invalid alpha specified for %s: %s - default 0.75" % (data[0], data[2])
                logprint(text)
            try:
                zoom  = float(data[3])
            except:
                zoom = 0
                text = "Invalid zoom specified for %s: %s - default 0" % (data[0], data[3])
                logprint(text)

            # track changes of predefined settings
            if feat in gff_feat_colors.keys():
                overwritten.add(data[0].lower()) 

            gff_feat_colors[feat] = (color, alpha, zoom)
    in_file.close()

    # default coloring for unknown annotations
    if not "others" in gff_feat_colors.keys():
        gff_feat_colors["others"] = ("grey", 0.5, 0)

    if verbose:
        # print configuration
        text = "\n\nGFF color specification:\n%s\n" % (60 * ".")
        for item in sorted(gff_feat_colors.keys()):
            text += "%-30s\t%-10s\t%-5s\t%s\n" % (item, str(gff_feat_colors[item][0]), str(gff_feat_colors[item][1]), str(gff_feat_colors[item][2]))
        logprint (text, printing=True) 

    # print overwritting feature type specifications
    if len(overwritten) != 0:
        text = "%d feature type specifications overwritten:" % len(overwritten)
        text += "\n\t"+ ", ".join(overwritten) + "\n"
        logprint(text, start=False, printing=True)

    text = "GFF color specification updated acc. to %s\n\t%s\n\n" % (gff_color_config_file, ", ".join(gff_feat_colors))
    logprint(text, start=False, printing=True)

    return gff_feat_colors

def read_gffs(input_gff_files, color_dict={"others": ("grey", 1, 0)}, type_nuc=True, prefix="", filetype='png', verbose=False):
    """
    create feature dictionary from input_gff
    sequence name as key and (feature type, start, stop) as value
    """
    if type(input_gff_files) != list:
        input_gff_files = [input_gff_files]

    # create dictionary with seq_name as key and (type, start and stop) as value 
    unknown_feats = set([])
    used_feats = set([])
    feat_dict = {}
    for input_gff in input_gff_files:
        text = "...reading " + input_gff
        logprint(text, start=False, printing=True)

        in_file = open(input_gff, 'rb')
        for line in in_file:
            if not line.startswith("#") and line.strip() != "":
                data = line.strip().split("\t")
                feat_type = data[2].lower()
                if data[6] == "-":
                    feat_type += "_rev"
                if not feat_type.lower() in color_dict.keys():
                    if feat_type.lower().replace("_rev", "") in color_dict.keys():
                        feat_type = feat_type.replace("_rev", "")
                    else:
                        unknown_feats.add(feat_type)
                        feat_type = "others"
                used_feats.add(feat_type)
                if not data[0] in feat_dict.keys():
                    feat_dict[data[0]] = [(feat_type, int(data[3]), int(data[4]))] # feature type, start, stop
                else:
                    feat_dict[data[0]].append((feat_type, int(data[3]), int(data[4]))) # feature type, start, stop
        if verbose:
            text = "\nAnnotations for: %s\n" % ", ".join(feat_dict.keys()[:10])
            if len(feat_dict.keys()) > 10:
                text = text[:-1] + ", ...\n"
            logprint(text, start=False, printing=True)
        in_file.close()

    # print feature types without specific shading settings 
    if len(unknown_feats) != 0:
        text = "Missing shading specification for %d feature type(s):\n\t%s\n" % (len(unknown_feats), ", ".join(sorted(unknown_feats)))
        logprint(text, start=False, printing=True)

    # create color legend
    colors, alphas = [], []
    for item in sorted(used_feats):
        colors.append(color_dict[item][0])
        alphas.append(color_dict[item][1])
    legend_figure(colors=colors, lcs_shading_num=len(used_feats), type_nuc=type_nuc, bins=sorted(used_feats), alphas=alphas, gff_legend=True, prefix=prefix, filetype=filetype)

    # print settings
    text = "GFF Feature Types: %s\nGFF Colors:        %s" % (", ".join(sorted(used_feats)), ", ".join(sorted(colors)))
    logprint(text, start=False, printing=True)

    return feat_dict

def read_matrix(matrix_file_name, delim="\t", symmetric=True, recursion=False, verbose=False):
    input_file = open(matrix_file_name, 'rb')

    # read sequence names from first column
    names = []
    for line in input_file:
        if not line.startswith("#") and not line.startswith(delim) and delim in line:
            names.append(line.strip().split(delim)[0])
    logprint("Delimiter '%s': %d names - %s\n" % (delim, len(names), ", ".join(names)))

    # check if names were found - otherwise try another delimiter
    if names == [] and not recursion:
        if delim == "\t":
            new_delim = ","
        else:
            new_delim = "\t"
        logprint("\nMatrix file not containing data delimited by '%s' - trying to read matrix with delimiter '%s'" % (delim.replace("\t", "\\t"), new_delim)) 
        info_dict = read_matrix(matrix_file_name, delim=new_delim, symmetric=symmetric, recursion=True, verbose=verbose)
        return info_dict
    elif names == []:
        logprint("Empty matrix file with alternative delimiter!")
        return info_dict
    input_file.close()

    input_file = open(matrix_file_name, 'rb')
    # read matrix entries as values in dictionary with tuple(names) as key
    info_dict = {}
    contradictory_entries = []
    for line in input_file:
        if not line.startswith("#") and not line.startswith(delim) and delim in line:
            data = line.strip().split(delim)
            for idx in range(len(data[1:])):
                # print tuple(sorted([data[0], names[idx]])), data[idx+1]
                if symmetric:
                    key = tuple(sorted([names[idx], data[0]]))
                else:
                    key = tuple(names[idx], data[0])
                if key in info_dict.keys():
                    if symmetric and info_dict[key] != data[idx+1] and data[idx+1] not in ["", "-"] and info_dict[key] not in ["", "-"]:
                        contradictory_entries.append(key)
                info_dict[key] = data[idx+1]
    input_file.close()

    if len(contradictory_entries) != 0:
        try:
            logprint("\nContradictory entries in matrix file %s:\n\t%s" % (matrix_file_name, ", ".join(contradictory_entries)))
        except:
            log_txt = "\nContradictory entries in matrix file %s:\n\t" % (matrix_file_name)
            for item in contradictory_entries:
                log_txt += str(item).replace("'", "") + ", "
            log_txt = log_txt[:-2]
            logprint(log_txt)
        logprint("Using value from bottom left triangle!")
    if verbose:
        logprint("\nMatrix information for Sequences named: " % ", ".join(names))

    return info_dict

def concatenate_files(file_list, combi_filename="temp_combined.fasta", verbose=False):
    """
    concatenate content of all files in file_list into a combined file named combi_filename
    """
    out_file = open(combi_filename, 'w')
    text = ""
    for item in file_list:
        if verbose:
            text += item + " "
            print item, 
        # read in_file linewise and write to out_file
        in_file = open(item, 'rb')
        for line in in_file:
            out_file.write(line.strip()+"\n")
        in_file.close()
    out_file.close()
    if verbose:
        logprint(text, start=False, printing=False)
    return combi_filename

def degap_fasta(input_fasta):
    """
    remove gaps from fasta - new degapped sequence file created
    """

    # degap all sequence files
    output_fastas = []
    if type(input_fasta) != list:
        input_fasta = list(input_fasta)
    for input_fas in input_fasta:
        output_fas = input_fas[:input_fas.rfind(".")] + "_degapped.fas"
        in_file = open(input_fas, 'rb')
        out_file = open(output_fas, 'w')
        for line in in_file:
            if line.startswith(">"):
                out_file.write(line.strip()+"\n")
            else:
                out_file.write(line.strip().replace("-", "")+"\n")
        out_file.close()
        in_file.close()
        output_fastas.append(output_fas)
    return output_fastas

def legend_figure(colors, lcs_shading_num, type_nuc=True, unit="%", filetype="png", max_len=None, min_len=0, bins=[], alphas=[], gff_legend=False, prefix="", verbose=False):
    """
    create figure color legend
    """
    max_legend_length_row = 8
    max_legend_length_col = 4

    # define output file
    if filetype not in ["png", "pdf", "svg"]:
        text = "Provide valid file type - png, pdf, or svg"
        logprint(text, start=False, printing=True)
        filetype="png"                   

    # check if length of information fit
    if not gff_legend and ((bins != [] and len(colors) != lcs_shading_num+1) or (bins != [] and len(colors) != len(bins)+1)):
        if bins != [] and len(colors) != lcs_shading_num+1:
            text = "**Attention**\nlcs_shading_num (%d) does not match number of colors (%d)!\n"% (lcs_shading_num, len(bins))
        elif bins != [] and len(colors) != len(bins)+1:
            text = "**Attention**\nnumber of LCS length bins (%d) does not match number of colors (%d)!\n" % (len(colors), len(bins))
        logprint(text, start=False, printing=True)
    elif gff_legend and len(bins) != len(colors):
        text = "**Attention**\nnumber of GFF Feature Types (%d) does not match number of colors (%d)!\n" % (len(colors), len(bins))
        logprint(text, start=False, printing=True)

    # set alpha values to opaque if none are provided
    if alphas == []:
        for item in colors:
            alphas.append(1)

    # legend data points
    data_points = range(len(colors))
    if not gff_legend:

        # specify intervals, if max_len provided
        if max_len != None:
            multi_factor = 100 # one digit
            if max_len <= 1:
                multi_factor = 1000 # two digits
            # len_interval_size = (max_len-min_len) * multi_factor *1. // lcs_shading_num * (1./ multi_factor)
            len_interval_size = (max_len-min_len) * 1. / lcs_shading_num
            len_pos = [float("%.2f" % (min_len))]
            # calculate interval positions
            for idx in range(lcs_shading_num):
                len_pos.append(float("%.2f" % (len_pos[-1] + len_interval_size)))

            if prefix.startswith("custom-matrix") and (0 <= max_len <= 100 and 0 <= min_len <= 100):
                unit = "%"
            elif prefix.startswith("custom-matrix"):
                unit = ""

            text = "\n%d Legend intervals from %.2f to %.2f: \n\t%s - number: %d, step: %.2f, unit: %s\n" % (lcs_shading_num+1, min_len, max_len, str(len_pos), len(len_pos), len_interval_size, unit)
            logprint(text, start=False, printing=True)
            pos = len_pos
            interval_size = len_interval_size
        # generate legend labels acc. to standard interval notation
        else:
            # use default max_len = 100 and min_len = 0
            len_interval_size = 100. / lcs_shading_num
            pos = [float("%.2f" % (0))]
            # calculate interval positions
            for idx in range(lcs_shading_num):
                pos.append(float("%.2f" % (pos[-1] + len_interval_size)))

            # interval_size = 100 // lcs_shading_num
            # pos = range(interval_size, 101+interval_size, interval_size)

        # remove unneccessary zeros in decimal places (i.e. if x.x00 in all entries)
        while True:
            last_digit_all_zero = True
            no_delim = False
            for idx in range(len(pos)):
                # only process if fraction with decimal places
                if not "." in str(pos[idx]):
                    no_delim = True
                    break
                # only process when all entries end in zero
                elif str(pos[idx])[-1] != "0":
                    last_digit_all_zero = False
                    break
            if not last_digit_all_zero or no_delim:
                break
            # remove last decimal place (== 0) from all entries
            else:
                temp_pos = pos[:]
                for idx in range(len(pos)):
                    if not str(pos[idx])[-2] == ".":
                        pos[idx] = float(str(pos[idx])[:-1])
                    else:
                        pos[idx] = int(str(pos[idx])[:-2])
                logprint("Shortening legend entries: %s - %s" % (temp_pos, pos))

        # eliminate fractions if unit == bp/aa
        if unit in ["aa", "bp"]:
            for idx in range(len(pos)):
                temp_pos = pos[:]
                rounded_unit = False
                if "." in str(pos[idx]):
                    rounded_unit = True
                    # round values up to next integer (keep integer, if not a fraction)
                    pos[idx] = int(pos[idx] / 1) + int(pos[idx] % 1 > 0)
                    if idx == len(pos) - 1 and pos[idx] == 101:
                        pos[idx] = 100
            if rounded_unit:
                logprint("Fractions not permitted for unit '%s': %s -> %s" % (unit, temp_pos, pos))

        if bins != []: # labels provided
            legend_labels = bins[:]
            legend_labels.append("max")
            legend_labels_lengths = []
            for item in bins:
                    legend_labels_lengths.append("[%d %s, %d %s)" % (item - min(bins), unit, item, unit))
            if len(bins) == len(colors) - 1:
                legend_labels_lengths.append("[%d %s, %s]" % (max(bins), unit, u"\u221E")) # infinite

        else:
            legend_labels = []
            legend_labels_lengths = []
            for idx in range(len(pos)):
                num = pos[idx]
                try:
                    legend_labels.append("[%d%%, %d%%)" % (num - interval_size, num))
                except:
                    legend_labels.append("[%d%%, %d%%)" % (num - len_interval_size, num))
                if max_len != None:
                    num = len_pos[idx]
                    # as int or float
                    if num == int(num) and int(len_interval_size) == len_interval_size:
                        legend_labels_lengths.append("[%d %s, %d %s)" % (num, unit, num + len_interval_size, unit))
                    else:
                        legend_labels_lengths.append("[%.2f %s, %.2f %s)" % (num, unit, num + len_interval_size, unit))
            legend_labels[-1] = "100" + unit
            if max_len != None:
                if num == int(num) and int(len_interval_size) == len_interval_size:
                    legend_labels_lengths[-1] = u"[%d %s, \u221E]" % (max_len, unit)
                else:
                    legend_labels_lengths[-1] = u"[%.2f %s, \u221E]" % (max_len, unit)

    # set labels and choose file name
    if gff_legend:
        label_text = bins[:]
        edge_col = None
        legend_file_name = "GFF_Shading_Legend_n%d." % lcs_shading_num + filetype
    elif max_len != None: 
        label_text = legend_labels_lengths[:]
        edge_col = "black"
        legend_file_name = "Polydotplot_LCS_Shading_Legend_max%d%s_n%d." % (max_len, unit, lcs_shading_num+1) + filetype
    elif bins != []:
        label_text = legend_labels_lengths[:]
        edge_col = "black"
        legend_file_name = "Polydotplot_LCS_Shading_Legend_%d%s_n%d." % (bins[0], unit, lcs_shading_num+1) + filetype
    else:
        label_text = legend_labels[:]
        edge_col = "black"
        legend_file_name = "Polydotplot_LCS_Shading_Legend_%%len_n%d." % (lcs_shading_num+1) + filetype

    if prefix != None and prefix != "":
        if not prefix.endswith("-"):
            prefix = prefix + "-"
        legend_type = "LCS"
        if prefix.startswith("custom-matrix"):
            prefix = prefix.replace("custom-matrix", "")[1:]
            legend_type = "CustomMatrix"
        legend_file_name = prefix + legend_file_name.replace("LCS", legend_type)

    # plot legend figure
    fig, ax = P.subplots(3, 1, figsize=(len(colors)*2, len(colors)*2))
    for idx in range(len(colors)):
        ax[0].bar(data_points[idx]+1, data_points[idx]+1, color=colors[idx], label=label_text[idx],
               alpha=alphas[idx], edgecolor=edge_col)
        ax[1].bar(data_points[idx]+1, 0, color=colors[idx], label=label_text[idx],
               alpha=alphas[idx], edgecolor=edge_col)
        ax[2].bar(data_points[idx]+1, 0, color=colors[idx], label=label_text[idx],
               alpha=alphas[idx], edgecolor=edge_col)
    ax[1].set_ylim(0,1)
    ax[2].set_ylim(0,1)
    ax[1].legend(ncol=((len(colors)-1)//max_legend_length_row)+1, framealpha=1) # vertical legend
    col_num = len(colors)
    if len(colors) > max_legend_length_col:
        remainder = 0
        if len(colors) % max_legend_length_col != 0:
            remainder = 1
        row_num = len(colors) // max_legend_length_col + remainder
        remainder = 0
        if len(colors) % row_num != 0:
            remainder = 1
        col_num = len(colors) // row_num + remainder
    ax[2].legend(ncol=col_num, framealpha=1) # horizontal legend

    P.savefig(legend_file_name)

    return legend_file_name


###############################
#     Analysis Functions      #
###############################

def wobble_replacement(sequence, general_ambiguity_code, verbose=False):
    """
    get all degenerated sequences for sequence with ambiguous residues
    (only residues considered that are keys in wobble_dictionary)
    """

    # get positions of ambiguous residues 
    wobble_pos = []
    for idx in range(len(sequence)):
        letter = sequence[idx]
        if letter in general_ambiguity_code.keys():
            wobble_pos.append(idx)

    if verbose: 
        text = "\t%d wobbles" % len(wobble_pos)
        logprint(text, start=False, printing=True)

    # replace one wobble through each iteration by all possible residues
    # repeat if still wobbles in new kmers
    kmer_variants = [sequence]
    while True:
        if verbose:
            text = "\t\t%d kmer variants" % len(kmer_variants)
            logprint(text, start=False, printing=True)
        temp_kmers = set([])
        for kmer in kmer_variants:
            for idx in wobble_pos:
                letter = kmer[idx]
                if letter in general_ambiguity_code.keys():
                    for base in general_ambiguity_code[kmer[idx]]:
                        newkmer = kmer[:idx] + base + kmer[idx+1:]
                        temp_kmers.add(newkmer)
        wobble = False
        for kmer in temp_kmers:
            for idx in range(len(kmer)):
                letter = kmer[idx]
                if letter in general_ambiguity_code.keys():
                    wobble = True
                    break
            if wobble: 
                   break
        kmer_variants = set(list(temp_kmers)[:])
        if not wobble: 
            break

    return kmer_variants

def split_diagonals(data, stepsize=1):
    """
    split array if point difference exceeds stepsize
    data = sorted list of numbers
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return longest

def lcs_from_x_values(x_values):
    """
    calculate length of longest common substring based on nested list of numbers
    """
    if len(x_values) == 0:
        return 0
    # get lengths of each subarray data
    lengths = np.array([len(i) for i in x_values])
    return max(lengths)


###############################
#     Matching Functions      #
###############################

def find_match_pos_diag(seq1, seq2, wordsize, report_lcs=False, rc_option=True, convert_wobbles=False, max_N_percentage=49, type_nuc=True, verbose=False):
    """
    find all matching positions with matches >= wordsize
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """
    global t1 # timer

    # look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc == True:
        any_residue = "N"
    else:
        any_residue = "X"

    # read sequences 
    seq_one = seq1.upper(); len_one = len(seq_one)
    seq_two = seq2.upper(); len_two = len(seq_two)

    # set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[2] # nucleotide_ambiguity_code or aminoacid_ambiguity_code

    # forward
    #################################
    kmer_pos_dict_one = {}; kmer_pos_dict_two = {} # dictionaries for both sequences

    # reverse complement
    #################################
    kmer_pos_dict_three = {}; kmer_pos_dict_four = {} # dictionaries for both sequences

    # create dictionaries with kmers (wordsize) and there position(s) in the sequence
    if rc_option:
        data_list = [(str(seq_one), kmer_pos_dict_one),
                     (str(seq_two), kmer_pos_dict_two),
                     (str(seq_one), kmer_pos_dict_three),
                     (str(seq_two.reverse_complement()), kmer_pos_dict_four)]
    else:
        data_list = [(str(seq_one), kmer_pos_dict_one),
                     (str(seq_two), kmer_pos_dict_two)]
    for (seq, kmer_pos_dict) in data_list:
        for i in range(len(seq)-wordsize+1):
            kmer = seq[i:i+wordsize]
            # discard kmer, if too many Ns included
            if kmer.count(any_residue)*100./wordsize <= max_N_percentage:
                if not convert_wobbles:
                    try:
                        kmer_pos_dict[kmer].append(i)
                    except KeyError:
                        kmer_pos_dict[kmer] = [i]
                else:
                    wobbles = False
                    for item in general_ambiguity_code.keys():
                        if item in kmer:
                            wobbles = True
                            break
                    if not wobbles:
                        try:
                            kmer_pos_dict[kmer].append(i)
                        except KeyError:
                            kmer_pos_dict[kmer] = [i]
                    else:
                        kmer_variants = wobble_replacement(kmer, general_ambiguity_code)
                        for new_kmer in kmer_variants:
                            # print "\t", new_kmer
                            try:
                                kmer_pos_dict[new_kmer].append(i)
                            except KeyError:
                                kmer_pos_dict[new_kmer] = [i]

    # find kmers shared between both sequences
    matches_for = set(kmer_pos_dict_one).intersection(kmer_pos_dict_two) # forward
    matches_rc  = set(kmer_pos_dict_three).intersection(kmer_pos_dict_four) # reverse complement

    if verbose:
        text = "[matches: %i for; %.i rc]" % (len(matches_for), len(matches_rc))
        logprint(text, start=False, printing=True)

    # create lists of x and y co-ordinates for scatter plot
    # keep all coordinates of all shared kmers (may match multiple times)
    diag_dict_for = {}
    diag_dict_rc  = {}
    for (match_list, pos_dict1, pos_dict2, diag_dict) in [(matches_for, kmer_pos_dict_one, kmer_pos_dict_two, diag_dict_for),
                                                          (matches_rc, kmer_pos_dict_three, kmer_pos_dict_four, diag_dict_rc)]:
        for kmer in match_list:
            for i in pos_dict1[kmer]:
                for j in pos_dict2[kmer]:
                    diag = i-j
                    points = set(range(i+1, i+wordsize+1))
                    if not diag in diag_dict.keys():
                        diag_dict[diag] = points
                    else:
                        diag_dict[diag].update(points)

    # convert coordinate points to line start and stop positions
    x1 = [] # x values reverse
    y1 = [] # y values forward
    for diag in diag_dict_for.keys():
        x_values = np.array(sorted(diag_dict_for[diag]))
        x1.extend(split_diagonals(x_values))
        y_values = split_diagonals(x_values - diag)
        y1.extend(y_values)

    x2 = [] # x values rc
    y2 = [] # y values rc
    if rc_option:
        for diag in diag_dict_rc.keys():
            factor = len_two + diag + 1 
            x_values = np.array(sorted(diag_dict_rc[diag]))
            x2.extend(split_diagonals(x_values))
            y_values = split_diagonals(factor - x_values, -1)
            y2.extend(y_values)

    if verbose:
        t1 = time_track(t1)

    if not report_lcs:
        return np.array(x1), np.array(y1), np.array(x2), np.array(y2) 
    else:
        # get length of longest common substring based on match lengths
        lcs_for = lcs_from_x_values(x1)
        lcs_rev = lcs_from_x_values(x2)
        return np.array(x1), np.array(y1), np.array(x2), np.array(y2), lcs_for, lcs_rev

def find_match_pos_regex(seq1, seq2, wordsize, substitution_count=0, report_lcs=False, rc_option=True, convert_wobbles=False, max_N_percentage=49, type_nuc=True, verbose=False):
    """
    find all matching positions with matches >= wordsize via regular expression search
    fuzzy matching - allow up to substitution_count substitutions  
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """
    global t1 # timer

    # read sequences 
    seq_one = seq1.upper(); len_one = len(seq_one)
    seq_two = seq2.upper(); len_two = len(seq_two)

    # set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[2] # nucleotide_ambiguity_code or aminoacid_ambiguity_code
    ambiguity_match_dict = alphabets(type_nuc)[3] 

    ambiq_residues = "[%s]" % "".join(general_ambiguity_code.keys())

     # look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc == True:
        any_residue = "N"
    else:
        any_residue = "X"

    # check for wobble presence
    if not (regex.search(ambiq_residues, str(seq_one)) == None and regex.search(ambiq_residues, str(seq_two)) == None):
        wobble_found = True
    else:
        wobble_found = False

    # dictionary for matches
    diag_dict_for = {}
    diag_dict_rc = {}
    counter = [0, 0]

    # one-way matching
    if rc_option:
        data_list = [(str(seq_one), str(seq_two), diag_dict_for, 0),
                     (str(seq_one), str(seq_two.reverse_complement()), diag_dict_rc, 1)]
    else:
        data_list = [(str(seq_one), str(seq_two), diag_dict_for, 0)]

    for seq_query, seq_target, diag_dict, counter_pos in data_list:
        # split query sequence into kmers
        if not rc_option and counter_pos == 1:
            break 

        for idx in range(len(str(seq_query))-wordsize+1):
            kmer = str(seq_query)[idx:idx+wordsize]

            # skip excessive N/X stretches (big black areas)
            if kmer.count(any_residue)*100./wordsize <= max_N_percentage:
                #  convert kmer to regular expression for wobble_matching 
                if convert_wobbles and wobble_found:
                    kmer_string = ""
                    # replace each residue with matching residues or wobbles
                    for jdx in range(len(kmer)):
                        kmer_string += ambiguity_match_dict[kmer[jdx]]
                else:
                    kmer_string = kmer

                # convert to regular expression tolerating substitution errors
                if type(substitution_count) == int and substitution_count != 0:
                    kmer_string = "(%s){s<=%d}" % (kmer_string, substitution_count)

                # search for regular expression in target sequence
                kdx = 0
                start = True
                if regex.search(kmer_string, seq_target[kdx:]) != None:
                    counter[counter_pos] += 1
                    while regex.search(kmer_string, seq_target[kdx:]) != None:
                        # search for regular expression pattern in target sequence
                        result = regex.search(kmer_string, seq_target[kdx:])

                        kmer2 = seq_target[kdx:][result.start():result.end()]

                        # skip excessive N/X stretches (big black areas)
                        if kmer2.count(any_residue)*100./wordsize <= max_N_percentage:
                            diag = idx-(kdx+result.start())
                            points = set(range(idx+1, idx+wordsize+1))
                            if not diag in diag_dict.keys():
                                diag_dict[diag] = points
                            else:
                                diag_dict[diag].update(points)

                        kdx += result.start() + 1
                        if kdx >= len(seq_target):
                            break
                        elif regex.search(kmer_string, seq_target[kdx:]) != None:
                            counter[counter_pos] += 1

    if verbose:
        text = "%5.i \tforward matches" % counter[0]
        text += "\n%5.i \treverse complementary matches" % counter[1]
        logprint(text, start=False, printing=True)

    # convert coordinate points to line start and stop positions
    x1 = [] # x values reverse
    y1 = [] # y values forward
    for diag in diag_dict_for.keys():
        x_values = np.array(sorted(diag_dict_for[diag]))
        x1.extend(split_diagonals(x_values))
        y_values = split_diagonals(x_values - diag)
        y1.extend(y_values)

    x2 = [] # x values rc
    y2 = [] # y values rc
    if rc_option:
        for diag in diag_dict_rc.keys():
            factor = len_two + diag + 1 
            x_values = np.array(sorted(diag_dict_rc[diag]))
            x2.extend(split_diagonals(x_values))
            y_values = split_diagonals(factor - x_values, -1)
            y2.extend(y_values)

    if verbose:
        t1 = time_track(t1)

    if not report_lcs:
        return np.array(x1), np.array(y1), np.array(x2), np.array(y2) 
    else:
        # get length of longest common substring based on match lengths
        lcs_for = lcs_from_x_values(x1)
        lcs_rev = lcs_from_x_values(x2)
        return np.array(x1), np.array(y1), np.array(x2), np.array(y2), lcs_for, lcs_rev


###############################
#     Dot Plot Functions      #
###############################

def selfdotplot(input_fasta, wordsize, prefix=None, plot_size=10, label_size=10, filetype='png', type_nuc=True, convert_wobbles=False, substitution_count=0, alphabetic_sorting=False, mirror_y_axis=False, title_length=float("Inf"), title_clip_pos="B", max_N_percentage=49, verbose=False, multi=True, ncols=4, nrows=5, gff_files=[], gff_color_dict={"others": ("grey", 1, 0)}):
    """
    self-against-self dotplot
    partially from biopython cookbook
    """

    # read sequences
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logprint("\nFailed to load sequences")
        return []

    if alphabetic_sorting:
        sequences = sorted(sequences)

    # check if at least one input sequence
    if len(sequences) == 0:
        text = "\n%s\n\nCreating %s selfdotplot images\n%s\n\n=>" % (50*"=", len(sequences), 28*"-") 
        text += " No sequences provided for selfdotplot!\n\nTerminating polydotplot!"
        logprint(text, start=False, printing=True)
        return
    elif len(sequences) == 1 and multi:
        text = "\n\nCreating collage output for single selfdotplot!"
        text += "\nRecommendation: Change to individual mode by using '--collage_output n'!\n\n"
        logprint(text, start=False, printing=True)

    if multi and (ncols == 0 or nrows == 0):
        ncols = max(ncols, 1)
        nrows = max(nrows, 1)
        text = "\n\nSelfdotplot Collage: Invalid collage - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)

    if multi and ncols > len(sequences):
        ncols = len(sequences)
        nrows = 1
        text = "\n\nSelfdotplot Collage: Few sequences - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)
    elif multi and ncols*(nrows-1) > len(sequences): 
        nrows = ((len(sequences)-1) // ncols) + 1
        text = "\n\nSelfdotplot Collage: Few sequences - correcting number of rows:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)

    if multi and not (nrows == 1 and ncols == 1) and plot_size <= label_size/2:
        label_size = plot_size * 3 // 2
        text = "Reducing label size for better visualization to %d\n" % label_size
        logprint(text, start=False, printing=True)

    # read gff annotation data if provided for shading 
    if gff_files != None and gff_files != []:
        text = "\n%s\n\nReading %s GFF annotation files\n%s\n\n=> %s\n" % (50*"=", len(gff_files), 28*"-", ", ".join(gff_files))
        logprint(text, start=False, printing=True)
        if prefix != None and prefix != "":
            legend_prefix = prefix + "-Selfdotplot"
        else:  legend_prefix = "Selfdotplot"
        feat_dict = read_gffs(gff_files, color_dict=gff_color_dict, type_nuc=type_nuc, prefix=legend_prefix, filetype=filetype, verbose=verbose)

    global t1

    print "\n%s\n\nCreating %s selfdotplot images\n%s\n\n=>" % (50*"=", len(sequences), 28*"-"),
    log_txt = "\n%s\n\nCreating %s selfdotplot images\n%s\n\n=>" % (50*"=", len(sequences), 28*"-")

    # preparations for file name
    name_graph = "Selfdotplots"
    if prefix != None:
        if not prefix[-1] == "-":
            prefix = prefix + "-"
    else:
        prefix = ""
    suffix = ""
    if convert_wobbles:
        suffix += "_wobbles"
    if substitution_count != 0:
        suffix += "_S%d" % substitution_count
    if multi:
        suffix += "_collage"

    # calculate fig ratios
    if not multi:
        ncols = 1
        nrows = 1
    figsize_x, figsize_y = calc_fig_ratio(ncols, nrows, plot_size)

    P.cla() # clear any prior graph
    if multi:
        fig = P.figure(figsize=(figsize_x, figsize_y))
        page_counter = 1
    list_of_png_names = []

    counter = 0
    for seq_name in sequences:
        print seq_name,
        log_txt += " " + seq_name

        counter += 1
        if not multi:
            P.cla() # clear any prior graph

        # read sequence
        seq_record = seq_dict[seq_name]
        name_seq   = seq_record.id
        seq_one    = seq_record.seq.upper()
        length_seq = len(seq_one)

        # get positions of matches
        if substitution_count != 0:
            # print "RE"
            x_lists, y_lists, x_lists_rc, y_lists_rc = find_match_pos_regex(seq_one, seq_one, wordsize, substitution_count=substitution_count, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, type_nuc=type_nuc, verbose=verbose)
        else:
            # print "DIAG", 
            x_lists, y_lists, x_lists_rc, y_lists_rc = find_match_pos_diag(seq_one, seq_one, wordsize, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, type_nuc=type_nuc, verbose=verbose)

        # plotting with matplotlib
        #################################

        # combined plotting
        if multi:
            # plotting subplot with matplotlib
            ax = P.subplot(nrows, ncols, counter) # rows, columns, plotnumber

            # shade annotated regions
            if gff_files != None and gff_files != []:
                if seq_name in feat_dict.keys():
                    features = feat_dict[seq_name]
                    for item in features:
                        feat_type, start, stop = item
                        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                        start = max(0, start - zoom - 0.5)
                        stop  = min(length_seq+1, stop + zoom + 0.5)
                        width = stop - start
                        ax.add_patch(patches.Rectangle((start, start), # (x,y)
                                                        width, width, # width, height
                                                        edgecolor=None, linewidth=line_width+zoom,
                                                        fill=True, facecolor=feat_color, 
                                                        alpha=strength))

            # collect lines
            lines = []
            color_list = []
            for (x_lines, y_lines, col) in [(x_lists_rc, y_lists_rc, line_col_rev), (x_lists, y_lists, line_col_for)]:
                if col != "white":
                    for ldx in range(len(x_lines)):
                        lines.append([(x_lines[ldx][0], y_lines[ldx][0]), (x_lines[ldx][-1], y_lines[ldx][-1])])
                        color_list.append(col)
            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            # print P.xticks()[0], P.yticks()[0]
            P.axis('scaled') # make images quadratic
            P.xlim(0, length_seq+1)
            if mirror_y_axis:
                P.ylim(0, length_seq+1) # rotate y axis (point upwards)
            else:
                P.ylim(length_seq+1, 0) # rotate y axis (point downwards)
            P.xlabel("[%s]" % aa_bp_unit, fontsize=label_size)
            P.ylabel("[%s]" % aa_bp_unit, fontsize=label_size)
            P.tick_params(axis='both', which='major', labelsize=label_size*.9)
            
            # # use same tick labels for x and y axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            P.title(unicode_name(shorten_name(name_seq, max_len=title_length, title_clip_pos=title_clip_pos)), fontsize=label_size, fontweight='bold')
            # P.title(unicode_name(name_seq), fontsize=label_size*1.3, fontweight='bold')

            # save figure and reinitiate if page is full
            if counter == ncols * nrows:

                # finalize layout - margins & spacing between plots  
                try:
                    P.tight_layout(h_pad=.02, w_pad=.02)
                except:
                    logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
                P.subplots_adjust(hspace=0.5, wspace=0.5) # space between rows - def 0.4

                # name and create output files (names derived from SEQNAME)
                fig_name = '%s%s_wordsize%i%s-%.3d.%s' % (prefix, name_graph, wordsize, suffix, page_counter, filetype)
                P.savefig(fig_name, bbox_inches='tight')
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)

                counter = 0
                page_counter += 1

                fig = P.figure(figsize=(figsize_x, figsize_y))

        # plotting separate figure files
        else: # not multi

            fig = P.figure(figsize=(plot_size, plot_size)) # figure size needs to be a square
            ax = P.subplot(1, 1, 1) # rows, columns, plotnumber

            # shade annotated regions
            if gff_files != None and gff_files != []:
                if seq_name in feat_dict.keys():
                    features = feat_dict[seq_name]
                    for item in features:
                        feat_type, start, stop = item
                        feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                        start = max(0, start - zoom - 0.5) 
                        stop  = min(length_seq+1, stop + zoom + 0.5)
                        width = stop - start
                        ax.add_patch(patches.Rectangle((start, start), # (x,y)
                                                        width, width, # width, height
                                                        edgecolor=None, linewidth=line_width+zoom,
                                                        fill=True, facecolor=feat_color, 
                                                        alpha=strength))

            # collect lines
            lines = []
            number = 0
            color_list = []
            for (x_lines, y_lines, col) in [(x_lists_rc, y_lists_rc, line_col_rev), (x_lists, y_lists, line_col_for)]:
                if col != "white":
                    for ldx in range(len(x_lines)):
                        lines.append([(x_lines[ldx][0], y_lines[ldx][0]), (x_lines[ldx][-1], y_lines[ldx][-1])])
                        color_list.append(col)

            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            P.axis('scaled') # make images quadratic
            P.xlim(0, length_seq+1)
            if mirror_y_axis:
                P.ylim(0, length_seq+1) # rotate y axis (point upwards)
            else:
                P.ylim(length_seq+1, 0) # rotate y axis (point downwards)
            P.xlabel("[%s]" % aa_bp_unit, fontsize=label_size)
            P.ylabel("[%s]" % aa_bp_unit, fontsize=label_size)
            P.tick_params(axis='both', which='major', labelsize=label_size*.9)
            
            # # use same tick labels for x and y axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            P.title(unicode_name(shorten_name(name_seq, max_len=title_length, title_clip_pos=title_clip_pos)), fontsize=label_size*1.3, fontweight='bold')

            # name and create output files (names derived from SEQNAME)
            fig_name = '%s%s-%d_%s_wordsize%i%s.%s' %(prefix, name_graph, counter, shorten_name(name_seq, max_len=title_length, title_clip_pos=title_clip_pos), wordsize, suffix, filetype)
            P.savefig(fig_name, bbox_inches='tight')

            P.close()
            P.cla() # clear any prior graph

            list_of_png_names.append(fig_name)

    if multi and counter >= 1:
        # finalize layout - margins & spacing between plots  
        try: 
            P.tight_layout(h_pad=.02, w_pad=.02)
        except:
            logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
        P.subplots_adjust(hspace=0.5, wspace=0.5) # space between rows - def 0.4

        # name and create output files (names derived from SEQNAME)
        fig_name = '%s%s_wordsize%i%s-%.3d.%s' %(prefix, name_graph, wordsize, suffix, page_counter, filetype)
        P.savefig(fig_name, bbox_inches='tight')
        P.close()
        P.cla() # clear any prior graph

        list_of_png_names.append(fig_name)

    print "\n\nDrawing selfdotplots done"
    log_txt += "\n\nDrawing selfdotplots done"
    logprint(log_txt, start=False, printing=False)

    return list_of_png_names

def pairdotplot(input_fasta, wordsize, prefix=None, plot_size=10, label_size=10, filetype='png', type_nuc=True, convert_wobbles=False, substitution_count=0, alphabetic_sorting=False, mirror_y_axis=False, title_length=float("Inf"), title_clip_pos="B", max_N_percentage=49, verbose=False, multi=True, ncols=4, nrows=5, x_label_pos_top=True, only_vs_first_seq=False, length_scaling=True, scale_delim_col="red"):
    """
    pairwise dotplot (all-against-all)
    """

    # read sequences
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logprint("\nFailed to load sequences")
        return []

    if alphabetic_sorting:
        sequences = sorted(sequences)

    # check if at least two input sequences
    if len(sequences) < 2:
        text = "\n%s\n\nCreating %d paired dotplot image \n%s\n\n=>" % (50*"=", len(sequences)*(len(sequences)-1)/2, 36*"-") 
        text += " Please provide at least two sequences for pairdotplot!\n\nTerminating paired dotplot!"
        logprint(text, start=False, printing=True)
        return
    elif len(sequences) == 2 and multi:
        text = "\n\nCreating collage output for single pairdotplot!"
        text += "\nRecommendation: Change to individual mode by using '--collage_output n'!\n\n"
        logprint(text, start=False, printing=True)

    if multi and (ncols == 0 or nrows == 0):
        ncols = max(ncols, 1)
        nrows = max(nrows, 1)
        text = "\n\nPairdotplot Collage: Invalid collage settings - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)

    if multi and ncols > len(sequences)*(len(sequences)-1):
        ncols = len(sequences)
        nrows = 1
        text = "\n\nPairdotplot Collage: Few sequences - correcting number of rows and columns:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)
    elif multi and ncols*(nrows-1) > len(sequences)*(len(sequences)-1): 
        nrows = ((len(sequences)-1) // ncols) + 1
        text = "\n\nPairdotplot Collage: Few sequences - correcting number of rows:\n\tncols=%d, nrows=%d\n" % (ncols, nrows)
        logprint(text, start=False, printing=True)

    if not only_vs_first_seq:
        text = "\n%s\n\nCreating %d paired dotplot image for\n%s\n\n=>" % (50*"=", len(sequences)*(len(sequences)-1)/2, 36*"-") 
        text += ", ".join(sequences) + "\n"
    else:
        text = "\n%s\n\nCreating %d paired dotplot images against 1st sequence '%s':\n%s\n\n=>" % (50*"=", len(sequences)-1, sequences[0], 36*"-") 
        text += ", ".join(sequences[1:]) + "\n"
    logprint(text, start=False, printing=True)

    if multi and not (nrows == 1 and ncols == 1) and plot_size <= label_size/2:
        label_size = plot_size * 3 // 2
        text = "Reducing label size for better visualization to %d\n" % label_size
        logprint(text, start=False, printing=True)

    y_label_rotation = "vertical"
    # for cartesian coordinate system with mirrored y-axis: plot x labels below plot
    if mirror_y_axis:
        x_label_pos_top = False

    # preparations for file name
    name_graph = "Pairdotplot"
    if prefix != None:
        if not prefix[-1] == "-":
            prefix = prefix + "-"
    else:
        prefix = ""
    suffix = ""
    if convert_wobbles:
        suffix += "_wobbles"
    if substitution_count != 0:
        suffix += "_S%d" % substitution_count
    if length_scaling:
        suffix += "_scaled"
    if multi:
        suffix += "_collage"

    # calculate fig ratios
    if not multi:
        ncols = 1
        nrows = 1
    figsize_x, figsize_y = calc_fig_ratio(ncols, nrows, plot_size)

    P.cla() # clear any prior graph
    list_of_png_names = []
    if multi:
        fig = P.figure(figsize=(figsize_x, figsize_y))
        page_counter = 1

    # prepare LCS data file
    lcs_data_file = open("%sPairdotplot_wordsize%d_lcs_data_file%s.txt" % (prefix, wordsize, suffix.replace("_scaled", "").replace("_collage", "")), 'w')
    lcs_data_file.write("\t".join(["#title1", "title2", "len_seq1", "len_seq2", "len_lcs_for", "%_min_seq_len", "len_lcs_rev", "%_min_seq_len"])+"\n")

    counter, seq_counter = 0, 0
    print "Drawing pairwise dotplot...", 
    log_txt = "Drawing pairwise dotplot..." 
    if verbose:
        seq_text = ""
    for idx in range(len(sequences)-1):
        if verbose:
            print "\n%d\t%s vs." % ((seq_counter+1), sequences[idx]),
            seq_text += "\n%d\t%s vs." % ((seq_counter+1), sequences[idx])
        rec_two  = seq_dict[sequences[idx]]
        name_two = rec_two.id
        seq_two  = rec_two.seq
        len_two  = len(seq_two)

        for jdx in range(idx+1, len(sequences)):
            rec_one  = seq_dict[sequences[jdx]]
            name_one = rec_one.id
            seq_one  = rec_one.seq
            len_one  = len(seq_one)

            counter     += 1
            seq_counter += 1
            if verbose:
                print sequences[jdx], 
                seq_text += " " + sequences[jdx]
            elif not seq_counter % 25:
                print seq_counter, 
                log_txt += " " + str(seq_counter) 

            # get positions of matches
            if substitution_count != 0:
                # print "RE"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_regex(seq_one, seq_two, wordsize, substitution_count=substitution_count, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, report_lcs=True, type_nuc=type_nuc, verbose=verbose)
            else:
                # print "DIAG"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(seq_one, seq_two, wordsize, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, report_lcs=True, type_nuc=type_nuc, verbose=verbose)          

            # write LCS data file            
            lcs_data_file.write("\t".join([name_one, name_two, str(len_one), str(len_two),
                                           str(lcs_for), str(round((lcs_for*100./min(len_one, len_two)), 3)),
                                           str(lcs_rev), str(round((lcs_rev*100./min(len_one, len_two)), 3))]) + "\n")


            # plotting with matplotlib
            #################################

            # combined plotting
            if multi:
                # plotting subplot with matplotlib
                ax = P.subplot(nrows, ncols, counter) # rows, columns, plotnumber

            else:
                # calculate figure size for separate figures
                if len_one >= len_two:
                    sizing = (plot_size, max(2, (plot_size)*len_two*1./len_one))
                    # sizing = (plot_size, min(plot_size, max(2, (plot_size-2)*len_two*1./len_one+2)))
                else:
                    sizing = (max(2, (plot_size)*len_one*1./len_two), plot_size)
                    # sizing = (min(plot_size, max(2, (plot_size-2)*len_one*1./len_two+2)), plot_size)
                fig = P.figure(figsize=(plot_size, plot_size))

                ax = P.subplot(1, 1, 1)

            # collect lines
            lines = []
            color_list = []
            for (x_lines, y_lines, col) in [(x2, y2, line_col_rev), (x1, y1, line_col_for)]:
                if col != "white":
                    for ldx in range(len(x_lines)):
                        lines.append([(x_lines[ldx][0], y_lines[ldx][0]), (x_lines[ldx][-1], y_lines[ldx][-1])])
                        color_list.append(col)
            color_list = np.array(color_list)

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
            ax.add_collection(lc)

            # format axes
            P.xlabel(unicode_name(shorten_name(name_one, max_len=title_length, title_clip_pos=title_clip_pos)) + " [%s]" % aa_bp_unit, fontsize=label_size, fontweight='bold', labelpad=4)
            P.ylabel(unicode_name(shorten_name(name_two, max_len=title_length, title_clip_pos=title_clip_pos)) + " [%s]" % aa_bp_unit, fontsize=label_size, fontweight='bold', labelpad=4)
            P.tick_params(axis='both', which='major', labelsize=label_size*.9)

            # P.axis('scaled') # make images scaled by size ### optional update ###
            if not multi:
                if length_scaling:
                    ax.set_aspect(aspect='equal', adjustable='box', anchor='NW')
                P.xlim(0, len_one+1)
                # xlimit = [0, len_one+1]
                if mirror_y_axis:
                    P.ylim(0, len_two+1) # rotate y axis (point upwards)
                else:
                    P.ylim(len_two+1, 0) # rotate y axis (point downwards)
            elif not length_scaling:
                P.xlim(0, len_one+1)
                # xlimit = [0, len_one+1]
                if mirror_y_axis:
                    P.ylim(0, len_two+1) # rotate y axis (point upwards)
                else:
                    P.ylim(len_two+1, 0) # rotate y axis (point downwards)
            else:
                max_len = max(len_one, len_two)
                P.xlim(0, max_len+1)
                # xlimit = [0, max_len+1]
                if mirror_y_axis:
                    P.ylim(0, max_len+1) # rotate y axis (point upwards)
                else:
                    P.ylim(max_len+1, 0) # rotate y axis (point downwards)

                # plot line deliminating shorter sequence
                if max_len != len_one:
                    ax.plot((len_one+1, len_one+1), (0, len_two), marker="", linestyle="--", color=scale_delim_col, markerfacecolor="r") 
                if max_len != len_two:
                    ax.plot((0, len_one), (len_two+1, len_two+1), marker="", linestyle="--", color=scale_delim_col, markerfacecolor="r") 

            # # use same tick labels for x and y axis
            # if P.xlim() == P.ylim():
            #     tick_locs, tick_labels = P.yticks()
            #     P.xticks(tick_locs)
            #     P.xlim(xlimit)

            # evtl. switch x axis position
            if x_label_pos_top:
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
            P.setp(ax.get_xticklabels(), fontsize=label_size*.9)
            P.setp(ax.get_yticklabels(), fontsize=label_size*.9)

            # save figure and reinitiate if page is full
            if multi and counter == ncols * nrows:

                # finalize layout - margins & spacing between plots  
                try: 
                    P.tight_layout(h_pad=.02, w_pad=.02)
                except:
                    logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
                if x_label_pos_top:
                    P.subplots_adjust(hspace=.5, wspace=.5, top=0.95) # space between rows - def 0.4
                else:
                    P.subplots_adjust(hspace=.5, wspace=.5, bottom=0.05) # space between rows - def 0.4

                # name and create output files (names derived from SEQNAME)
                fig_name = '%s%s_wordsize%i%s-%.3d.%s' %(prefix, name_graph, wordsize, suffix, page_counter, filetype)
                P.savefig(fig_name, bbox_inches='tight')
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)

                counter = 0
                page_counter += 1

                fig = P.figure(figsize=(figsize_x, figsize_y))

            # plotting separate figure files
            elif not multi:

                # finalize layout - margins & spacing between plots  
                try:
                    P.tight_layout(h_pad=.02, w_pad=.02)
                except:
                    logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
                if y_label_rotation == "horizontal":
                    if x_label_pos_top:
                        P.subplots_adjust(hspace=0.02, wspace=0.02, left=0.13, top=0.95) # space between rows - def 0.4
                    else:
                        P.subplots_adjust(hspace=0.02, wspace=0.02, left=0.13, bottom=0.05) # space between rows - def 0.4
                else:
                    P.subplots_adjust(hspace=0.02, wspace=0.02) # space between rows - def 0.4

                # name and create output files
                fig_name = '%s%s-%d_wordsize%i%s.%s' % (prefix, name_graph, counter, wordsize, suffix, filetype)
                P.savefig(fig_name)
                P.close()
                P.cla()

                list_of_png_names.append(fig_name)
                fig = P.figure()

        if only_vs_first_seq:
            break

    # save figure
    if multi and counter >= 1:

        # finalize layout - margins & spacing between plots  
        try:
            P.tight_layout(h_pad=.02, w_pad=.02)
        except:
            logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
        if x_label_pos_top:
            P.subplots_adjust(hspace=0.5, wspace=0.5, top=0.95) # space between rows - def 0.4
        else:
            P.subplots_adjust(hspace=0.5, wspace=0.5, bottom=0.05) # space between rows - def 0.4

        # name and create output files (names derived from SEQNAME)
        fig_name = '%s%s_wordsize%i%s-%.3d.%s' %(prefix, name_graph, wordsize, suffix, page_counter, filetype)
        P.savefig(fig_name, bbox_inches='tight')
        P.close()
        P.cla()

        list_of_png_names.append(fig_name)

    if not verbose:
        print seq_counter, "done"
        log_txt += str(seq_counter) + " done"
    else:
        print "\n%d done" % seq_counter
        log_txt += "\n%d done" % seq_counter
    logprint(log_txt, start=False, printing=False)

    if verbose:
        print
        logprint(seq_text, start=False, printing=False)

    return list_of_png_names

def polydotplot(input_fasta, wordsize, prefix=None, plot_size=10, label_size=10, filetype='png', type_nuc=True, convert_wobbles=False, substitution_count=0, alphabetic_sorting=False, mirror_y_axis=False, title_length=float("Inf"), title_clip_pos="B", max_N_percentage=49, verbose=False, gff_files=[], gff_color_dict={"others": ("grey", 1, 0)}, x_label_pos_top=True, lcs_shading=True, lcs_shading_ref=0, lcs_shading_interval_len=100, lcs_shading_ori=0, lcs_shading_num=5, spacing=0.04, input_user_matrix_file="", user_matrix_print=True, rotate_labels=False):
    """
    all-against-all dotplot
    derived from dotplot function

    lcs_shading_refs:
        0 color relative to maximum lcs observed in dataset [default]
        1 color by coverage of shorter sequence (e.g. lcs = 70% of seq1)
    lcs_shading_ori
        0 forward only
        1 reverse only
        2 both orientations (in opposite plot)
    """

    # read sequences
    seq_dict, sequences = read_seq(input_fasta)
    if seq_dict == {}:
        logprint("\nFailed to load sequences")
        return []

    if alphabetic_sorting:
        sequences = sorted(sequences)

    if len(sequences) == 0:
        text = "\n%s\n\nCreating %dx%d polydotplot image\n%s\n\n=>" % (50*"=", len(sequences), len(sequences), 30*"-") 
        text += " No sequences provided for polydotplot!\n\nTerminating polydotplot!"
        logprint(text, start=False, printing=True)
        return
    elif len(sequences) == 1:
        text = "\n\nCreating polydotplot for single sequence!"
        text += "\nRecommendation: Use selfdotplot via '--plotting_mode 0'!\n\n"
        logprint(text, start=False, printing=True)

    text = "\n%s\n\nCreating %dx%d polydotplot image\n%s\n\n=>" % (50*"=", len(sequences), len(sequences), 30*"-") 
    text += " " + " ".join(sequences) + "\n"
    logprint(text, start=False, printing=True)

    # read gff annotation data if provided for shading 
    if gff_files != None and gff_files != []:
        text = "\n%s\n\nReading %s GFF annotation files\n%s\n\n=> %s\n" % (50*"=", len(gff_files), 28*"-", ", ".join(gff_files))
        logprint(text, start=False, printing=True)
        if prefix != None and prefix != "":
            legend_prefix = prefix + "-Polydotplot"
        else:  legend_prefix = "Polydotplot"
        feat_dict = read_gffs(gff_files, color_dict=gff_color_dict, type_nuc=type_nuc, prefix=legend_prefix, filetype=filetype, verbose=verbose)

    if lcs_shading and not type_nuc:
        if lcs_shading_ori != 0:
            lcs_shading_ori = 0
            text = "Protein shading does not support reverse complementary matching!\n"
            logprint(text, start=False, printing=True)

    # read custom shading matrix & match names of sequences to fasta
    if input_user_matrix_file != "" and input_user_matrix_file != None:
        logprint("Reading user matrix file: %s" % input_user_matrix_file)
        # lcs_shading_ori = 2
        custom_dict = read_matrix(input_user_matrix_file)
        if custom_dict != {}:
            custom_shading = True
            custom_similarity_dict = {}
            invalid_entries = []
            custom_max = 0
            custom_min = float("Inf")
            for key in custom_dict.keys():
                number_key = []

                # convert number into float
                try: 
                    value = float(custom_dict[key])
                    if not "." in custom_dict[key]:
                        value = int(custom_dict[key])
                    custom_max = max(custom_max, value)
                    custom_min = min(custom_min, value)
                except: 
                    value = custom_dict[key]
                    if value == "": 
                        value = None
                    invalid_entries.append(key)
                # match matrix names with sequence names
                for item in key:
                    if item in sequences:
                        number_key.append(sequences.index(item))
                    else:
                        number_key.append(-1)
                # dictionary with tuple of sorted sequence indices as key and number as value
                custom_similarity_dict[tuple(sorted(number_key))] = value
            if len(invalid_entries) != 0:
                text = "No valid number in custom similarity matrix for %d entries: \n\t" % (len(invalid_entries))
                for key in invalid_entries:
                    text += str(key) + " - " + str(custom_dict[key]) + "; "
                logprint(text[:-2]+"\n")

        text = "Custom user matrix given: min %.2f, max %.2f\n" % (custom_min, custom_max)

        # artificially rounding intervals if likely identity/divergence percentages
        if 0 <= custom_min < 1 and 0 < custom_max <= 1:
            rounding_factor = 5
            multi_factor = 100
            text += " > artificially rounding custom shading intervals: old (%.2f, %.2f) - " % (custom_min, custom_max)
            custom_min = max(0, (multi_factor*custom_min // rounding_factor) * (1.*rounding_factor/multi_factor))
            custom_max = min((multi_factor*custom_max // rounding_factor) * (1.*rounding_factor/multi_factor), 1)
            text += "new (%.2f, >%2f)\n" % (custom_min, custom_max)

        elif 0 <= custom_min < 100 and 0 < custom_max <= 100:
            rounding_factor = 5
            text += " > artificially rounding custom shading intervals: old (%.2f, %.2f) - " % (custom_min, custom_max)
            custom_min = max(0, (custom_min // rounding_factor) * rounding_factor)
            custom_max = min((custom_max // rounding_factor) * rounding_factor, 100)
            text += "new (%d, >%d)\n" % (custom_min, custom_max)

        logprint(text)

    else:
        custom_shading = False

    name_graph = "Polydotplot"
    suffix = ""
    if convert_wobbles:
        suffix += "_wobbles" 
    if substitution_count != 0:
        suffix += "_S%d" % substitution_count
    if custom_shading:
        suffix += "_matrix" 
    if lcs_shading:
        suffix += "_%dshades_ref%d_ori%s" % (lcs_shading_num+1, lcs_shading_ref, lcs_shading_ori)
        if "ref2" in suffix and type_nuc:
            suffix = suffix.replace("ref2", "%dbp" % lcs_shading_interval_len)
        elif "ref2" in suffix:
            suffix = suffix.replace("ref2", "%daa" % lcs_shading_interval_len)


    # name and create output files (names derived from SEQNAME)
    if prefix != None and str(prefix) != "":
        prefix = str(prefix) + "-"
    else:
        prefix = ""

    # preparations for background shading
    if lcs_shading or custom_shading:
        # create color range white to grey
        colors = create_color_list(lcs_shading_num+1, color_map=None, logging=True)
        colors_2 = create_color_list(lcs_shading_num+1, color_map="OrRd", logging=True)

        if custom_shading:
            text = "Custom Matrix Colors: " + ", ".join(colors_2)

    # write lcs lengths to file
    lcs_data_file = open("%sPolydotplot_lcs_data_file%s.txt" % (prefix, suffix.replace("_scaled", "").replace("_collage", "")), 'w')
    lcs_data_file.write("\t".join(["#title1", "title2", "len_seq1", "len_seq2", "len_lcs_for", "%_min_seq_len", "len_lcs_rev", "%_min_seq_len"])+"\n")

    # compare sequences pairwise - save lcs and line information in dictionary for plotting
    data_dict = {} # keys = tuple(idx, jdx), value = x1, y1, x2, y2 (line positions)
    lcs_dict  = {} # keys = tuple(idx, jdx), value = length of lcs: lcs_len or (lcs_for, lcs_rev)
    for_lcs_set = set([]) # keep lengths to calculate max (excluding self comparisons)
    rev_lcs_set = set([]) # keep lengths to calculate max (all)

    text = "\nTotal plot count:   %d" % (len(sequences)*(len(sequences)))
    text += "\nTotal calculations: %d" % (len(sequences)*(len(sequences)+1)/2)
    logprint(text, start=False, printing=True)

    print "\nCalculating shared regions and lengths of longest_common_substring...", 
    log_txt = "\nCalculating shared regions and lengths of longest_common_substring..."
    # determine  matches and length of lcs by comparing all sequence pairs
    if verbose:
        seq_text = ""
    counter = 0
    for idx in range(len(sequences)):
        if verbose:
            print "\n%d\t%s vs." % ((counter+1), sequences[idx]),
            seq_text += "\n%d\t%s vs." % ((counter+1), sequences[idx])
        rec_two  = seq_dict[sequences[idx]]
        name_two = rec_two.id
        seq_two  = rec_two.seq
        len_two  = len(seq_two)

        for jdx in range(idx, len(sequences)):
            rec_one  = seq_dict[sequences[jdx]]
            name_one = rec_one.id
            seq_one  = rec_one.seq
            len_one  = len(seq_one)

            counter += 1
            if verbose:
                print sequences[jdx], 
                seq_text += " " + sequences[jdx]
            elif len(sequences) < 5:
                print "\t%s (%d %s), %s (%d %s)" % (name_one, len_one, aa_bp_unit, name_two, len_two, aa_bp_unit)
                log_txt += "\t%s (%d %s), %s (%d %s)\n" % (name_one, len_one, aa_bp_unit, name_two, len_two, aa_bp_unit)
            else:
                if not counter % 25:
                    print counter,
                    log_txt += str(counter)

            # get positions of matches &  length of longest common substring based on match lengths
            if substitution_count != 0:
                # print "RE"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_regex(seq_one, seq_two, wordsize, substitution_count=substitution_count, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, report_lcs=True, type_nuc=type_nuc, verbose=verbose)
            else:
                # print "DIAG"
                x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(seq_one, seq_two, wordsize, convert_wobbles=convert_wobbles, max_N_percentage=max_N_percentage, report_lcs=True, type_nuc=type_nuc, verbose=verbose) 
            data_dict[(idx, jdx)] = x1[:], y1[:], x2[:], y2[:]
            lcs_dict[idx, jdx] = lcs_for, lcs_rev

            if idx != jdx:
                for_lcs_set.add(lcs_for)
            rev_lcs_set.add(lcs_rev)

            lcs_data_file.write("\t".join([name_one, name_two, str(len_one), str(len_two),
                                           str(lcs_for), str(round((lcs_for*100./min(len_one, len_two)), 3)),
                                           str(lcs_rev), str(round((lcs_rev*100./min(len_one, len_two)), 3))]) + "\n")

    if not verbose:
        print len(sequences)*(len(sequences)+1)/2, " done\n"
        log_txt += str(len(sequences)*(len(sequences)+1)/2) + " done\n"
    else:
        print "\n%d done" % (len(sequences)*(len(sequences)+1)/2)
        log_txt += "\n%d done" % (len(sequences)*(len(sequences)+1)/2)
    logprint(log_txt, start=False, printing=False)

    if verbose:
        logprint ("\n\nlcs_dict\n" + str(lcs_dict))
        if custom_shading:
            logprint ("\ncustom_dict\n" + str(custom_dict))
            logprint ("\ncustom_similarity_dict\n\n" + str(custom_similarity_dict))

    if verbose:
        print
        logprint(seq_text+"\n", start=False, printing=False)

    if lcs_shading_ref == 2:
        color_bins = []
        text = "\nLCS lengh bins: "
        for idx in range(lcs_shading_num):
            color_bins.append(lcs_shading_interval_len*(idx+1))
            text += " " + str(lcs_shading_interval_len*(idx+1))
        logprint(text, start=False, printing=True)

    # calculate maximum lcs length
    if lcs_shading_ori == 0: # forward only
        if len(for_lcs_set) != 0:
            max_lcs = max(for_lcs_set)
        else:
            max_lcs = None
    elif lcs_shading_ori == 1: # reverse complement only
        if len(rev_lcs_set) != 0:
            max_lcs = max(rev_lcs_set)
        else:
            max_lcs = None
    else: # both orientations
        if len(rev_lcs_set) != 0 and len(for_lcs_set) != 0:
            max_lcs = max(max(rev_lcs_set), max(for_lcs_set))
        elif len(rev_lcs_set) != 0:
            max_lcs = max(rev_lcs_set)
        elif len(for_lcs_set) != 0:
            max_lcs = max(for_lcs_set)
        else:
            max_lcs = None

    if not max_lcs == None:
        text =   "Maximum LCS:          %d %s" % (max_lcs, aa_bp_unit)
        logprint(text, start=False, printing=True)
    if custom_shading:
        text = "Maximum custom value: %d\n" % custom_max
        logprint(text, start=False, printing=True)

    # count sequences
    ncols = len(sequences); nrows = len(sequences)

    # get sequence lengths to scale plot widths and heights accordingly
    size_ratios = []
    for item in sequences:
        size_ratios.append(len(seq_dict[item].seq))  

    P.cla() # clear any prior graph
    # use GridSpec to resize plots according to sequence length
    if mirror_y_axis:
        height_ratios = size_ratios[::-1]
    else:
        height_ratios = size_ratios[:]
    gs = gridspec.GridSpec(nrows, ncols,
                           width_ratios=size_ratios,
                           height_ratios=height_ratios)
    fig = P.figure(figsize=(plot_size, plot_size))

    # for cartesian coordinate system with mirrored y-axis: plot x labels below plot
    if mirror_y_axis and representation == 1:
        x_label_pos_top = True
    elif mirror_y_axis or representation == 2:
        x_label_pos_top = False

    # print y labels on the right, if upper right triangle is displayed
    if (representation == 1 and not mirror_y_axis) or (representation == 2 and mirror_y_axis): 
        y_label_pos = 0 # last column
    else: # left y label 
        y_label_pos = 1 # first column

    # determine label orientations
    if len(sequences) > 5 or rotate_labels:
        x_label_rotation = 45
        y_label_rotation = "horizontal"
        if x_label_pos_top:
            xhalign = 'left'
            xvalign = 'bottom'
        else:
            xhalign = 'right'
            xvalign = 'top'
        yhalign = "right"
    else:
        x_label_rotation = "horizontal"
        y_label_rotation = "vertical"
        xvalign = "center"
        xhalign = "center"
        yhalign = "center"
    yvalign = 'center'

    # check combination of shading parameters for triangular output
    if representation != 0 and lcs_shading and custom_shading: # both directions in triangle
        logprint("\nAttention: For triangular output custom-shading and LCS shading cannot be combined!\n")
    elif representation != 0 and lcs_shading and lcs_shading_ori == 2: # both directions in triangle
        logprint("\nAttention: For triangular output LCS shading for both orientations is combined to max of both orientations!\n")

    print "\nDrawing polydotplot...",
    log_txt = "\nDrawing polydotplot..."

    # draw subplots
    if verbose:
        if lcs_shading and custom_shading:
            lcs_text = "\n" + "\t".join(["#Seq1", "Seq2", "LCS for [%s]" %aa_bp_unit, "LCS for [%s]" %aa_bp_unit, "Custom matrix value", "Matrix color index", "LCS color index"]) + "\n"
        elif lcs_shading:
            lcs_text = "\n" + "\t".join(["#Seq1", "Seq2", "LCS for [%s]" %aa_bp_unit, "LCS for [%s]" %aa_bp_unit, "LCS color index for", "LCS color index rev"]) + "\n"
        elif custom_shading:
            lcs_text = "\n" + "\t".join(["#Seq1", "Seq2", "Custom matrix value", "Color index for", "Color index rev"]) + "\n"

    if verbose:
        seq_text = ""
    counter, seq_counter = 0, 0
    for idx in range(len(sequences)):
        if verbose:
            print "\n%d\t%s vs." % ((seq_counter+1), sequences[idx]), 
            seq_text += "\n%d\t%s vs." % ((seq_counter+1), sequences[idx])
        rec_two = seq_dict[sequences[idx]]
        len_two = len(rec_two.seq)
        name_two = rec_two.id

        for jdx in range(idx, len(sequences)):
            rec_one = seq_dict[sequences[jdx]]
            len_one = len(rec_one.seq)
            name_one = rec_one.id

            counter += 1
            seq_counter += 1
            if verbose:
                print sequences[jdx], 
                seq_text += " " + sequences[jdx]
            elif not seq_counter % 25:
                print seq_counter,
                log_txt += str(seq_counter)

            # optional shade background according to length of LCS and/or user matrix
            #########################################################################

            # get interval based on LCS
            background_colors = [None, None]
            if lcs_shading and (lcs_shading_ref==1 or lcs_shading_ref==2 or max_lcs!=None): # self plot max_lcs_for == None
                lcs_len = lcs_dict[(idx, jdx)]
                l1 = lcs_len[0] # forward
                l2 = lcs_len[1] # reverse complement

                lcs_shading_bool = True

                # calculate shading acc. to chosen option
                if lcs_shading_ref == 1: # percentage of shorter sequence
                    color_idx0 = min(len(colors)-1, l1*lcs_shading_num // min(len_one, len_two))
                    color_idx1 = min(len(colors)-1, l2*lcs_shading_num // min(len_one, len_two))
                elif lcs_shading_ref == 2: # by given interval size
                    color_idx0 = min(len(colors)-1, l1 // lcs_shading_interval_len)
                    color_idx1 = min(len(colors)-1, l2 // lcs_shading_interval_len)
                    if color_idx0 >= len(colors):
                        color_idx0 = len(colors)
                    if color_idx1 >= len(colors):
                        color_idx1 = len(colors)
                else: # percentage of maximum lcs length
                    color_idx0 = min(len(colors)-1, l1*lcs_shading_num // max_lcs)
                    color_idx1 = min(len(colors)-1, l2*lcs_shading_num // max_lcs)
            else:
                lcs_shading_bool = False

            # get interval based on custom matrix
            if custom_shading:
                # matrix value
                try:
                    custom_value = custom_similarity_dict[(idx, jdx)]
                except:
                    custom_value = ""

                # bottom left triangle = LCS forward/reverse or best of both
                if lcs_shading_bool:
                    if lcs_shading_ori == 0: # forward
                        color_idx1 = color_idx0
                    elif lcs_shading_ori == 2: # both directions
                        color_idx1 = max(color_idx0, color_idx1)

                # top right triangle = custom value (not colored if text matrix provided)
                if type(custom_value) == int or type(custom_value) == float: 
                    color_idx0 = int((custom_value-custom_min)*lcs_shading_num // (custom_max-custom_min))
                # no color if string is proviced
                else:
                    color_idx0 = 0

            # use best LCS of both orientations for coloring triangle with two-ori-LCS
            if representation != 0 and lcs_shading_ori == 2: # both directions in triangle
                color_idx0, color_idx1 = max(color_idx0, color_idx1), max(color_idx0, color_idx1) 

            # set colors dependent on lcs dependent on orientation
            if lcs_shading_bool and not custom_shading:
                if idx != jdx:   
                    if lcs_shading_ori == 0:
                        color_idx1 = color_idx0
                    elif lcs_shading_ori == 1:
                        color_idx0 = color_idx1
                    background_colors[0] = colors[color_idx0] 
                    background_colors[1] = colors[color_idx1] 
                # for selfcomparison, only color reverse complement
                elif lcs_shading_ori != 0 and not custom_shading:
                    background_colors[0] = colors[color_idx1] 
            # set different colors for shading by LCS + user matrix
            elif lcs_shading_bool and custom_shading:
                # print colors, background_colors, color_idx0, color_idx1
                background_colors[0] = colors_2[color_idx0] 
                background_colors[1] = colors[color_idx1]
            # set grey color range for user matrix if no LCS shading
            elif custom_shading: 
                background_colors[0] = colors[color_idx0] 
                background_colors[1] = colors[color_idx0] 

            if verbose:
                if custom_shading and lcs_shading_bool:
                    lcs_text += "\t".join([name_one, name_two, str(lcs_len[0]), str(lcs_len[1]), str(custom_value), str(color_idx0), str(color_idx1)]) + "\n"
                elif lcs_shading_bool:
                    lcs_text += "\t".join([name_one, name_two, str(lcs_len[0]), str(lcs_len[1]), str(color_idx0), str(color_idx1)]) + "\n"
                elif custom_shading:
                    lcs_text += "\t".join([name_one, name_two, str(custom_value), str(color_idx0), str(color_idx1)]) + "\n"

            # calculate figure position in polyplot
            # diagonal (self-dotplots)
            if idx == jdx:
                if mirror_y_axis:
                    seq_num = sequences.index(name_one)+1
                    counter1 = seq_num + len(sequences) * (len(sequences)-seq_num)
                    counter  = counter + (counter - 1) // (nrows)
                else:
                    # skip positions below diagonal
                    counter1 = counter + (counter - 1) // (nrows) # + row_pos
                    counter = counter1
                counters = [counter1]

            # draw both graphs at once (due to symmetry)
            else:
                if mirror_y_axis:
                    col_pos = sequences.index(name_two)+1
                    row_pos = len(sequences) - (sequences.index(name_one)+1)
                    counter1 = row_pos * ncols + col_pos 
                    counter2 = (ncols - col_pos) * ncols + ncols - row_pos
                else:
                    counter1 = counter
                    col_pos = (counter - 1) % ncols
                    row_pos = (counter - 1) // (nrows)
                    counter2 = col_pos * ncols + row_pos + 1
                counters = [counter1, counter2] # lower, upper

            if len(counters) == 2:
                seq_counter += 1
                if not verbose and not seq_counter % 25:
                    print seq_counter,
                    log_txt += str(seq_counter)

            x_lists, y_lists, x_lists_rc, y_lists_rc = data_dict[(idx, jdx)]

            # plot diagram(s)
            for kdx in range(len(counters)):

                if representation == 0 or len(counters) == 1 or (representation == 1 and kdx == 0) or (representation == 2 and kdx == 1):

                    fig_pos = counters[kdx]
                    # plotting subplot with matplotlib
                    ax = P.subplot(gs[fig_pos-1]) # rows, columns, plotnumber

                    # shade annotated regions if gff file(s) provided
                    if idx == jdx and gff_files != None and gff_files != []:
                        if name_one in feat_dict.keys():
                            features = feat_dict[name_one]
                            if len_two != len_one:
                                logprint("Polydot GFF shading for diagonal fields - nequal length error!")
                                return
                            for item in features:
                                feat_type, start, stop = item
                                feat_color, strength, zoom = gff_color_dict[feat_type.lower()]
                                start = max(0, start - zoom - 0.5)
                                stop  = min(len_one+1, stop + zoom + 0.5)
                                width = stop - start
                                ax.add_patch(patches.Rectangle((start, start), # (x,y)
                                                                width, width, # width, height
                                                                edgecolor=None, linewidth=line_width+zoom,
                                                                fill=True, facecolor=feat_color, 
                                                                alpha=strength))

                    # if custom matrix value printed into upper matrix triangle, skip data plotting
                    # text print in top triangle
                    if user_matrix_print and custom_shading and kdx==0 and idx!=jdx: 
                        data_plotting = False
                    # dotplot in bottom triangle
                    else: 
                        data_plotting = True

                    # mirror plot, if plotting below diagonal
                    if kdx == 0:
                        l1, l2 = len_one, len_two
                        n1, n2 = name_one, name_two
                        x1, y1 = x_lists, y_lists
                        x2, y2 = x_lists_rc, y_lists_rc
                    else:
                        l2, l1 = len_one, len_two
                        n2, n1 = name_one, name_two
                        x1, y1 = y_lists, x_lists
                        x2, y2 = y_lists_rc, x_lists_rc

                    if mirror_y_axis:
                        x1, y1, x2, y2 = y1, x1, y2, x2
                        n1, n2 = n2, n1

                    if data_plotting:
                        # collect lines
                        lines = []
                        color_list = []
                        for (x_lines, y_lines, col) in [(x2, y2, line_col_rev), (x1, y1, line_col_for)]:
                            if col != "white":
                                for ldx in range(len(x_lines)):
                                    lines.append([(x_lines[ldx][0], y_lines[ldx][0]), (x_lines[ldx][-1], y_lines[ldx][-1])])
                                    color_list.append(col)
                        color_list = np.array(color_list)

                        # draw lines
                        lc = cllct.LineCollection(lines, colors=color_list, linewidths=line_width)
                        ax.add_collection(lc) 

                    # plot value provided by customer instead of dotplot
                    else: 
                        alignment = {'horizontalalignment': 'center', 'verticalalignment': 'center'}
                        # P.text(0.5, 0.5, custom_value, size='medium', transform=ax.transAxes, **alignment)
                        P.text(0.5, 0.5, custom_value, size=label_size*1.5, transform=ax.transAxes, **alignment)
                        # P.text(0.5, 0.5, custom_value, size=label_size*1.5, transform=ax.transAxes, 
                               # horizontalalignment='center', verticalalignment='center', color="black")

                    if custom_shading:
                        # omit diagonal
                        if idx == jdx:
                            ax.set_facecolor("white")
                        # use white background for text fields (top right triangle only [kdx 0])
                        elif type(custom_value) != int and type(custom_value) != float and kdx == 0:
                            ax.set_facecolor("white")
                        else:
                            ax.set_facecolor(background_colors[kdx])
                    # set background color if lcs shading
                    elif lcs_shading_bool and background_colors[kdx] != None:
                        ax.set_facecolor(background_colors[kdx])

                    # set axis limits
                    # P.xlim(0, l1+1)
                    if mirror_y_axis:
                        P.xlim(0, l2+1)
                        P.ylim(0, l1+1) # rotate y axis (point upwards)
                    else:
                        P.xlim(0, l1+1)
                        P.ylim(l2+1, 0) # rotate y axis (point downwards)

                    ## axis labelling
                    ##################

                    # determine axis positions
                    if x_label_pos_top:
                        ax.xaxis.tick_top()
                        ax.xaxis.set_label_position('top')
                        x_label_bool = fig_pos <= ncols
                        x_tick_bool  = fig_pos > ncols*(ncols-1)
                    else:
                        x_label_bool = fig_pos > ncols*(ncols-1)
                        x_tick_bool  = fig_pos <= ncols

                    # settings for y labels on right side
                    if y_label_pos == 0: # right label
                        ax.yaxis.tick_right()
                        ax.yaxis.set_label_position("right")
                        label_dist = 30
                    else:
                        label_dist = 8

                    # x axis labels dependent on plot position/number 
                    if x_label_bool: # x title and labels on top or bottom
                        P.xlabel(unicode_name(shorten_name(n1, max_len=title_length, title_clip_pos=title_clip_pos)), fontsize=label_size, rotation=x_label_rotation, verticalalignment=xvalign, horizontalalignment=xhalign, fontweight='bold', labelpad=8) # axis naming
                        if not x_label_rotation in ["horizontal", "vertical"]:
                            P.setp(ax.get_xticklabels(), fontsize=label_size*.9, rotation="vertical")
                        else:
                            P.setp(ax.get_xticklabels(), fontsize=label_size*.9, rotation=x_label_rotation)
                    elif x_tick_bool and x_label_pos_top: # x ticks on bottom row
                        ax.xaxis.tick_bottom() # ticks without labels on bottom
                        P.setp(ax.get_xticklabels(), fontsize=label_size, rotation=x_label_rotation, visible=False)
                    elif x_tick_bool: # x ticks on top row
                        ax.xaxis.tick_top() # # ticks without labels on top
                        P.setp(ax.get_xticklabels(), fontsize=label_size, rotation=x_label_rotation, visible=False) # inner diagrams without labelling
                    elif idx == jdx and representation != 0:
                        if not mirror_y_axis and representation == 1:  # upper
                            ax.xaxis.tick_bottom()
                        elif mirror_y_axis and representation == 2:  # lower
                            ax.xaxis.tick_top()
                        elif mirror_y_axis and representation == 1:  # upper
                            ax.xaxis.tick_bottom()
                        elif not mirror_y_axis and representation == 2:  # lower
                            ax.xaxis.tick_top()
                        P.setp(ax.get_xticklabels(), visible=False) # inner diagrams without labelling
                    else: # no x ticks on internal rows
                        ax.axes.get_xaxis().set_visible(False)

                    # y axis labels dependent on plot position/number 
                    if fig_pos % ncols == y_label_pos or (ncols == 1 and nrows == 1): # y title and labels in 1st column
                        P.ylabel(unicode_name(shorten_name(n2, max_len=title_length, title_clip_pos=title_clip_pos)), fontsize=label_size, rotation=y_label_rotation, verticalalignment=yvalign, horizontalalignment=yhalign, fontweight='bold', labelpad=label_dist)
                        P.setp(ax.get_yticklabels(), fontsize=label_size*.9) # axis naming
                    elif fig_pos % ncols == 0: # y ticks in last column
                        ax.yaxis.tick_right()
                        P.setp(ax.get_yticklabels(), visible=False) # inner diagrams without labelling
                    elif idx == jdx and representation != 0:
                        if not mirror_y_axis and representation == 1:  # upper
                            ax.yaxis.tick_left()
                        elif mirror_y_axis and representation == 2:  # lower
                            ax.yaxis.tick_left()
                        elif mirror_y_axis and representation == 1:  # upper
                            ax.yaxis.tick_right()
                        elif not mirror_y_axis and representation == 2:  # lower
                            ax.yaxis.tick_right()
                        P.setp(ax.get_yticklabels(), visible=False) # inner diagrams without labelling
                    else:
                        ax.axes.get_yaxis().set_visible(False)

    if not verbose:
        print seq_counter, "done"
        log_txt += str(seq_counter) + " done"
    else:
        print "\n%d done" % seq_counter
        log_txt += "\n%d done" % seq_counter
    logprint(log_txt, start=False, printing=False)

    if verbose:
        try:
            logprint(lcs_text, start=False, printing=True)
        except:
            pass

    # finalize layout - margins & spacing between plots  
    P.tick_params(axis='both', which='major', labelsize=label_size*.9)
    try:
        P.tight_layout(h_pad=.02, w_pad=.02)
    except:
        logprint("Attention - pylab.tight_layout failed! Please check sequence names and layout settings!")
    # gs.tight_layout(fig, h_pad=.02, w_pad=.02) # less overlapping tick labels, but also disturbingly large spacing
    if y_label_rotation == "horizontal":
        if x_label_pos_top:
            P.subplots_adjust(hspace=spacing, wspace=spacing, left=0.13, top=0.87) # space between rows - def 0.4
        else:
            P.subplots_adjust(hspace=spacing, wspace=spacing, left=0.13, bottom=0.13) # space between rows - def 0.4
    else:
        P.subplots_adjust(hspace=spacing, wspace=spacing) # space between rows - def 0.4

    # save figure and close instance
    fig_name = '%s%s_wordsize%i%s.%s' % (prefix, name_graph, wordsize, suffix, filetype)
    P.savefig(fig_name)
    P.close()
    P.cla()


    # create figure color legend
    if lcs_shading:
        if lcs_shading_ref == 1: # percentage of shorter sequence
            legend_file_name = legend_figure(colors, lcs_shading_num, unit="%", filetype=filetype, prefix=prefix)
        elif lcs_shading_ref == 2: # interval sizes
            legend_file_name = legend_figure(colors, lcs_shading_num, unit=aa_bp_unit, filetype=filetype, prefix=prefix, bins=color_bins)
        else: # relative of maximum lcs
            legend_file_name = legend_figure(colors, lcs_shading_num, unit=aa_bp_unit, filetype=filetype, prefix=prefix, max_len=max_lcs)

    if custom_shading:
        custom_prefix = "custom-matrix-" + prefix
        legend_file_name_custom = legend_figure(colors_2, lcs_shading_num, unit="%", filetype=filetype, prefix=custom_prefix, max_len=custom_max, min_len=custom_min)

    if lcs_shading and custom_shading:
        return [fig_name, legend_file_name, legend_file_name_custom]
    elif lcs_shading:
        return [fig_name, legend_file_name]
    elif custom_shading:
        return [fig_name, legend_file_name_custom]
    else:
        return [fig_name]


###############################
#        Function Call        #
###############################

def main(seq_list, wordsize, modes=[0, 1, 2], prefix=None, plot_size=10, label_size=10, filetype="png", type_nuc=True, convert_wobbles=False, substitution_count=0, rc_option=True, alphabetic_sorting=False, only_vs_first_seq=False, gff=None, multi=True, ncols=1, nrows=1, lcs_shading=True, lcs_shading_num=5, lcs_shading_ref=0, lcs_shading_interval_len=100, lcs_shading_ori=0, gff_color_config_file="", input_user_matrix_file="", user_matrix_print=False, length_scaling=True, title_length=50, title_clip_pos="B", spacing=0.04, max_N_percentage=49, mirror_y_axis=False, verbose=False):

    global t1, line_col_rev

    # check input variables
    if convert_wobbles and max_N_percentage > 49:
        max_N_percentage = 49
        if type_nuc:
            ambiq_res = "N"
        else:
            ambiq_res = "X"
        text = "Provide valid max_N_percentage, kmers with >50%% %ss are ignored\n" % (ambiq_res)
        logprint(text, start=False, printing=True)

    if filetype not in ["png", "pdf", "svg"]:
        text = "Provide valid file type - png, pdf, or svg - given:%s\n" % filetype
        logprint(text, start=False, printing=True)
        filetype = "png"

    # read gff color config file if provided
    if len(input_gff_files) != 0 and input_gff_files != None:
        if gff_color_config_file not in ["", None]:
            text = "\n%s\n\nReading GFF color configuration file\n%s\n\n=> %s\n" % (50*"=", 28*"-", gff_color_config_file)
            logprint(text, start=False, printing=True)
        gff_feat_colors = read_gff_color_config(gff_color_config_file)
    else:
        gff_feat_colors = {}
        if gff_color_config_file not in ["", None]:
            text = "Please provide GFF annotation files to use configuration file", gff_color_config_file 
            logprint(text, start=False, printing=True)

    # if color is set to white, reverse complementary matches are skipped        
    if not rc_option:
        line_col_rev = "white" # reverse matches not calculated
    elif not type_nuc:
        logprint("Reverse complement deactivated for proteins!")
        line_col_rev = "white" # reverse matches not calculated

    mode_text = []
    for item in modes:
        mode_text.append(str(item))
    text = "%s\n\nRunning plotting modes %s" % (50*"=", ", ".join(mode_text))
    logprint(text, start=False, printing=True)


    # create dotplots
    ##########################################

    # self dotplots
    t1 = time.time()
    if 0 in modes:
        list_of_png_names = selfdotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, multi=multi, ncols=ncols, nrows=nrows, gff_files=gff, gff_color_dict=gff_feat_colors, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
        t1 = time_track(t1)
        if list_of_png_names != [] and list_of_png_names != None:
            text = "-> Image file(s): %s\n" % ", ".join(list_of_png_names)
        else:
            text = "No image files were created!\n"
        logprint(text, start=False, printing=True)
        logprint(50*"=")

    # paired dotplots
    if 1 in modes:
        if multi:
            list_of_png_names = pairdotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, only_vs_first_seq=only_vs_first_seq, multi=multi, ncols=ncols, nrows=nrows, length_scaling=length_scaling, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
            t1 = time_track(t1)
        else:
            if not length_scaling:
                    text = "\nPairwise dotplot with individual output files scaled by sequence length automatically!"
                    logprint(text, start=False, printing=True)
            list_of_png_names = pairdotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, only_vs_first_seq=only_vs_first_seq, multi=multi, ncols=ncols, nrows=nrows, length_scaling=True, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
            t1 = time_track(t1)
        if list_of_png_names != [] and list_of_png_names != None:
            text = "-> Image file(s): %s\n" % ", ".join(list_of_png_names)
        else:
            text = "No image files were created!\n"
        logprint(text, start=False, printing=True)
        logprint(50*"=")

    # all-against-all dotplot
    if 2 in modes:
        list_of_png_names = polydotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, lcs_shading=lcs_shading, lcs_shading_num=lcs_shading_num, lcs_shading_ref=lcs_shading_ref, lcs_shading_interval_len=lcs_shading_interval_len, lcs_shading_ori=lcs_shading_ori, input_user_matrix_file=input_user_matrix_file, user_matrix_print=user_matrix_print, spacing=spacing, gff_files=gff, gff_color_dict=gff_feat_colors, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
        t1 = time_track(t1)
        if list_of_png_names != [] and list_of_png_names != None:
            text = "-> Image file(s): %s\n" % ", ".join(list_of_png_names)
        else:
            text = "No image files were created!\n"
        logprint(text, start=False, printing=True)
        logprint(50*"=")

    text = "\n" + 50 * "#" + "\n" + 50 * "#"
    text += "\n\nThank you for using FlexiDot!\n"
    logprint(text, start=False, printing=True)


load_modules()

# testing mode for debugging
trial_mode = True
trial_mode = False

# parameters = check_input(sys.argv)
parameters = check_input(sys.argv, trial_mode=trial_mode)

# read out parameters
commandline, auto_fas, input_fasta, output_file_prefix, collage_output, m_col, n_row, filetype, type_nuc, input_gff_files, gff_color_config_file, wordsize, plotting_modes, wobble_conversion, substitution_count, rc_option, alphabetic_sorting, only_vs_first_seq, lcs_shading, lcs_shading_num, lcs_shading_ref, lcs_shading_interval_len, lcs_shading_ori, input_user_matrix_file, user_matrix_print, plot_size, line_width, line_col_for, line_col_rev, x_label_pos_top, label_size, spacing, length_scaling, title_length, title_clip_pos, max_N_percentage, mirror_y_axis, representation, verbose = parameters

# evtl. overwrite parameters for testing purposes in trial mode
if trial_mode:
    input_fasta        = ["Inversionen_IDs_v2_test2.fas"]
    # input_fasta        = ["Inversionen_IDs_v2_test3.fas"]
    # input_fasta        = ["test-sequences-8.fas"]
    # input_gff_files    = ["Seq2_annotations.gff3"]
    # input_user_matrix_file = "matrix.txt"
    # user_matrix_print  = True
    output_file_prefix = "#Test"
    plot_size = 10
    plotting_modes  = [0,1,2]
    plotting_modes  = [2]
    plotting_modes  = [0]
    lcs_shading     = False
    lcs_shading     = True
    lcs_shading_ref = 2
    lcs_shading_num = 4
    lcs_shading_ori = 0
    lcs_shading_interval_len = 15
    wordsize        = 10
    wordsize        = 7
    x_label_pos_top = True
    filetype        = "pdf"
    filetype        = "png"
    mirror_y_axis   = False
    mirror_y_axis   = True

    output_file_prefix = "#R-upper"
    representation  = 0 # both 
    representation  = 1 # upper 
    representation  = 2 # lower 

    wobble_conversion = False
    wobble_conversion = True

    substitution_count = 0
    substitution_count = 1

    rc_option     = True
    rc_option     = False
    label_size    = 10

    verbose = False
    verbose = True

if auto_fas:
    path = os.path.dirname(os.path.abspath(__file__))
    files_long = glob.glob(path+"/*.fasta")
    files_long.extend(glob.glob(path+"/*.fas"))
    files_long.extend(glob.glob(path+"/*.fa"))
    files_long.extend(glob.glob(path+"/*.fna"))
    input_fasta = []
    for i in files_long:
        if not "combined" in i:
            filename = i[i.rfind('\\')+1:]
            input_fasta.append(filename)

if trial_mode:
    # start logging file
    logprint(commandline, start=True, printing=False, prefix=output_file_prefix)



######################
# FlexiDot Execution #
######################

main(input_fasta, wordsize, modes=plotting_modes, prefix=output_file_prefix, plot_size=plot_size, label_size=label_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=wobble_conversion, substitution_count=substitution_count, rc_option=rc_option, alphabetic_sorting=alphabetic_sorting, only_vs_first_seq=only_vs_first_seq, gff=input_gff_files, multi=collage_output, ncols=m_col, nrows=n_row, lcs_shading=lcs_shading, lcs_shading_num=lcs_shading_num, lcs_shading_ref=lcs_shading_ref, lcs_shading_interval_len=lcs_shading_interval_len, lcs_shading_ori=lcs_shading_ori, gff_color_config_file=gff_color_config_file, input_user_matrix_file=input_user_matrix_file, user_matrix_print=user_matrix_print, length_scaling=length_scaling, title_length=title_length, title_clip_pos=title_clip_pos, spacing=spacing, max_N_percentage=max_N_percentage, mirror_y_axis=mirror_y_axis, verbose=verbose)


