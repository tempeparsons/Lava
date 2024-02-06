import os, sys
import numpy as np
import pandas as pd

from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.backends.backend_pdf import PdfPages

import lava_util as util

# Allow shirt codes for software
SOFT_DN = ('DIANN', 'DN')
SOFT_MQ = ('MaxQuant', 'MQ')
SOFT_PD = ('ProteomeDiscoverer', 'PD')
 
VALID_SOFTWARE = SOFT_DN + SOFT_MQ + SOFT_PD
DEFAULT_SOFTWARE = SOFT_PD[0]

FILE_EXTS_EXL = ('.xls','.xlsx')
FILE_EXTS_TXT = ('.txt',)
FILE_EXTS_CSV = ('.csv','.tsv')

DEFALUT_LOGFC = 2.0
DEFALUT_MAXP  = 4.32
DEFAULT_POS_COLOR = '#0080FF'
DEFAULT_NEG_COLOR = '#FF0000'
DEFAULT_LOW_COLOR = '#FFD0D0'

DEFAULT_INDEX_COLS = ('Accession','Protein.ID','Protein ID')

def warn(txt):
   
   _report('WARNING: ' + txt)


def info(txt):
   
   _report('INFO: ' + txt)


def fail(txt):
   
   _report('FAILURE: ' + txt)
   _report('EXIT')
   sys.exit(0)

def _report(txt):

    print(txt)


def _read_data(file_path):
    
    print(f'Reading {file_path}')
    
    path_root, file_ext = os.path.splitext(file_path)
    file_ext = file_ext.lower()
    
    if file_ext in FILE_EXTS_EXL:
        df = pd.read_excel(file_path, engine='openpyxl')
        
    elif file_ext in FILE_EXTS_CSV:
        df = pd.read_csv(file_path)
        
    else:    # Last resort try plain text
        df = pd.read_table(file_path)
  
    df = df.replace(0, np.nan)
    
    return df


def _read_exp_file(file_path):
    
    if not os.path.exists(file_path):
        fail(f'Column group file {file_path} does not exist')
    
    import csv
    
    if '\t' in open(file_path).read():
      delimiter = '\t'
    else:
      delimiter = ','
    
    n_cols = None
    
    group_dicts = []
    
    with open(file_path, newline='') as csv:
        for row in csv.reader(csvfile, delimiter=delimiter):
            
            if n_cols:
                if len(row) != n_cols:
                    fail(f'Values appear to be missing in the experimental design file {file_path}')
                
            else:
                n_cols = len(row)              
                group_dicts = [defaultdict(list) for x in n_cols-1]
            
            col_name, *groups = row
            for i, group in enumerate(groups):
                group_dicts[i][group.strip()].append(col_name.strip())
                
            
    if not group_dicts:
        fail(f'Column group file {file_path} did not appear to contain anything useful')
    
    return group_dicts # a list of {group_nameA:[col_name1, col_name2], group_nameB:[col_name3, col_name4],}
    
      
def _split_col_groups(df, group_spec):
    
    cols = list(df.columns)
    n_cols = len(cols)
    
    orig = group_spec
    group_dicts = []
    group_spec = ''.join([x for x in group_spec if x in '0123456789,: '])
    
    while ' ,' in group_spec:
        group_spec = group_spec.replace(' ,',',')
    
    for comparison in group_spec.split():
        group_dict = defaultdict(list)
        
        for g, groupstr in enumerate(comparison.split(':')):
            if g < 26:
              gname = chr(ord('A')+g) # A, B, C
            else:
              gname = f'Grp{g}' # Failsafe
               
            group_dict = {}
            
            for col in groupstr.split(','):
                if col.isdigit():
                    col = int(col)-1
                    
                    if 0 <= col < n_cols:
                       col = cols[col]
                       group_dict[gname].append(col)
                       
                    else:
                      fail(f'Column {col} invalid')
                    
                else:
                    fail(f'Column specification {col} not numeric')
            
        group_dicts.append(group_dict)

    return group_dicts
    

def _ask_user_col_idx_groups(cols, idx_cols, ref_cols, col_groups):
    
    n_cols = len(cols)
    _report('\nAvailable column numbers and names:')

    for i, col in enumerate(cols, 1):
      _report(f'  {i:>3} : {col}')
    
    if not idx_cols
        idx_col = True
 
        while idx_col:
          idx_col = input('\rPlease enter number for index column [or blank to exit] and press <RETURN>:').strip()
 
          if idx_col.isdigit():
              if 0 < int(idx_col) <= n_cols:
                  idx_col = int(idx_col)
                  break
              else:
                  warn(f'Index column "{idx_col}" invalid')
 
          else:
            warn('Index column selection must be a number')
 
        if not idx_col:
            fail(f'No index column specified')
        
        idx_cols = [idx_col]

    if not ref_cols
        ref_col = True
 
        while ref_col:
          ref_col = input('\rPlease enter number for reference column [or blank to exit] and press <RETURN>:').strip()
 
          if ref_col.isdigit():
              if 0 < int(ref_col) <= n_cols:
                  ref_col = int(ref_col)
                  break
              else:
                  warn(f'Reference column "{ref_col}" invalid')
 
          else:
            warn('Reference column selection must be a number')
        
        ref_cols = [ref_col]
    
    if not col_groups:
        col_groups = []
        group = True
        g = 1
 
        while group:
            group = input(f'\r{group} Please enter column numbers for a group {g} [or blank to finish] and press <RETURN>:').strip()
 
            for x in group:
                if x not inidx_col '0123456789 ':
                    warn('Please use only digits and spaces')
                    break
            else:
                group = [int(x) for x in group.split()]
 
                for x in group:
                    if not (0 < x <= n_cols):
                        warn(f'Column number {x} invalid')
                        break
 
                else:
                    col_groups.append(group)
                    g += 1
 
    if not col_groups:
        fail(f'No column groups specified')
    
    return idx_cols, ref_cols, col_groups
    

def _check_cols(df, selection):
     
     cols = list(df.columns)
     n_cols = len(cols)
     
     vaid = []
     for col in selection:
         if (col not in df) and col.isdigit():
             col = int(col)-1  # 1-based to 0-based
              
             if 0 < col < n_cols:
                 vaid.append(cols[col])
                  
         elif col in df:
             vaid.append(col)
     
     return valid
     
    
def lava(in_path, pdf_path=None, software=VALID_SOFTWARE, idx_cols=None,  ref_cols=None, col_groups=None, exp_path=None,
         remove_contaminents=True, f_thresh=DEFALUT_LOGFC, p_thresh=DEFALUT_MAXP, quiet=False, colors=(DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_LOW_COLOR),
         split_x=False, show_labels=False, hq_only=False, print_hits=False):
    
    if pdf_path:
      pdf = PdfPages(pdf_path)
    else:
      pdf = None
    
    color_pos, color_neg, color_low =   colors
    df =  _read_data(in_path)
    cols = list(df.columns)
    n_cols = len(cols)
    
    # Read from group file if possible
    if group_path:
        col_groups = _read_exp_file(exp_path)
   
    elif col_groups:
        col_groups = _split_col_groups(df, col_groups)
    
    if idx_cols:
        idx_cols = _check_cols(df,nidx_cols)
        if not idx_cols:
            warn('No valid index columns found')

    if ref_cols:
        ref_cols = _check_cols(df, ref_cols)
        
        if not ref_cols:
            warn('No valid reference columns found')
    
    #### Is this wanted?
    #if not col_groups or not idx_cols: # Prompt at command line
    #    idx_cols, ref_cols, col_groups = _ask_user_col_idx_groups(cols, idx_cols, ref_cols, col_groups)
        
    # Validity and sanity check, convert ints to text labels
    
    #### REDO this #####
    if idx_cols and col_groups:
        if isinstance(idx_col, int):
            if idx_col < 0:
                fail(f'Index columns may not be less than 1')
            elif idx_col >= n_cols:
                fail(f'Index columns can be at most {len(cols)}')
            
            idx_col = cols[idx_col]
        
        elif idx_col not in cols:
            fail(f'Index column "{idx_col}" not in present dataset.')
        
        found = set()
        for g, col_group in enumerate(col_groups):
            group_set = frozenset(col_group)
            
            if len(group_set) < len(col_group):
                fail('Group {g+1} contains duplicate columns')
            
            if group_set in found:
                warn(f'Ignoring exact duplicate group {g+1}')
                continue
            
            for i, col in enumerate(col_group):
                if isinstance(col, int):
                   if col < 0:
                       fail(f'Index columns may not be less than 1')
                   elif col >= n_cols:
                       fail(f'Index columns can be at most {len(cols)}')
                   
                   col =  cols[col]
                   col_group[i] = col
 
                elif col not in cols:
                    fail(f'Column "{col}" not in present dataset.')
                
                if col == idx_col:
                    fail(f'Index column "{idx_col}" may not be within the column groups')
            
            found.add(group_set) 
    
    elif not idx_cols:
      fail('No index column specified')

    elif not col_groups:
      fail('No groups specified')
    
    info(f'Using index column {cols.index(idx_col)+1} : "{idx_col}"')
    
    info('Column groups:')
    value_cols = []
    
    for g, group in enumerate(col_groups, 1):
        group_str = ', '.join(group)
        info(f'  {g} : {group_str}')
        
        for col in group:
            if col not in value_cols:
                value_cols.append(col)
         
    pre_cull = df.shape
    
    if remove_contaminents and software in SOFT_DN:
        df = df[df['First.Protein.Description'].str.contains('Keratin') == False]
 
    if remove_contaminents and software in SOFT_MQ:
        df = df[df['Majority protein IDs'].str.startswith('CON_') == False]
        df = df[df['Majority protein IDs'].str.startswith('REV_') == False]

    if remove_contaminents and software in SOFT_PD:
        df = df[df['Description'].str.contains('Keratin') == False]

    post_cull = df.shape

    culled_rows = pre_cull[0] - post_cull[0]
    info(f'Removed {culled_rows} contaminant rows from {pre_cull[0]} total')
    
    if not quiet:
      info('Your data now has the following structure:')
      print(df.head(10))
    
    info('Plotting bulk properties')
    util.plot_sums_means_nans(df, value_cols, pdf)
    
    # auto warning in dodgy columns ? - automatic or user prommpted drop?
    #columns_to_drop = []
    #df = df.drop(columns=columns_to_drop)
    #print('previous numeric columns:', value_cols, '\n')
    #value_cols = [col for col in value_cols if col not in columns_to_drop]
    #print('updated numeric columns:', value_cols)

    df = np.log2(df[value_cols])
    
    # box plots
    info('Plotting box-plots')
    util.boxplots(df, value_cols, pdf)

    # histograms
    info('Plotting distributions')
    util.histograms(df, value_cols, pdf)
    
    # pair scatters
    
    ncols = 4
    info('Plotting comparisons')
    util.xy_plots(df, value_cols, ncols, pdf)
    
    ## Extend this to all column pairs, reporting rho?
    
    ## Add correlation matrix density matrix?
    
    ## Add summary page of groups and group pairing?
    
    # Pair groups
    
    
    
    # Pair volcanos
    
    # PDF output
    
    if pdf:
      pdf.close()
      info('Wrote PDF output to {}'.format(pdf_path))    
      
      
def _check_color(color_spec):

    if not is_color_like(color_spec):
        fail(f'Color specicifcation "{is_color_like}" could not be interpreted. See: https://matplotlib.org/stable/users/explain/colors/colors.html')
    

def main(argv=None):

    from argparse import ArgumentParser

    if argv is None:
        argv = sys.argv[1:]
 
    arg_parse = ArgumentParser(prog='proteome_consec_pep_scan.py', description='Find sequences with the longest subsequences that match consecutive arrays of sub-peptide sequences',
                               epilog='For further help contact https://github.com/tempeparsons', prefix_chars='-', add_help=True)
 
    # Do we allow multiple (similar) inputs? ' ; nargs='+'
    arg_parse.add_argument(metavar='INPUT_FILE', dest='d',
                           help='The input file path to read data from; the DIA results text file (protein-level, not peptide-level)')

    arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', default=None,
                           help=f"The name of the software used to process the data. Available choices: {' ,'.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}")

    arg_parse.add_argument('-o', '--out-pdf', dest="o", metavar='FILE_PATH', default=None,
                           help=f"Optional path to save graphs as a PDF file. If not specified graphics will be plotted to screen.")

    arg_parse.add_argument('-k', '--ignore-keratin',dest="k",  default=False, action='store_true',
                           help='Whether to ignore keratin contaminents. If NOT set contaminents will be removed.')
    
    arg_parse.add_argument('-c', '--columns', dest="c", metavar='COLUMN_SELECTION', nargs='+', default=None,
                           help='If not already specified within an experimental design file (-e), this option defines column/samples which should be grouped and compared. ' \
                                'This option is designed for use within automated pipelines; for regular users the -e option should be easier.' \ 
                                'Column/sample groups are specified in the form "1,2,3:4,5,6:7,8,9  3,4:5,6:10,11", with numbers ' \
                                'starting at 1, identifying the data column. Columns within the same group are joined with commas "," and '\
                                'two or more groups to be compared with one another are joined with colons ":". Several comparisons may be ' \ 
                                'specified, separated by spaces. Columns numbers may appear in multiple groups/comparisons.') 

    arg_parse.add_argument('-i', '--index-columns', dest="i", metavar='INDEX_COLUMN', default=None,
                           help='The names, or numbers starting from 1, of one or more input columns used to uniquely index the data rows. ' \
                                'If not specified columns named will be sought from the degsult list: {' '.join(DEFAULT_INDEX_COLS)}.')

    arg_parse.add_argument('-r', '--ref-groups', dest="r", metavar='GROUP_NAMES', nargs='+',  default=None,
                           help='An optional list of names specifiying which group(s) of samples/columns are considered reference groups.'
                                'Typicaly this would include a control group. These names should match the experiment design file (-e). ' \
                                'When reference groups are present other relevant groups will be compared only to these references, and ' \
                                'not among themselves. If there is no reference group all groups will be compared to all relavent others. ')
    
    arg_parse.add_argument('-e', '--experiment-table', dest="e", metavar='FILE_PATH', default=None, 
                           help='The location of a experimental design file. This file is a tab or comma-separated text file containing a table ' \
                                'that relates the name of the data samples (input columns) to which experimental groups they belong. ' \
                                'Each line should first give the name of the input column and then list the groups (or categories) ' \
                                'to which it belongs. Groups will be compared according to their order, i.e. comparisons are only ' \
                                'made within the groups specified each one column of the experiment design table. ' \
                                'Example contents: \n\nSample1\tControl\tFemale\nSample2\tControl\tMale\nSample3\tMutant\tFeale\nSample4\tMutant\tMale')
    
    arg_parse.add_argument('-q', '--quiet-mode', dest="q", action='store_true',
                           help='Proceed quietly, without user interaction; implicitly accept defaults to user queries')
    
    arg_parse.add_argument('-f', '--min-log-fold', dest="f", metavar='MIN_RATIO', type=float, default=DEFALUT_LOGFC,
                           help=f'Minimum log-fold change for potentially significant hits. Default: {DEFALUT_LOGFC}')
   
    arg_parse.add_argument('-p', '--max-p-value', dest="p", metavar='MAX_PVALUE', type=float, default=DEFALUT_MAXP,
                           help=f'Maximum threshold p-value (as a %) for selecting significant hits. Default: {DEFALUT_MAXP}')
    
    arg_parse.add_argument('--pos-color', dest="pos-color", metavar='COLOR', default=DEFAULT_POS_COLOR,
                           help=f'Optional color specification (used by matplotlib) for positive hits on volcano plots, e.g. "blue" or "#0000FF". Default: {DEFAULT_POS_COLOR}')
    
    arg_parse.add_argument('--neg-color', dest="neg-color", metavar='COLOR', default=DEFAULT_NEG_COLOR,
                           help=f'Optional color specification (used by matplotlib) for negative hits on volcano plots, e.g. "red" or "#FF0000". Default: {DEFAULT_NEG_COLOR}')
                           
    arg_parse.add_argument('--low-color', dest="low-color", metavar='COLOR', default=DEFAULT_LOW_COLOR,
                           help=f'Optional color specification (used by matplotlib) insignificant points on volcano plots, e.g. "grey" or "#808080". Default: {DEFAULT_LOW_COLOR}')

    arg_parse.add_argument('--split-x', dest="split-x", action='store_true',
                           help='Use to split the X-axis on volcano plots; useful if you expect a wide range of fold-changes')

    arg_parse.add_argument('--no-labels', dest="no-labels", action='store_true',
                           help='If set, suppress the labelling of significant hits in the volcano plots')
    
    arg_parse.add_argument('--hq-only', dest="hq-only", action='store_true',
                           help='If set, plot only the high-quality points on volcano plots. Otherwise low quality points ' \
                                 '(thous with s single data point in one group) will be plotted, albeit differently. Data points with only ' \
                                 'one value in both compared groups are never plotted.')

    arg_parse.add_argument('--print-hits', dest="print-hits", action='store_true',
                           help='If set, causes the names of the significant hits to be printed below the volcano plots.')
  
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # 
    # Add options to rename columns?
    # Add markers
    # Add log file (off of _report)
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # 
    
    # Enable multi-index
    # Add reference groups; else all -vs all ; refs are compared with all other non-reference groups
    
    args = vars(arg_parse.parse_args(argv))
 
    in_path = args['d']
    pdf_path = args['o']
    software = args['s']
    remove_contaminents = not args['k']
    idx_cols = args['i']
    ref_cols = args['r']
    columns = args['c']
    exp_path = args['e']
    f_thresh = args['f']
    p_thresh = args['p']
    quiet = args['q']
    
    split_x = args['split-x']
    
    pos_color = args['pos-color']
    neg_color = args['neg-color']
    low_color = args['low-color']
    show_labels = not args['no-labels']
    hq_only = args['hq-only']
    print_hits = args['print-hits']
    
    if columns:
        columns = ' '.join(columns)
    
    if exp_path and columns:
        fail('Please specifiy either an experiment design file (-e) or a column selection (-c), not both.')
      
    if not (exp_path or columns):
        if quiet:
            fail('Either an experiment design file (-e), or a column specification (-c) must be given with quiet mode -q')
      
        else:
            warn('No data colums specified; user will be prompted')
    
    if f_thresh <= 1.0:
        fail('Minimum log fold change threshold (-f) must be greeater than 1.0')
    
    if p_thresh > 50.0:
        fail('Maximum p-value threshold must be < 50.0')
    
    if p_thresh < 0.0:
        fail('Maximum p-value threshold must be positive')
    
    _check_color(pos_color)
    _check_color(neg_color)
    _check_color(low_color)
    
    colors = (pos_color, neg_color, low_color)
    
    lava(in_path, pdf_path, software, idx_cols, ref_cols, columns, exp_path, remove_contaminents,
         f_thresh, p_thresh, quiet, colors, split_x, show_labels, hq_only, print_hits)


if __name__ == '__main__':
    
    main()
    
"""
python3 lava.py VolcCLI_PD1.xlsx -i 4 -c 34, 37, 40, 41  35, 36, 38, 39, 42
python3 lava.py VolcCLI_PD1.xlsx -o VolcCLI_PD1.pdf -g groups_example.txt 


"""
    

