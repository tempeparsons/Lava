import os, sys
import numpy as np
import pandas as pd

from matplotlib.colors import is_color_like

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


def _read_group_file(file_path):
    
    if not os.path.exists(file_path):
        fail(f'Column group file {file_path} does not exist')
    
    groups = []
    
    with open(file_path) as file_obj:
        idx_col = file_obj.readline().strip()
 
        for line in file_obj:
            cols = line.strip().split('\t')
            cols = [x.strip() for x in cols] # remove any unintended edge whitespace
            print(line, cols)
            groups.append(cols)
            
    if not groups:
        fail(f'Column group file {file_path} did not appear to contain anything useful')
    
    return idx_col, groups
    
      
      
def _split_col_groups(group_spec):
    
    orig = group_spec
    col_groups = []
    
    group_spec = ''.join([x for x in group_spec if x in '0123456789, '])
    
    while ' ,' in group_spec:
      group_spec = group_spec.replace(' ,',',')
    
    for groupstr in group_spec.split():
        group = []
        
        for x in groupstr.split(','):
            if x.isdigit():
                group.append(int(x)-1)
        
        if group:
          col_groups.append(group)

    return col_groups
    

def _ask_user_col_idx_groups(cols):
    
    col_groups = []
    n_cols = len(cols)
    _report('\nAvailable column numbers and names:')

    for i, col in enumerate(cols, 1):
      _report(f'  {i:>3} : {col}')
 
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
        warn('Index column must be a number')
 
    if not idx_col:
        fail(f'No index column specified')
 
    group = True
    g = 1
 
    while group:
        group = input(f'\r{group} Please enter column numbers for a group {g} [or blank to finish] and press <RETURN>:').strip()
 
        for x in group:
            if x not in '0123456789 ':
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
    
    return idx_col, col_groups
    
    
    
def lava(in_path, software=VALID_SOFTWARE, remove_contaminents=True, idx_col=None, col_groups=None, group_path=None,
         f_thresh=DEFALUT_LOGFC, p_thresh=DEFALUT_MAXP, quiet=False, colors=(DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_LOW_COLOR),
         split_x=False, show_labels=False, hq_only=False, print_hits=False):
    
    color_pos, color_neg, color_low =   colors
    df =  _read_data(in_path)
    cols = list(df.columns)
    n_cols = len(cols)
    
    # Read from group file if possible
    if group_path:
        idx_col, col_groups = _read_group_file(group_path)
         
    elif col_groups and idx_col:
        if idx_col.isdigit():
          idx_col = int(idx_col)-1 # 1-based to 0-based
        
        col_groups = _split_col_groups(col_groups)
    
    else: # Prompt at command line
        idx_col, col_groups = _ask_user_col_idx_groups(cols)
        
        
    # Validity and sanity check, convert ints to text labels
    if idx_col and col_groups:
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
      print('your data now has the following structure:')
      df.head(10)
    
    # sum, mean, nan summary
    
    util.plot_sums_means_nans(df, value_cols)
    
    # auto warning in dodgy columns
    
    # box plots
    
    # histograms
    
    # pair scatters
    
    
    # pair volcanos
    
    # PDF output or screen
    
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

    arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', nargs='+', default=None,
                           help=f"The name of the software used to process the data. Available choices: {' ,'.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}")

    arg_parse.add_argument('-k', '--ignore-keratin',dest="k",  default=False, action='store_true',
                           help='Whether to ignore keratin contaminents. If NOT set contaminents will be removed.')
    
    arg_parse.add_argument('-c', '--columns', dest="c", metavar='COLUMN_SELECTION', nargs='+', default=None,
                           help='The input data column selection and grouping. If not specified here, column slections may be specified via a group file (-g), else the user ' \
                                'will be prompted on the command line (unless -y is specified). Column order and group selections are specified in the form "9,10,11 3,4,5", with numbers ' \
                                'starting at 1, grouping related columns with commas, and separating different groups with spaces. This specification may only contain numeric ' \
                                'digits, brackets, commas and spaces.')
    
    arg_parse.add_argument('-i', '--index-column', dest="i", metavar='INDEX_COLUMN', default=None,
                           help='The name, or number (starting from 1) of the input column used to uniquely index the data rows. Must be specified if the -c option is used.')
    
    arg_parse.add_argument('-g', '--group-file', dest="g", metavar='FILE_PATH', default=None,
                           help='The location of a group file; stating the index column, data column order and groupings. This is a TAB-separated file giving' \
                                'the index column (to uniquely refer to each row) on the first line and then, on subsequent lines, grouped column names (or numbers, from 1) ' \
                                'so that columns that go together atr all on the same line. Example contents: \n\nAcession\nWT-1\tWT-2\tWT-3\nMUT-1\tMUT-2\tMUT-3\n\n')
    
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
    # Add markers
    # Add secondary index
    # Add whether dual or alternative index
    # Add options to rename columns
    # Move away from strict TSV for groups file
    # Add log file (off of _report)
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # 
    
    args = vars(arg_parse.parse_args(argv))
 
    in_path = args['d']
    software = args['s']
    remove_contaminents = not args['k']
    idx_col = args['i']
    columns = args['c']
    group_path = args['g']
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
    
    if group_path and columns:
      fail('Please specifiy either a group file (-g) or a column selection (-c), not both.')
    
    if columns and not idx_col:
      fail('No index column (-i) was specified to go with column slection (-c).')
    
    if not (group_path or columns):
      if quiet:
        fail('Either a group file (-g), or a column specification (-c) and index column (-i), must be given with quiet mode -q')
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
    
    lava(in_path, software, remove_contaminents, idx_col, columns, group_path, f_thresh, p_thresh, quiet, colors, split_x, show_labels, hq_only, print_hits)


if __name__ == '__main__':
    
    main()
    
"""
python3 lava.py VolcCLI_PD1.xlsx -i 4 -c 34, 37, 40, 41  35, 36, 38, 39, 42
python3 lava.py VolcCLI_PD1.xlsx -g groups_example.txt


"""
    

