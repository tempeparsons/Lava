import os, sys
import numpy as np
import pandas as pd

from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.backends.backend_pdf import PdfPages

import lava_util as util
import lava_plots as plots

VERSION = 0.1

# Allow shirt codes for software
SOFT_DN = ('DIANN', 'DN')
SOFT_MQ = ('MaxQuant', 'MQ')
SOFT_PD = ('ProteomeDiscoverer', 'PD')
 
VALID_SOFTWARE = SOFT_DN + SOFT_MQ + SOFT_PD
DEFAULT_SOFTWARE = SOFT_PD[0]

FILE_EXTS_EXL = ('.xls','.xlsx')
FILE_EXTS_TXT = ('.txt',)
FILE_EXTS_CSV = ('.csv','.tsv')

DEFALUT_FC = 2.0
DEFALUT_MAXP  = 4.32
DEFAULT_POS_COLOR = '#0080FF'
DEFAULT_NEG_COLOR = '#FF0000'
DEFAULT_LOW_COLOR = '#FFD0D0'

DEFAULT_INDEX_COLS = ('Accession','Protein.ID','Protein ID')

LOG_FILE = None

def warn(txt):
   
    _report('WARNING: ' + txt)


def info(txt):
   
    _report('INFO: ' + txt)


def fail(txt):
   
    _report('FAILURE: ' + txt)
    _report('EXIT')
    sys.exit(0)


def _report(txt):
    
    if LOG_FILE:
        LOG_FILE.write(txt + '\n')
    
    print(txt)


def _read_data(file_path):
    
    info(f'Reading {file_path}')
    
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


def _read_exp_file(df, file_path):
    
    if not os.path.exists(file_path):
        fail(f'Column group file {file_path} does not exist')
    
    import csv
    
    if '\t' in open(file_path).read():
        delimiter = '\t'
    else:
        delimiter = ','
    
    n_cols = None
    
    group_dicts = []
    
    with open(file_path, newline='') as file_obj:
        for row in csv.reader(file_obj, delimiter=delimiter):
            if not row:
                continue
            
            if row[0].startswith('#'):
                continue
                
            if n_cols:
                if len(row) != n_cols:
                    fail(f'Values appear to be missing in the experimental design file {file_path}')
                
            else:
                n_cols = len(row)              
                group_dicts = [defaultdict(list) for x in range(n_cols-1)]
            
            col_name, *groups = row
            col_name = col_name.strip()
            
            if col_name not  in df:
              fail(f'Sample column named {col_name} is not present in input data')
            
            for i, group in enumerate(groups):
                group_dicts[i][group.strip()].append(col_name)
                
            
    if not group_dicts:
        fail(f'Experimental design file {file_path} did not appear to contain anything useful')
    
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
                      fail(f'Sample column number {col+1} invalid')
                    
                else:
                    fail(f'Sample column specification {col} not numeric')
            
        group_dicts.append(group_dict)

    return group_dicts
    

def _ask_user_col_idx_groups(cols, idx_cols, ref_groups, col_groups):
    
    n_cols = len(cols)
    _report('\nAvailable column numbers and names:')

    for i, col in enumerate(cols, 1):
        _report(f'  {i:>3} : {col}')
    
    if not idx_cols:
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

    if not ref_groups:
        ref_group = True
 
        while ref_group:
          ref_group = input('\rPlease enter number for reference column [or blank to exit] and press <RETURN>:').strip()
 
          if ref_group.isdigit():
              if 0 < int(ref_group) <= n_cols:
                  ref_group = int(ref_group)
                  break
              else:
                  warn(f'Reference column "{ref_group}" invalid')
 
          else:
              warn('Reference column selection must be a number')
        
        ref_groups = [ref_group]
    
    if not col_groups:
        col_groups = []
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
    
    return idx_cols, ref_groups, col_groups
    

def _check_cols(df, selection):
     
     cols = list(df.columns)
     n_cols = len(cols)
     
     valid = []
     for col in selection:
         if (col not in df) and col.isdigit():
             col = int(col)-1  # 1-based to 0-based
              
             if 0 < col < n_cols:
                 valid.append(cols[col])
             else:   
                 warn('Column number {col} not valid') 
                    
         elif col in df:
             valid.append(col)
     
     return valid
     
    
def lava(in_path, pdf_path=None, software=VALID_SOFTWARE, idx_cols=None,  ref_groups=None, col_groups=None, exp_path=None,
         remove_contaminents=True, f_thresh=DEFALUT_FC, p_thresh=DEFALUT_MAXP, quiet=False, colors=(DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_LOW_COLOR),
         split_x=False, hit_labels=False, hq_only=False, print_hits=False):
    
    info(f'Lava version {VERSION}')
    
    if pdf_path:
        pdf = PdfPages(pdf_path)
    else:
        pdf = None
    
    color_pos, color_neg, color_low =   colors
    df =  _read_data(in_path)
    cols = list(df.columns)
    n_cols = len(cols)
    
    # Read from group file if possible
    if exp_path:
        group_dicts = _read_exp_file(df, exp_path)
   
    elif col_groups:
        group_dicts = _split_col_groups(df, col_groups)
    
    if idx_cols:
        idx_cols = _check_cols(df, idx_cols)
        if not idx_cols:
            fail('No valid index columns found')
    
    else:
        idx_cols = [x for x in DEFAULT_INDEX_COLS if x in df]
      
        if not idx_cols:
            fail('No valid index columns found after falling back to default; {' '.join(DEFAULT_INDEX_COLS)}')
    
    #### Change to Tkinter GUI
    #if not col_groups or not idx_cols: # Prompt at command line
    #    idx_cols, ref_groups, col_groups = _ask_user_col_idx_groups(cols, idx_cols, ref_groups, col_groups)
        
    if not idx_cols:
        fail('No index column specified')

    elif not group_dicts:
        fail('No groups specified')
    
    for idx_col in idx_cols:
        info(f'Using index column {cols.index(idx_col)+1} : "{idx_col}"')
    
    df.set_index(idx_cols, inplace=True)
    
    info('Sample column groups for comparison:')
    value_cols = []
    
    groups = set()
    group_cols = {}
    for i, group_dict in enumerate(group_dicts):
       info(f'Comparison {i+1}:')
       
       for gname in sorted(group_dict):
           cols = group_dict[gname]
           
           if gname in group_cols:
               fail(f'Sample group "{gname}" is repeated in a different context; group names should only be present in one experimental design column.')
           
           if len(cols) < 2:
               warn(f'Sample group "{gname}" contains only one sample; cannot be used')
               del group_dict[gname]
               continue
               
       
           groups.add(gname)
           info(f'    Group "{gname}": {" ".join(cols)}')
           
           for col in cols:
              if col not in value_cols:
                  if col in idx_cols:
                      fail(f'Ovelap between index columns and sample columns not permitted; "{col}" is an index column.')
                
                  value_cols.append(col)
    
       group_cols.update(group_dict)
       
    if ref_groups:
        ref_groups = set(ref_groups)
        info(f'Reference sample groups:  {" ".join(ref_groups)}')
        
        for ref_group in ref_groups:
           if ref_group not in groups:
               fail('Reference sample group "{ref_group}" not found in experimental design. Available group names: {' '.join(sorted(groups))}')
    
    else:
        info(f'No reference groups specified')
    
    info(f'The comparisons will be made between the folloing pairs of sample groups:')
    pairs = []
    
    for group_dict in group_dicts:
        groups = set(group_dict.keys())
        refs = groups & ref_groups
        
        if refs:
            non_ref = sorted(groups - ref_groups)
            for g1 in sorted(refs):
                for g2 in non_ref:
                    info(f'   {g1}  -  {g2}')
                    pairs.append((g1, g2))
          
        else:
            groups = sorted(groups)
            for i, g1 in enumerate(groups[:-1]):
                for g2 in groups[i+1:]:
                    info(f'   {g1}  -  {g2}')
                    pairs.append((g1, g2))
         
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
    plots.plot_sums_means_nans(df, value_cols, pdf)
    
    ### auto warning in dodgy columns ? - automatic or user prommpted drop?
    #columns_to_drop = []
    #df = df.drop(columns=columns_to_drop)
    #print('previous numeric columns:', value_cols, '\n')
    #value_cols = [col for col in value_cols if col not in columns_to_drop]
    #print('updated numeric columns:', value_cols)

    df = np.log2(df[value_cols])
    
    # box plots
    info('Plotting box-plots')
    plots.boxplots(df, value_cols, pdf)

    # histograms
    info('Plotting distributions')
    plots.histograms(df, value_cols, pdf)

    info('Plotting correlations')
    plots.correlation_plot(df, value_cols, pdf)
    
    # pair scatters
    ncols = 4
    info('Plotting comparisons')
    plots.xy_plots(df, value_cols, ncols, pdf)
    
    
    FCdict = {}
    Zdict = {}

    dftest = []
    for g1, g2 in pairs:
        cols1 = group_cols[g1]
        cols2 = group_cols[g2]
        df2 = df[cols1+cols2].copy()
        
        key = f'{g1}:{g2}'
        FCdict[key] = [df2, cols1, cols2]
        Zdict[key]  = [df2, cols1, cols2]
 
 
    overFCthresh = {}
    dfnames = []
    FCdfs = []
    for k, (df2, cols1, cols2) in FCdict.items():
        grp1 = [c for c in df2.columns if c in cols1]
        grp2 = [c for c in df2.columns if c in cols2]
 
        df2['nobs_grp1'] = df2.loc[:,grp1].count(axis=1)
        df2['nobs_grp2'] = df2.loc[:,grp2].count(axis=1)
        df2 = util.remove_rows(df2, len(grp1), len(grp2))
        df2 = df2.drop(columns = ['nobs_grp1', 'nobs_grp2'])

        df2['mean_grp1'] = df2.loc[:,grp1].mean(axis=1)
        df2['mean_grp2'] = df2.loc[:,grp2].mean(axis=1)
        df2['grp1grp2_FC'] = df2['mean_grp1'] - df2['mean_grp2']
 
        FCdfs.append(df2)
        dfnames.append(k)
 
        p, n = util.prpn_log2FCs_over_threshold(df2, f_thresh)
        overFCthresh[k] = (p, n)
 
 

    Zdfs = []
    q95s = []
    for k, (df, cols1, cols2) in Zdict.items():
        grp1 = [c for c in df.columns if c in cols1]
        grp2 = [c for c in df.columns if c in cols2]
 
        df = util.make_znorm(df)
        df['nobs_grp1'] = df.loc[:,grp1].count(axis=1)
        df['nobs_grp2'] = df.loc[:,grp2].count(axis=1)
        df = util.remove_rows(df, len(grp1), len(grp2))
        df = df.drop(columns = ['nobs_grp1', 'nobs_grp2'])
 
        df['zmean_grp1'] = df.loc[:,grp1].mean(axis=1)
        df['zstd_grp1'] = df.loc[:,grp1].std(axis=1)
        
        q95grp1 = df.zstd_grp1.quantile(q=0.95)
        df['zstd_grp1_q95'] = df['zstd_grp1'].fillna(q95grp1)
        df['znobs_grp1'] = df.loc[:,grp1].count(axis=1)
        df['znobs_grp1_q95'] = np.where(df['znobs_grp1'] == 1, 2, df['znobs_grp1'])
        df['zmean_grp2'] = df.loc[:,grp2].mean(axis=1)
        df['zstd_grp2'] = df.loc[:,grp2].std(axis=1)
        
        q95grp2 = df.zstd_grp2.quantile(q=0.95)
        df['zstd_grp2_q95'] = df['zstd_grp2'].fillna(q95grp2)
        df['znobs_grp2'] = df.loc[:,grp2].count(axis=1)
        df['znobs_grp2_q95'] = np.where(df['znobs_grp2'] == 1, 2, df['znobs_grp2'])
 
        df = util.ttest_from_stats_eqvar(df)
 
        q95 = (q95grp1, q95grp2)
        q95s.append(q95)
 
        Zdfs.append(df)

    FZdfs = []
    for i, df in enumerate(FCdfs):
        FZdf = pd.concat([df, Zdfs[i]], axis=1)
        FZdfs.append(FZdf[['grp1grp2_FC', 'pvals', 'Tpvals', 'Tpvals_q95', 'zstd_grp1_q95', 'zstd_grp2_q95']])
 
    plotdict = dict(zip(dfnames, FZdfs))
    q95dict = dict(zip(dfnames, q95s))

    ##review that you're about the plot the expected pairs
    ##review what percentage of proteins will be higher-lower than the positive-negative fold change limits

    #for k, v in plotdict.items():
    #    info(f'you are plotting these pairs: {k}')
    
    info(f'Using log2 fold-change threshold: {f_thresh}')
    for k, (pos, neg) in overFCthresh.items():
        info(f'  - pair {k} proteins over theshold: pos; {pos}% neg; {neg}%')

    # if you got too high or low a %age of proteins over the FC threshold, type your new number instead of 'same'
    # to keep it the same, just leave it as is.

    # now examine the distrbutions of p-values and fold changes for all your pairs.
    # the majority of fold changes should be normally distrbuted around zero
    # the distrbution of p-values will depend on your experiment

    plots.pvalFC_hists(plotdict, pdf)

    #try:
    #    splitxlist = util.get_splitxlist(key)
    #
    #except IndexError:
    #    fail('x-axis range is not large enough for split x-axis. Remove the -x option')
    
    for a, b in pairs:
        pair_name = f'{a}:{b}'
        plots.volcano(pdf, pair_name, plotdict[pair_name], q95dict[pair_name], f_thresh, p_thresh, colors,
                      hq_only, hit_labels, print_hits, lw=0.25, ls='--')
 
    
    
    """

    if split_x:
        fig, ax = plt.subplots(1, 3, figsize = (10,8), sharey=True)
        fig.subplots_adjust(wspace=0.05)

        a0 = ax.flat[0]
        a1 = ax.flat[1]
        a2 = ax.flat[2]

        a0.set_xlim(splitxlist[0], splitxlist[1])
        a1.set_xlim(splitxlist[2], splitxlist[3])
        a2.set_xlim(splitxlist[4], splitxlist[5])

        a1.axvline(FClim, ls=ls, c='k', alpha=0.5, lw=lw)
        a1.axvline(-FClim, ls=ls, c='k', alpha=0.5, lw=lw)
        for a in ax.flat:
            a.axhline(pthresh,ls=ls, c='k', alpha=0.5, lw=lw)

        a1.set_title(df_pair, fontsize=18)

        a0.spines['right'].set_visible(False)
        a1.spines['left'].set_visible(False)
        a1.spines['right'].set_visible(False)
        a2.spines['left'].set_visible(False)

        d = .8
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                      linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        a0.plot([1,1], [1,0], transform=a0.transAxes, **kwargs)
        a1.plot([0,0], [1,0], [1,1], [1,0], transform=a1.transAxes, **kwargs)
        a2.plot([0,0], [0,1], transform=a2.transAxes, **kwargs)
        a0.yaxis.tick_left()
        a0.set_ylabel('-log2 transformed p-value')
        a1.set_xlabel('log2 fold change')
        a2.yaxis.tick_right()
        a1.tick_params(axis='y', which='both',length=0)
        a1.set_yticks([])
        a2.yaxis.set_label_position("right")
        a2.set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
 
        for a in ax.flat:
            a.scatter(plotlist[1], plotlist[2], s=8, color='w', edgecolors='darkgrey', lw=lw)

        marker_text = []

        for i, name in enumerate(plotlist[0]):

            if name in markers:
                marker_text.append((plotlist[1][i], plotlist[2][i], name))
                for a in ax.flat:
                    a.scatter(plotlist[1][i], plotlist[2][i], s=3, color='k')

            if plotlist[1][i] >= FClim and plotlist[2][i] >= pthresh:
                pos_overFC.append(name)
                for a in ax.flat:
                    a.scatter(plotlist[1][i], plotlist[2][i], s=10, color=color_pos)
                if label_significant_proteins:
                    a1.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)
                    a2.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)

            if plotlist[1][i] <= -FClim  and plotlist[2][i] >= pthresh:
                neg_overFC.append(name)
                for a in ax.flat:
                    a.scatter(plotlist[1][i], plotlist[2][i], s=10, color=color_neg)
                if label_significant_proteins:
                    a0.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)
                    a1.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)

            if not plot_high_quality_only:
                if plotlist[3][i] == q951 or plotlist[4][i] == q952:
                    for a in ax.flat:
                        a.scatter(plotlist[1][i], plotlist[2][i], s=6, color=color_low)


        for tx in marker_text:
            for a in ax.flat:
                a.annotate(tx[2], xy=(tx[0], tx[1]),  xytext=(tx[0], tx[1]),
                        arrowprops=dict(arrowstyle='-', fc="k", ec="k", lw=0.5, relpos=(0.25, 0.5)),
                        bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
                        ha="left", va="center", size=6)

        legend_elements = [Line2D([0], [0], marker='o', color='w', label='greater in grp1', markerfacecolor=color_pos, markersize=8, linestyle=''),
                           Line2D([0], [0], marker='o', color='w', label='greater in grp2', markerfacecolor=color_neg, markersize=8, linestyle='')]
        a1.legend(handles=legend_elements)

    """
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
 
    arg_parse = ArgumentParser(prog='lava.py', description='Find sequences with the longest subsequences that match consecutive arrays of sub-peptide sequences',
                               epilog='For further help contact https://github.com/tempeparsons', prefix_chars='-', add_help=True)
 
    # Do we allow multiple (similar) inputs? ' ; nargs='+'
    arg_parse.add_argument(metavar='INPUT_FILE', dest='d',
                           help='The input file path to read data from; the DIA results text file (protein-level, not peptide-level)')
    
    arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', default=None,
                           help=f"The name of the software used to process the data. Available choices: {' ,'.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}")

    arg_parse.add_argument('-o', '--out-pdf', dest="o", metavar='FILE_PATH', default=None,
                           help=f"Optional path to save graphs as a PDF file. If not specified graphics will be plotted to screen.")

    arg_parse.add_argument('-k', '--ignore-contaminants',dest="k",  default=False, action='store_true',
                           help='Whether to ignore kerating and other contaminants. If NOT set contaminents will be removed.')
    
    arg_parse.add_argument('-c', '--columns', dest="c", metavar='COLUMN_SELECTION', nargs='+', default=None,
                           help='If not already specified within an experimental design file (-e), this option defines column/samples which should be grouped and compared. ' \
                                'This option is designed for use within automated pipelines; for regular users the -e option should be easier. ' \
                                'Column/sample groups are specified in the form "1,2,3:4,5,6:7,8,9  3,4:5,6:10,11", with numbers ' \
                                'starting at 1, identifying the data column. Columns within the same group are joined with commas "," and '\
                                'two or more groups to be compared with one another are joined with colons ":". Several comparisons may be ' \
                                'specified, separated by spaces. Columns numbers may appear in multiple groups/comparisons.') 

    arg_parse.add_argument('-i', '--index-columns', dest="i", metavar='INDEX_COLUMN', default=None,
                           help=f'The names, or numbers starting from 1, of one or more input columns used to uniquely index the data rows. ' \
                                'If not specified columns named will be sought from the default list: {' '.join(DEFAULT_INDEX_COLS)}.')

    arg_parse.add_argument('-r', '--ref-groups', dest="r", metavar='GROUP_NAME', nargs='+',  default=None,
                           help='An optional list of names specifiying which group(s) of samples/columns are considered reference groups.'
                                'Typicaly this would include a control group. These names should match the group/category names in the experiment design file (-e). ' \
                                'When reference groups are present other relevant groups will be compared only to these references, and ' \
                                'not among themselves. If there is no reference group all groups will be compared to all relavent others. ')
    
    arg_parse.add_argument('-e', '--experiment-table', dest="e", metavar='FILE_PATH', default=None, 
                           help='The location of a experimental design file. This file is a tab or comma-separated text file containing a table ' \
                                'that relates the name of the data samples (input columns) to which experimental groups they belong. ' \
                                'Each line should first give the name of the input column and then list the groups (or categories) ' \
                                'to which it belongs. Groups will be compared according to their order, i.e. comparisons are only ' \
                                'made within the groups specified each one column of the experiment design table. ' \
                                'Example contents: \n\nSample1\tControl\tFemale\nSample2\tControl\tMale\nSample3\tMutant\tFemale\nSample4\tMutant\tMale')
 
    arg_parse.add_argument('-l', '--log-status', dest="l", action='store_true',
                           help=f"When set, writes status information to a log file; th elog file path will be based on the input path.")
    
    arg_parse.add_argument('-q', '--quiet-mode', dest="q", action='store_true',
                           help='Proceed quietly, without user interaction; implicitly accept defaults to user queries')
    
    arg_parse.add_argument('-f', '--min-log-fold', dest="f", metavar='MIN_RATIO', type=float, default=DEFALUT_FC,
                           help=f'Minimum log2-fold change for potentially significant hits. Default: {DEFALUT_FC}')
   
    arg_parse.add_argument('-p', '--max-p-value', dest="p", metavar='MAX_PVALUE', type=float, default=DEFALUT_MAXP,
                           help=f'Maximum threshold p-value (as a percentage) for selecting significant hits. Default: {DEFALUT_MAXP}')
    
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
    # Export proteins in each quadrant ; 4 lists for each comparison
    # Export table of log2FC, FC and p-values
    # Add parameters/setup page, summary of groups and group pairing
    # Add CLI text colors
    
    # (M) Add file path detection in col names ; auto-truncation
    # (M) Add markers : list of index col ids to highlight
    # (M) Add grouping info to xy correlation plots
    
    # (E) Tweak volcanos; unlogged axis labels 
    # (H) Split volcanos
            
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # 
    
    args = vars(arg_parse.parse_args(argv))
 
    in_path = args['d']
    pdf_path = args['o']
    software = args['s']
    remove_contaminents = not args['k']
    idx_cols = args['i']
    ref_groups = args['r']
    columns = args['c']
    exp_path = args['e']
    f_thresh = args['f']
    p_thresh = args['p']
    quiet = args['q']
    log = args['l']
    
    split_x = args['split-x']
    
    pos_color = args['pos-color']
    neg_color = args['neg-color']
    low_color = args['low-color']
    hit_labels = not args['no-labels']
    hq_only = args['hq-only']
    print_hits = args['print-hits']
    
    if log:
        file_root = os.path.splitext(in_path)[0]
        i = 0
        log_path = f'{file_root}_log{i}.txt'
        while os.path.exists(log_path):
          i += 1
          log_path = f'{file_root}_log{i}.txt'
        
        global LOG_FILE
        LOG_FILE = open(log_path, 'w')
    
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
        fail('Minimum fold change threshold (-f) must be greeater than 1.0')
    
    f_thresh = np.log2(f_thresh)   
        
    if p_thresh > 50.0:
        fail('Maximum p-value threshold must be < 50.0')
    
    if p_thresh < 0.0:
        fail('Maximum p-value threshold must be positive')
    
    _check_color(pos_color)
    _check_color(neg_color)
    _check_color(low_color)
    
    colors = (pos_color, neg_color, low_color)
    
    lava(in_path, pdf_path, software, idx_cols, ref_groups, columns, exp_path, remove_contaminents,
         f_thresh, p_thresh, quiet, colors, split_x, hit_labels, hq_only, print_hits)


if __name__ == '__main__':
    
    main()
    
"""
python3 lava.py VolcCLI_PD1.xlsx -o VolcCLI_PD1.pdf -e groups_example.txt -l 
python3 lava.py VolcCLI_PD1.xlsx -o VolcCLI_PD1.pdf -e groups_example.txt -l -r A


"""
    

