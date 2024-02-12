import os, sys, textwrap, re
import numpy as np
import pandas as pd

from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patheffects as PathEffects
from matplotlib.colors import LinearSegmentedColormap

import lava_util as util
import lava_plots as plots

VERSION = '0.2.0'

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

ANSI_ESC= {'END':'\033[0m', 'BOLD':'\033[1m',
           'ITALIC':'\033[2m', 'UNDERLINE':'\033[3m',
           'BLACK':'\033[30m', 'RED':'\033[31m',
           'GREEN':'\033[32m', 'YELLOW':'\033[33m',
           'BLUE':'\033[34m', 'MAGENTA':'\033[35m',
           'CYAN':'\033[36m', 'WHITE':'\033[37m',
           'GREY':'\033[90m', 'LT_RED':'\033[91m',
           'LT_GREEN':'\033[92m',  'LT_YELLOW':'\033[93m',
           'LT_BLUE':'\033[94m', 'LT_MAGENTA':'\033[95m',
           'LT_CYAN':'\033[96m', 'LT_WHITE':'\033[97m'}


def _color(txt, cname='RED'):
  
    return f"{ANSI_ESC.get(cname, 'RED')}{txt}{ANSI_ESC['END']}"
    

def warn(txt):
   
    _report(_color('WARNING: ','YELLOW') + txt)


def info(txt):
   
    _report(_color('INFO: ', 'BLUE') + txt)


def fail(txt):
   
    _report(_color('FAILURE: ', 'RED') + txt)
    _report(_color('EXIT', 'RED'))
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
    
    # Check for paths as column names
    
    cols = list(df.columns)
    n_map = {}
    rep_count = defaultdict(int)
    
    for col in cols:
        
        if len(col) > 32:
            orig = col
            match_obj_nix = re.search('\/.*?\.[\w:]+', col)
            match_obj_win = re.search(r'\\.*?\.[\w:]+', col)
            
            if match_obj_nix or match_obj_win:
               if match_obj_win:
                 dirsep = '\\'
               else:
                 dirsep = '/'  
               
               parts = col.split(dirsep)[::-1] # last first
               i = 0
               col = parts[0]
               
               while (len(col) < 8) and (i+1 < len(parts)):
                 i += 1
                 col = parts[i] + dirsep + col
               
               prefix = col
               j = 0
               while cols.count(col) > 1: # duplicate
                   j += 1
                   col = f'{prefix}_{j}'              
               
               n_map[orig] = col
    
    if n_map:
       df.rename(columns=nmap, inplace=True)
    
    return df

    
def save_volcano_table(pairs, plotdict,  save_path, f_thresh, p_thresh):
   
   # Rename more col headings, mark q95
   # Colour excel table
   
   nmap = {'grp1grp2_FC':'-log2_fold_change', 'Tpvals':'-log2_pvalue',
           'pvals':'p-value', 'Tpvals_q95':'-log2_Q95_pvalue',
           'zstd_grp1_q95':'Q95_sigma_A','zstd_grp1_q95':'Q95_sigma_B',}
   
   quad_map = {(True,True):'HIT_pos', (True,False):'HIT_neg', (False,True):'fail_pos', (False,False):'fail_neg', }
   
   path_root, file_ext = os.path.splitext(save_path)
   file_ext = file_ext.lower()
   
   keys = [f'{a}:::{b}' for a,b in pairs]
   for key in keys:
       df = plotdict[key]
       lfc = np.array(df['grp1grp2_FC'])
       pvs = np.array(df['Tpvals'])
       n = len(lfc)
       
       quad_cats = [quad_map[((pvs[i] >= p_thresh and (abs(lfc[i]) >= f_thresh)), lfc[i] >= 0)] for i in range(n)]
       
       fc = 2.0 ** (-lfc)
       df.insert(0, 'fold_change', fc)
       df.insert(0, 'hit_class', quad_cats)
       df.sort_values(by=['pvals',], ascending=True, inplace=True)
       df.rename(columns=nmap, inplace=True)
   
   if file_ext in FILE_EXTS_EXL:
       with pd.ExcelWriter(save_path, engine='openpyxl') as writer:
           for key in keys:
               plotdict[key].to_excel(writer, sheet_name=key.replace(':::', '_vs_'))
 
       
       info(f'Saved tables to {save_path}')
   
   else:
       if file_ext in ('.csv'):
           sep = ','
       else:
           sep = '\t'
           
       for key in keys:
           name = key.replace(':::', '_vs_')
           file_path = f'{path_root}_{name}{file_ext}'
           plotdict[key].to_csv(file_path, sep=',', na_rep='nan', quotechar='"')
           info(f'Saved table to {file_path}')


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
     
 
def _check_path(file_path, should_exist=True):
    
    if file_path is not None:
       file_path = os.path.abspath(file_path)
 
       if should_exist:
          if not os.path.exists(file_path):
              fail(f'File "{file_path}" does not exist')
 
          elif not os.path.isfile(file_path):
              fail(f'Location "{file_path}" is not a regular file')
 
          elif os.stat(file_path).st_size == 0:
              fail(f'File "{file_path}" is of zero size')
 
          elif not os.access(file_path, os.R_OK):
              fail(f'File "{file_path}" is not readable')
    
    return file_path

    
def lava(in_path, exp_path=None, software=DEFAULT_SOFTWARE, pdf_path=None, table_path=None,
         idx_cols=None, ref_groups=None, markers=None, col_groups=None, remove_contaminents=True,
         f_thresh=DEFALUT_FC, p_thresh=DEFALUT_MAXP, quiet=False,
         colors=(DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_LOW_COLOR),
         split_x=False, hit_labels=False, hq_only=False):
    
    in_path = _check_path(in_path)
    exp_path = _check_path(exp_path)
    pdf_path = _check_path(pdf_path)
    table_path = _check_path(table_path)
    
    if pdf_path:
        pdf = PdfPages(pdf_path)
    else:
        pdf = None
    
    info(f'Lava version {VERSION}')    
    
    option_report = [(f'Lava Report version {VERSION}', None),
                     ('Input options', ' '.join(sys.argv[1:])),
                     ('Input data file',in_path),
                     ('Exp. design file',exp_path),
                     ('PDF output file',pdf_path),
                     ('Table output file',table_path),
                     ('Input software', software),
                     ('Min. fold-change',f_thresh),
                     ('Max. p-value', f'{p_thresh:.2f}%'),
                     ('Remove contaminents',remove_contaminents),
                     ('Split volcano X-axis',split_x),
                     ('Show hit labels',hit_labels),
                     ('Show high quality',hq_only)]
        
    color_pos, color_neg, color_low =   colors
    df =  _read_data(in_path)
    cols = list(df.columns)
    n_cols = len(cols)
    
    if exp_path:  # Read from exp design file if possible
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
            
    if not idx_cols:
        fail('No index column specified')

    elif not group_dicts:
        fail('No groups specified')
    
    for idx_col in idx_cols:
        info(f'Using index column {cols.index(idx_col)+1} : "{idx_col}"')
    
    option_report.append(('Index column(s)', ' '.join(idx_cols)))
    
    # Markers
    
    if markers:
        avail_idx = set()
        for idx_col in idx_cols:
            avail_idx.update(df[idx_col])
        
        valid = []    
        for marker in markers:
           if marker not in avail_idx:
               warn(f'Marker "{marker}" is not present in index')
           else:
               valid.append(marker)
        
        markers = valid
    
    df.set_index(idx_cols, inplace=True)
        
        
    msg = 'Sample group membership'
    info(msg)
    option_report.append((msg, None))
    value_cols = []
    
    groups = set()
    group_cols = {}
    for i, group_dict in enumerate(group_dicts):
       info(f'Comparison {i+1}:')
       option_report.append(('Comparison', f'{i+1}'))
       
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
           option_report.append((f'Group "{gname}"', " ".join(cols)))
           
           for col in cols:
              if col not in value_cols:
                  if col in idx_cols:
                      fail(f'Ovelap between index columns and sample columns not permitted; "{col}" is an index column.')
                
                  value_cols.append(col)
    
       group_cols.update(group_dict)
       
    if ref_groups:
        ref_groups = set(ref_groups)
        info(f'Reference sample groups:  {" ".join(ref_groups)}')
        option_report.append(('Ref. sample groups', " ".join(ref_groups)))


        for ref_group in ref_groups:
           if ref_group not in groups:
               fail('Reference sample group "{ref_group}" not found in experimental design. Available group names: {' '.join(sorted(groups))}')
    
    else:
        info(f'No reference groups specified')
        option_report.append(('Reference sample groups', 'None specified'))
    
    info(f'The comparisons will be made between the following pairs of sample groups:')
    option_report.append(('Group pairs compared', None))
    
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
                    option_report.append((f'Pair {len(pairs)}', f'{g1}  vs  {g2}'))
          
        else:
            groups = sorted(groups)
            for i, g1 in enumerate(groups[:-1]):
                for g2 in groups[i+1:]:
                    info(f'   {g1}  -  {g2}')
                    pairs.append((g1, g2))
                    option_report.append((f'Pair {len(pairs)}', f'{g1}  vs  {g2}'))
         
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
    option_report.append((f'Contaminent rows', pre_cull[0]))
    
    if not quiet:
        info('Your data now has the following structure:')
        print(df.head(10))
    
    if pdf:
        font_height = 0.16 # = 11 pt
        wrap_len = 90
        
        n_lines = 0
        y = font_height
        y_offsets = []
        for i, (head, value) in enumerate(option_report):
          y_head = y + font_height
          
          if isinstance(value, str) and len(value) > wrap_len:
              sub_lines = textwrap.wrap(value, width=wrap_len)
              m = len(sub_lines)
              n_lines += m
              option_report[i] = (head, '\n'.join(sub_lines))
              y += font_height * m
          
          else:
              n_lines += 1
              y += font_height
          
          if value is None:
              y += 0.3 * font_height
              y_head += 0.3 * font_height
          
          y_offsets.append((y, y_head))
                 
        fig_height = y_offsets[-1][0] + 2 * font_height
        fig, ax = plt.subplots(figsize=(8.0, fig_height))
        margin = font_height/fig_height
        
        for i, (head, value) in enumerate(option_report):
        
            if value is None:
              text = f'{head}'
              hcolor = '#0080FF'
              
            else:
              hcolor= '#400000'
            
              if isinstance(value, bool):
                  value = 'Yes' if value else 'No'
 
              elif isinstance(value, int):
                  value = f'{value:,}'
 
              elif isinstance(value, float):
                  value = f'{value:.2g}'
              
              text = f'{head} : {value}'
            
            y, y_head = y_offsets[i]
              
            y = 1.0 - (y/fig_height)
            y_head = 1.0 - (y_head/fig_height)
            
            if value:
               ax.text(0.2 - 0.5 * margin, y_head, head + ':', transform=ax.transAxes, color=hcolor, va='bottom', ha='right', fontsize=9.0)
               ax.text(0.2, y, value, transform=ax.transAxes, va='bottom', ha='left', fontsize=9.0)
            else:
                
               if i == 0:
                   ax.text(0.0, y_head, head, transform=ax.transAxes, color=hcolor, va='bottom', ha='left', fontsize=16.0)
                   ax.axhline(y_head, linewidth=1, color='#808080', alpha=0.5)         
                
               else:
                   ax.text(0.0, y_head, head, transform=ax.transAxes, color=hcolor, va='bottom', ha='left', fontsize=11.0)
        
        ax.axhline(margin, linewidth=1, color='#808080', alpha=0.5)         
        ax.axis('off')
                
        plt.subplots_adjust(top=1.0-margin, bottom=margin, left=margin, right=1.0-margin)  
        pdf.savefig(dpi=300)
         
    
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
        
        key = f'{g1}:::{g2}'
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

    for a, b in pairs:
        pair_name = f'{a}:::{b}'
        plots.volcano(pdf, pair_name, plotdict[pair_name], q95dict[pair_name], f_thresh, p_thresh, colors,
                      split_x, hq_only, hit_labels, markers, lw=0.25, ls='--')
    
    if table_path:
        save_volcano_table(pairs, plotdict, table_path, f_thresh, p_thresh)
    
    
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
    
    arg_parse.add_argument('-e', '--experiment-table', dest="e", metavar='FILE_PATH', default=None, 
                           help='The location of a experimental design file. This file is a tab or comma-separated text file containing a table ' \
                                'that relates the name of the data samples (input columns) to which experimental groups they belong. ' \
                                'Each line should first give the name of the input column and then list the groups (or categories) ' \
                                'to which it belongs. Groups will be compared according to their order, i.e. comparisons are only ' \
                                'made within the groups specified each one column of the experiment design table. ' \
                                'Example contents: \n\nSample1\tControl\tFemale\nSample2\tControl\tMale\nSample3\tMutant\tFemale\nSample4\tMutant\tMale')
   
    arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', default=DEFAULT_SOFTWARE,
                           help=f"The name of the software used to process the data present in the input file. Available choices: {' ,'.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}")

    arg_parse.add_argument('-i', '--index-columns', dest="i", metavar='INDEX_COLUMN', default=None,
                           help=f'The names, or numbers starting from 1, of one or more input columns used to uniquely index the data rows. ' \
                                'If not specified columns named will be sought from the default list: {' '.join(DEFAULT_INDEX_COLS)}.')

    arg_parse.add_argument('-r', '--ref-groups', dest="r", metavar='GROUP_NAME', nargs='+',  default=None,
                           help='An optional list of names specifiying which group(s) of samples/columns are considered reference groups.'
                                'Typicaly this would include a control group. These names should match the group/category names in the experiment design file (-e). ' \
                                'When reference groups are present other relevant groups will be compared only to these references, and ' \
                                'not among themselves. If there is no reference group all groups will be compared to all relavent others. ')
    
    arg_parse.add_argument('-m', '--marker-ids', dest="m", metavar='MARKER_ID', nargs='+',  default=None,
                           help='An optional list of marker IDs/accessions to label on plots; must be space separated abd match values in the index columns (-i).')
 
    arg_parse.add_argument('-g', '--graphics-pdf', dest="g", metavar='PDF_FILE_PATH', default=None,
                           help=f"Optional path to save graphs as a PDF file. If not specified graphics will be plotted to screen.")
 
    arg_parse.add_argument('-o', '--out-table', dest="o", metavar='TABLE_FILE_PATH', default=None,
                           help=f"Optional save file path for volcano plot results. The output format is determined by the file extansion; " \
                                "'.xlsx' or '.xls' saves as Excel spreadsheets, '.csv' saves as comma-separated text, and anything else as tab-separated text " \
                                "For text formats, different comparison/volcano plots are saved in separate files, labelled with the groups compared.")

    arg_parse.add_argument('-k', '--ignore-contaminants',dest="k",  default=False, action='store_true',
                           help='Whether to ignore kerating and other contaminants. If NOT set contaminents will be removed.')
    
    arg_parse.add_argument('-c', '--columns', dest="c", metavar='COLUMN_SELECTION', nargs='+', default=None,
                           help='If not already specified within an experimental design file (-e), this option defines column/samples which should be grouped and compared. ' \
                                'This option is designed for use within automated pipelines; for regular users the -e option should be easier. ' \
                                'Column/sample groups are specified in the form "1,2,3:4,5,6:7,8,9  3,4:5,6:10,11", with numbers ' \
                                'starting at 1, identifying the data column. Columns within the same group are joined with commas "," and '\
                                'two or more groups to be compared with one another are joined with colons ":". Several comparisons may be ' \
                                'specified, separated by spaces. Columns numbers may appear in multiple groups/comparisons.') 
 
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

    arg_parse.add_argument('-sx', '--split-x', dest="sx", metavar='XGAP_WIDTH', type=float, default=0.0,
                           help='Optionally split the X-axis (fold change) on volcano plots if the largest gap between data points on ' \
                                'the postive and/or negative side is larger than the value (in axis units) specified here; useful if ' \
                                'you expect a wide range of fold-changes')

    arg_parse.add_argument('--no-labels', dest="no-labels", action='store_true',
                           help='If set, suppress the labelling of significant hits in the volcano plots')
    
    arg_parse.add_argument('--hq-only', dest="hq-only", action='store_true',
                           help='If set, plot only the high-quality points on volcano plots. Otherwise low quality points ' \
                                 '(thous with s single data point in one group) will be plotted, albeit differently. Data points with only ' \
                                 'one value in both compared groups are never plotted.')
    
    args = vars(arg_parse.parse_args(argv))
 
    in_path = args['d']
    pdf_path = args['g']
    table_path = args['o']
    software = args['s']
    remove_contaminents = not args['k']
    idx_cols = args['i']
    ref_groups = args['r']
    markers = args['m']
    columns = args['c']
    exp_path = args['e']
    f_thresh = args['f']
    p_thresh = args['p']
    quiet = args['q']
    log = args['l']
    
    split_x = args['sx']
    
    pos_color = args['pos-color']
    neg_color = args['neg-color']
    low_color = args['low-color']
    hit_labels = not args['no-labels']
    hq_only = args['hq-only']
    
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
    
    lava(in_path, exp_path, software, pdf_path, table_path, idx_cols, ref_groups, markers, columns,
         remove_contaminents, f_thresh, p_thresh, quiet, colors, split_x, hit_labels, hq_only)
  
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # #    
    # Export table mark q95, check q95 rendering on volcano
    # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # 


if __name__ == '__main__':
    
    main()
    
"""
Example: 
python3 lava.py VolcCLI_PD1.xlsx -o VolcCLI_PD1_out.xlsx -g VolcCLI_PD1.pdf -e groups_example.txt -l -r A -m P00918

"""
    

