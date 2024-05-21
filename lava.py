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

VERSION = '0.3.1'

plots.VERSION = VERSION

# Allow shirt codes for software
SOFT_DN = ('DIANN', 'DN')
SOFT_MQ = ('MaxQuant', 'MQ')
SOFT_PD = ('ProteomeDiscoverer', 'PD')
 
VALID_SOFTWARE = SOFT_DN + SOFT_MQ + SOFT_PD
DEFAULT_SOFTWARE = SOFT_PD[0]
DEFAULT_MIN_PEPS = 2

FILE_EXTS_EXL = ('.xls','.xlsx')
FILE_EXTS_TXT = ('.txt',)
FILE_EXTS_CSV = ('.csv','.tsv')

DEFALUT_FC = 2.0
DEFALUT_MAXP  = 0.05
DEFAULT_POS_COLOR = '#D00000'
DEFAULT_NEG_COLOR = '#0080FF'
DEFAULT_LOW_COLOR = '#A0A000'

DEFAULT_INDEX_COLS = ('Accession','Majority protein IDs','Protein.ID','Protein ID')

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
        df = pd.read_table(file_path, na_values=[' '])
  
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

    
def save_volcano_table(pairs, plotdict,  save_path, f_thresh, p_thresh, min_peps):
   
   nmap = {'grp1grp2_FC':'log2_fold_change','Tpvals':'-log2_pvalue',
           'pvals':'p-value', 'nobs_grp1':'nobs_A', 'nobs_grp2':'nobs_B',
           'mean_grp1' : 'log2mean_A', 'mean_grp2': 'log2mean_B',
           'zmean_grp1': 'zmean_A', 'zmean_grp2': 'zmean_B',
           'zstd_grp1': 'zstd_A',  'zstd_grp2': 'zstd_B',           
           }
   
   quad_map = {(True,True):'POS', (True,False):'NEG', (False,True):'fail_pos', (False,False):'fail_neg', }
   
   path_root, file_ext = os.path.splitext(save_path)
   file_ext = file_ext.lower()
   
   keys = [f'{a}:::{b}' for a,b in pairs]
   color_dict = {}
   for key in keys:
       df = plotdict[key]
       lfc = np.array(df['grp1grp2_FC'])
       pvs = np.array(df['Tpvals'])
       n = len(lfc)
       
       if 'npeps' in df:
           npeps = np.array(df['npeps'])
       else:
           npeps = np.full(n, min_peps)

       cats = []
       cat_colors = []
       for i in range(n):
           if npeps[i] < min_peps:
               klass = 'low_pep'
               color = 'background-color: #E0E0E0'
 
           elif pvs[i] >= p_thresh:
               if lfc[i] >= f_thresh:
                  klass = 'POS'
                  color = 'background-color: #FFD0D0'
               elif lfc[i] <= -f_thresh:
                  klass = 'NEG'
                  color = 'background-color: #D0D0FF'
               else:
                  klass = 'fail'
                  color = 'background-color: #FFFFD0'
 
           else:
                klass = 'fail'
                color = 'background-color: #FFFFD0'
 
           cat_colors.append(color)
           cats.append(klass)
       
       df = df.rename(columns=nmap)
       #df.drop(nmap.keys(), axis=1)
       
       j = list(df.columns).index('labels') +1
       
       fc = 2.0 ** lfc
       df.insert(j, 'fold_change', fc)
       df.insert(j, 'hit_class', cats)
       
       sort_cols = [-pvs]
       sort_cols.append(~((lfc <= -f_thresh) & (pvs >= p_thresh)))
       sort_cols.append(~((lfc >= f_thresh) & (pvs >= p_thresh)))
       
       if 'npeps' in df:
         sort_cols.append(np.array(df['npeps'] < min_peps))
       
       sort_idx = np.lexsort(tuple(sort_cols))
       color_dict[key] = [cat_colors[i] for i in sort_idx]
       plotdict[key] = df.iloc[sort_idx]
       
   def grp1_col_bg(vals):
   
     return ['background-color: #EEFFEE'] * len(vals)

   def grp2_col_bg(vals):
   
     return ['background-color: #FFEEFF'] * len(vals)
   
   def color_klass(vals):
     
     return ['background-color: #FFEEFF'] * len(vals)
   
   def color_fc(vals):
     
     styles = ['font-weight: bold' if abs(v) > f_thresh else None for v in vals]
     
     return styles
   
   op_thresh = 2.0 ** (-p_thresh)
   
   def color_pv(vals):
     
     styles = ['font-weight: bold' if v < op_thresh else None for v in vals]
     
     return styles

   def color_class(vals):
     
     styles = []
     
     for v in vals:
         if v > f_thresh:
             styles.append('color: red')
 
         elif v < -f_thresh:
             styles.append('color: blue')
 
         else:
             styles.append(None)
     
     return styles

   def color_nobs(vals):
     
     styles = ['color: red' if v < 2 else None for v in vals]

     return styles

   def color_npep(vals):
     
     styles = ['color: red' if v < min_peps else None for v in vals]

     return styles
     
   if file_ext in FILE_EXTS_EXL:
       
       with pd.ExcelWriter(save_path, engine='xlsxwriter') as writer:
           for key in keys:
               df = plotdict[key]

               in_colsA = [x for x in list(df.columns) if x.startswith('inA_') or x.startswith('npepA_')]
               in_colsB = [x for x in list(df.columns) if x.startswith('inB_') or x.startswith('npepB_')]
               sheet_name=key.replace(':::', '_v_')[:31]
               
               if 'npeps' in df:
                   df.style.apply(grp1_col_bg, axis=1, subset=in_colsA)\
                           .apply(grp2_col_bg, axis=1, subset=in_colsB)\
                           .apply(lambda v:color_dict[key], axis=0, subset=['hit_class', 'fold_change','log2_fold_change', 'p-value', '-log2_pvalue'])\
                           .apply(color_nobs, axis=1, subset=['nobs_A','nobs_B'])\
                           .apply(color_npep, axis=1, subset=['npeps'])\
                           .apply(color_fc, axis=1, subset=['log2_fold_change'])\
                           .apply(color_pv, axis=1, subset=['p-value'])\
                           .to_excel(writer, sheet_name=sheet_name)
               
               else:
                   df.style.apply(grp1_col_bg, axis=1, subset=in_colsA)\
                           .apply(grp2_col_bg, axis=1, subset=in_colsB)\
                           .apply(lambda v:color_dict[key], axis=0, subset=['hit_class', 'fold_change','log2_fold_change', 'p-value', '-log2_pvalue'])\
                           .apply(color_nobs, axis=1, subset=['nobs_A','nobs_B'])\
                           .apply(color_fc, axis=1, subset=['log2_fold_change'])\
                           .apply(color_pv, axis=1, subset=['p-value'])\
                           .to_excel(writer, sheet_name=sheet_name)
               
        
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
    

def _read_label_file(label_path):

    import csv
    
    if '\t' in open(label_path).read():
        delimiter = '\t'
    else:
        delimiter = ','
    
    label_dict = {}
    with open(label_path, newline='') as file_obj:
        for row in csv.reader(file_obj, delimiter=delimiter):
            if not row:
                continue
            
            prot_ids, label = row[:2]
            label = label.strip()
            
            for prot_id in prot_ids.split(';'):
                prot_id = prot_id.strip()
                label_dict[prot_id] = label
            
    return label_dict


def _read_bg_file(df, bg_path, group_dicts, bg_groups):

    import csv
    
    if '\t' in open(bg_path).read():
        delimiter = '\t'
    else:
        delimiter = ','
    
    col_names = set()
    for group_dict in group_dicts:
        for group, cols in group_dict.items():
            col_names.update(cols)
    
    bg_dict = {}
    with open(bg_path, newline='') as file_obj:
        for row in csv.reader(file_obj, delimiter=delimiter):
            if not row:
                continue
            
            if row[0].startswith('#'):
                continue
            
            col_name, bg_name = row[:2]
            
            if col_name not in col_names:
               valid = ', '.join(sorted(col_names))
               fail(f'Sample column name "{col_name}" from background file {bg_path} does not refer to a short name in the experimental design. Valid column names: {valid}')
            
            if bg_name not in bg_groups:
               valid = ', '.join(sorted(bg_groups))
               fail(f'Background name "{bg_name}" from background file {bg_path} does not refer to a group in the experimental design. Valid background groups: {valid}')

            bg_dict[col_name] = bg_name
            
    return bg_dict
            
            
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
    col_rename = None
    
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
            if col_name not in df:
              fail(f'Sample column named {col_name} is not present in input data')
            
            for i, group in enumerate(groups):
                group = group.strip()
                if group:
                    group_dicts[i][group].append(col_name)
                
            
    if group_dicts:
        col_rename = None
        for i, group_dict in enumerate(group_dicts):
             n_single = len([grp for grp in group_dict if len(group_dict[grp]) == 1])
             
             if n_single == len(group_dict):
                 info('Renaming columns according to exp. design file:')
                 col_rename = {group_dict[grp][0]:grp for grp in group_dict}
                 
                 for orig in sorted(col_rename):
                     info(f'  from "{orig}" to "{col_rename[orig]}"')
                 
                 renamed_group_dicts = []
                 for j, group_dict in enumerate(group_dicts):
                     if j == i:
                         continue
                     
                     for grp in group_dict:
                         group_dict[grp] = [col_rename[x] for x in group_dict[grp]]
                    
                     renamed_group_dicts.append(group_dict)
                 
                 group_dicts = renamed_group_dicts
                 df.rename(columns=col_rename, inplace=True)
                 break
                 
        
    else:
        fail(f'Experimental design file {file_path} did not appear to contain anything useful')
    
    bg_groups = {}
    compare_groups = []
    
    for group_dict in group_dicts:
      if len(group_dict) == 1:
        bg_groups.update(group_dict)
      else:
        compare_groups.append(group_dict)
    
    return compare_groups, bg_groups # a list of {group_nameA:[col_name1, col_name2], group_nameB:[col_name3, col_name4],}
    
      
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
     
     if isinstance(selection, str):
         selection = [selection]
      
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

    
def lava(in_path, exp_path=None, bg_path=None, software=DEFAULT_SOFTWARE, pdf_path=None, table_path=None,
         idx_cols=None, ref_groups=None, markers=None, col_groups=None, min_peps=DEFAULT_MIN_PEPS, pep_col=None,
         take_descripts=False, label_file=None, label_cols=None, extra_cols=None, remove_contaminents=True, 
         f_thresh=DEFALUT_FC, p_thresh=DEFALUT_MAXP, quiet=False, colors=(DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_LOW_COLOR),
         split_x=False, hit_labels=False, hq_only=False, znorm_fc=False, quant_scale=False, do_pp=True, do_ff=False,
         intermed_out=False):
    
    in_path = _check_path(in_path, should_exist=True)
    exp_path = _check_path(exp_path, should_exist=True)
    bg_path = _check_path(bg_path, should_exist=False)
    pdf_path = _check_path(pdf_path, should_exist=False)
    table_path = _check_path(table_path, should_exist=False)
    
    if pdf_path:
        pdf = PdfPages(pdf_path)
    else:
        pdf = None
    
    info(f'Lava version {VERSION}')    
    
    cmd_args = [f'"{x}"' if ' ' in x else x for x in sys.argv[1:]]
    
    option_report = [(f'Lava Report version {VERSION}', None),
                     ('Input options', ' '.join(cmd_args)),
                     ('Input data file', in_path),
                     ('Exp. design file', exp_path),
                     ('Background file', bg_path or 'None specified'),
                     ('PDF output file', pdf_path),
                     ('Table output file', table_path),
                     ('Input software', software),
                     ('Min. fold-change', 2.0 ** f_thresh),
                     ('Max. p-value', p_thresh),
                     ('Remove contaminents', remove_contaminents),
                     ('Split volcano X-axis', split_x),
                     ('Show hit labels', hit_labels),
                     ('Show high quality', hq_only),
                     ('Z-norm fold-changes', znorm_fc),
                     ('Quantile spot size', quant_scale),
                     ]
        
    color_pos, color_neg, color_low = colors
    df =  _read_data(in_path)
    cols = list(df.columns)
    n_cols = len(cols)
    group_dicts = None
    extra_id_col = None
    
    if exp_path:  # Read from exp design file if possible
        group_dicts, bg_groups = _read_exp_file(df, exp_path)
   
    elif col_groups:
        group_dicts = _split_col_groups(df, col_groups)
        bg_groups = {}
        
    if bg_path:
        bg_dict = _read_bg_file(df, bg_path, group_dicts, bg_groups)
    else:
        bg_dict = {}

    if label_file:
        label_dict = _read_label_file(label_file)
    else:
        label_dict = {}
    
    # Remove any non-comparison groups, e.g. backgrounds
        
    if idx_cols:
        idx_cols = _check_cols(df, idx_cols)
        
        if not idx_cols:
            fail('No valid index columns found')
        
        for col in DEFAULT_INDEX_COLS:
            if (col in df) and (col not in idx_cols):
                extra_id_col = col
                break
                
    else:
        idx_cols = [x for x in DEFAULT_INDEX_COLS if x in df]
      
        if not idx_cols:
            fail(f'No valid index columns found after falling back to default; {' '.join(DEFAULT_INDEX_COLS)}')
        
        extra_id_col = None
            
    if not idx_cols:
        fail('No index column specified')

    elif not group_dicts:
        fail('No groups specified')
    
    if pep_col:
        if pep_col not in df:
            fail(f'Peptide count column "{pep_col}" not found in input file')
        else:
            info(f'Using peptide count column "{pep_col}"')
    
    else:
       min_peps = 1
        
    for idx_col in idx_cols:
        info(f'Using index column {cols.index(idx_col)+1} : "{idx_col}"')
    
    option_report.append(('Index column(s)', ' '.join(idx_cols)))
    
    if take_descripts:
        info(f'Descriptions will be kept for output')

    option_report.append(('Keep descriptions', take_descripts))
    
    if label_dict:
        idx = list(df[idx_cols[0]])
        n_missing = 0
        labels = []
        
        for i, x in enumerate(idx):
            if x in label_dict:
                label = label_dict[x]
            else:
                n_missing += 1
                #label = x
                label = f'Q{n_missing}'
                
            labels.append(label)
       
        if n_missing:
            warn(f'A total of {n_missing:,} rows were missing labels in file {label_file}')
        
        label_cols = [f'From file "label_file"']
        
    elif label_cols:
        label_cols = _check_cols(df, label_cols)
        info(f'Specified label columns ' + '+'.join(label_cols))
    
        if len(label_cols) == 1:
            labels = list(df[label_cols[0]])
        else:
            vals = [[str(x) for x in df[col]] for col in label_cols]
            labels = [':'.join(x) for x in zip(*vals)]
        
    else:
        if 'Gene Symbol' in df:
            labels = list(df['Gene Symbol'])
            label_cols = ['Gene Symbol']

        elif 'Gene Symbol_1' in df:
            labels = list(df['Gene Symbol_1'])
            label_cols = ['Gene Symbol_1']
 
        elif 'Fasta headers' in df:
            gene_names = []
            for i, head in enumerate(df['Fasta headers']):
                if not isinstance(head, str):
                    gene_names.append(df[idx_cols[0]][i])
                    continue
 
                gene_name = re.findall('GN=(\S+)', head)
 
                if not gene_name:
                    gene_name = re.findall('\|(\w+)\|', head)
 
                if gene_name:
                    gene_names.append(';'.join(sorted(set(gene_name))))
 
                else:
                    gene_names.append(df[idx_cols[0]][i])

            labels = gene_names
            label_cols = ['Fasta GN']
 
        elif 'Description' in df:
            gene_names = []
 
            for i, desc in enumerate(df['Description']):
                gene_name = []
                match = re.search('GN=(\S+)', desc)
 
                if match:
                    gene_name = match.group(1)
                else:
                    gene_name = df[idx_cols[0]][i]
 
                gene_names.append(gene_name)
 
            labels = gene_names
            label_cols = ['Description GN']
 
        else:
            for col in DEFAULT_INDEX_COLS:
                if col in df:
                    break
                
            else:
                col = idx_cols[0] 
           
            labels = list(df[col])
            label_cols = [col]
    
    labels = [df[idx_cols[0]][i] if x in (np.nan,'',None) else x for i, x in enumerate(labels)]
    labels = [f'{x[:16]}...' if len(x) > 16 else x for x in labels]
    df['labels'] = labels
 
    if extra_cols:
        extra_cols = _check_cols(df, extra_cols)
    else:
        extra_cols = []
    
    option_report.append(('Label column', '+'.join(label_cols)))
    option_report.append(('Label file', label_file or 'None specified'))
    option_report.append(('Peptide column', pep_col or 'None specified'))
    option_report.append(('Min hit peptides', min_peps))
    
    # Descriptions
    
    if take_descripts:
        if 'Description' in df:
           df['description'] = list(df['Description'])
        
        elif 'Fasta headers' in df:
           descriptions = []
           for i, head in enumerate(df['Fasta headers']):
               if (not head) or (head is np.nan):
                 descriptions.append('None')
                 continue
               
               description = re.findall('\|\S+\s+(.+)\s+OS=', head)
               
               if description:
                   descriptions.append(description[0])
               else:
                   descriptions.append('None')

           df['description'] = descriptions
    
    # Markers
    
    if markers:
        avail_idx = {}
        for idx_col in idx_cols:
            for idx, label in zip(df[idx_col], df['labels']):
                avail_idx[idx] = label

        valid = []    
        for marker in markers:
           marker = marker.strip(',')
           
           if marker not in avail_idx:
               warn(f'Marker "{marker}" is not present in index')
           else:
               valid.append(avail_idx[marker])
        
        marker_text = ' '.join([f'{a};{b}' for a, b in zip(markers, valid)])
        markers = valid
        info(f'Markers:  {marker_text}')
        option_report.append(('Markers', marker_text))

    df.set_index(idx_cols, inplace=True)

    msg = 'Sample group membership'
    info(msg)
    option_report.append((msg, None))
    value_cols = []
    bg_cols = []
    for bg_grp, cols in bg_groups.items():
      bg_cols += cols
      
    bg_cols = sorted(set(bg_cols))
    
    groups = set()
    group_cols = {}
    for i, group_dict in enumerate(group_dicts):
       info(f'Comparison {i+1}:')
       option_report.append(('Comparison', f'{i+1}'))
       
       for gname in sorted(group_dict):
           cols = group_dict[gname]
           
           #if gname in group_cols:
           #    fail(f'Sample group "{gname}" is repeated in a different context; group names should only be present in one experimental design column.')
           
           if len(cols) < 2:
               warn(f'Sample group "{gname}" contains only one sample; cannot be used')
               del group_dict[gname]
               continue
               
       
           groups.add(gname)
           info(f'    Group "{gname}": {", ".join(cols)}')
           option_report.append((f'Group "{gname}"', ", ".join(cols)))
           
           for col in cols:
              if col not in value_cols:
                  if col in idx_cols:
                      fail(f'Overlap between index columns and sample columns not permitted; "{col}" is an index column.')
                
                  value_cols.append(col)
    
       group_cols.update(group_dict)
    
    if bg_groups:
        names = sorted(bg_groups)
        info(f'Background groups:  {", ".join(names)}')
        option_report.append(('Background groups', ", ".join(names)))
    
    else:
        info(f'No background groups specified')
        option_report.append(('Background groups', 'None specified'))
         
       
    if ref_groups:
        ref_groups = set(ref_groups)
        info(f'Ref. sample groups:  {", ".join(ref_groups)}')
        option_report.append(('Ref. sample groups', ", ".join(ref_groups)))

        for ref_group in ref_groups:
           if ref_group not in groups:
               avail = ' '.join(sorted(groups))
               fail(f'Reference sample group "{ref_group}" not found in experimental design. Available group names: {avail}')
    
    else:
        ref_groups = set()
        info(f'No reference groups specified')
        option_report.append(('Ref. sample groups', 'None specified'))
    
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
                    msg = f'{g1}  vs  {g2}'
                    info('   ' + msg)
                    pairs.append((g1, g2))
                    option_report.append((f'Pair {len(pairs)}', msg))
          
        else:
            groups = sorted(groups)
            for i, g1 in enumerate(groups[:-1]):
                for g2 in groups[i+1:]:
                    msg = f'{g1}  vs  {g2}'
                    info('   ' + msg)
                    pairs.append((g1, g2))
                    option_report.append((f'Pair {len(pairs)}', msg))
    
    if bg_dict:
         bgs = set(bg_dict.values())
         info(f'The following background corrections will be employed:')
         for bg in bgs:
            col_names = [x for x in bg_dict if bg_dict[x] == bg]
            col_names = ', '.join(col_names)
            info(f'   With background {bg}: samples {col_names}')
    
    pre_cull = df.shape
    
    if remove_contaminents and software in SOFT_DN:
        df = df[df['First.Protein.Description'].str.contains('Keratin') == False]
 
    if remove_contaminents and software in SOFT_MQ:
        if 'Majority protein IDs' in idx_cols:
          df = df[df.index.str.startswith('CON_') == False]
          df = df[df.index.str.startswith('REV_') == False]
        
        else:
          df = df[df['Majority protein IDs'].str.startswith('CON_') == False]
          df = df[df['Majority protein IDs'].str.startswith('REV_') == False]
  
    if remove_contaminents and software in SOFT_PD:
        if 'Description' in df:
            df = df[df['Description'].str.contains('Keratin') == False]
            df = df[df['Description'].str.contains('Trypsin') == False]

    if remove_contaminents:
        for col in DEFAULT_INDEX_COLS:
            if col == idx_cols[0] and len(idx_cols) == 1:
                df = df[df.index.str.startswith('cRAP') == False]
            elif col in df:
                df = df[df[col].str.startswith('cRAP') == False]
    
    post_cull = df.shape

    culled_rows = pre_cull[0] - post_cull[0]
    info(f'Removed {culled_rows} contaminant rows from {pre_cull[0]}, total')
    option_report.append((f'Input rows', pre_cull[0]))
    option_report.append((f'Contaminent rows', culled_rows))
    option_report.append((f'Useful rows', post_cull[0]))
    
    if pep_col and min_peps > 1:
       n1 = df.shape[0]
       nlow = np.count_nonzero(df[pep_col] < min_peps)
       info(f'Found {nlow} low peptide (<{min_peps}) rows from {n1:,} total')
       option_report.append((f'Low peptide rows', nlow))
       
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
        plots._watermark(fig, f'Lava v{VERSION} by H.T. Parsons')
        
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
    
    orig_df = df.copy()
    df = np.log2(df[value_cols + bg_cols])
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
    
    plotdict = {}
    overFCthresh = {}
    
    # Z normalise value distributions; used in T-test and optionally for fold-change

    norm_cols = util.make_znorm(df, value_cols, bg_dict, bg_groups) # Now keeps original cols, now only done once
    
    if bg_dict:
        info('Plotting corrected correlations')
        plots.correlation_plot(df, norm_cols, pdf, title="Corrected sample correlation matrix")
        
    info('Plotting PCA')
    if bg_dict:
        plots.plot_dual_pca(df, value_cols, norm_cols, group_cols, pdf)    
    
    else:
        plots.plot_pca(df, value_cols, group_cols, pdf)    
    
    if znorm_fc or bg_dict:
      info('Plotting distribution normalization')
      plots.plot_norm(df, value_cols, norm_cols, pdf)
      
    info(f'Using log2 fold-change threshold: {f_thresh}')
    for g1, g2 in pairs:
        key = f'{g1}:::{g2}'
        grp1 = group_cols[g1]
        grp2 = group_cols[g2]
      
        # Columns for Z normalise value distributions; used in T-test and optionally for fold-change
        zgrp1 = ['znorm_' + x for x in grp1]
        zgrp2 = ['znorm_' + x for x in grp2]
        
        df2 = df[grp1+grp2+zgrp1+zgrp2].copy()
                
        # Add labels and original, untransformed data to output
        df2['labels'] = orig_df['labels']
        
        for col1 in grp1:
            df2['inA_' + col1] = orig_df[col1]
            
            if intermed_out:
                df2['log2_A_' + col1] = df[ col1]
                df2['normA_' + col1] = df['znorm_' + col1]
 
                if bg_cols:
                    df2['corrA_' + col1] = df['bgcorr_' +col1]

        for col2 in grp2:
            df2['inB_' + col2] = orig_df[col2]
 
            if intermed_out:
                df2['log2_B_' + col2] = df[ col2]
                df2['normB_' + col2] = df['znorm_' + col2]

                if bg_cols:
                    df2['corrB_' + col2] = df['bgcorr_' + col2]
        
        if extra_id_col:
            df2['protein_id'] = orig_df[extra_id_col]
        
        if pep_col:
            df2['npeps'] = orig_df[pep_col]
        
        # Count non-Nan observations
        df2['nobs_grp1'] = df2.loc[:,grp1].count(axis=1)
        df2['nobs_grp2'] = df2.loc[:,grp2].count(axis=1)
        
        # Remove rows , for this sample pair, with insuffcicient data
        df2 = util.remove_rows(df2, len(grp1), len(grp2)) # Now only removes rows doesn't change zeros
         
        # Identify single and all-zero rows before any modifications, zero-filling etc.
        azeros1 = df2['nobs_grp1'] == 0
        azeros2 = df2['nobs_grp2'] == 0
        singles1 = df2['nobs_grp1'] == 1
        singles2 = df2['nobs_grp2'] == 1
        
        # Any pure zeros are real zeros not NaN; set before taking means
        df2.loc[azeros1, grp1] = 0
        df2.loc[azeros2, grp2] = 0
        df2.loc[azeros1, zgrp1] = 0
        df2.loc[azeros2, zgrp2] = 0
        
        # Any singlular data has nobs set to 2 for t-test
        
        # Means ignoring sparse NaNs
        df2['mean_grp1'] = df2.loc[:,grp1].mean(axis=1)
        df2['mean_grp2'] = df2.loc[:,grp2].mean(axis=1)
        df2['zmean_grp1'] = df2.loc[:,zgrp1].mean(axis=1)
        df2['zmean_grp2'] = df2.loc[:,zgrp2].mean(axis=1)
        
        # Get minimum non-zero/non-nan means; these are used to replace zeros for FC calc
        min_mu1 = np.quantile(np.array(df2['mean_grp1'])[~azeros1], 0.05)
        min_mu2 = np.quantile(np.array(df2['mean_grp2'])[~azeros2], 0.05)
        min_zmu1 = np.quantile(np.array(df2['zmean_grp1'])[~azeros1], 0.05)
        min_zmu2 = np.quantile(np.array(df2['zmean_grp1'])[~azeros2], 0.05)
        
        # Replace all-zero means with a lower bound
        df2['mean_grp1'] = np.where(azeros1, min_mu1, df2['mean_grp1'])
        df2['mean_grp2'] = np.where(azeros2, min_mu2, df2['mean_grp2'])
        df2['zmean_grp1'] = np.where(azeros1, min_zmu1, df2['zmean_grp1'])
        df2['zmean_grp2'] = np.where(azeros2, min_zmu2, df2['zmean_grp2'])
        
        # Calc fold changes
        if znorm_fc:
            df2['grp1grp2_FC'] = df2['zmean_grp1'] - df2['zmean_grp2']
        else:
            df2['grp1grp2_FC'] = df2['mean_grp1'] - df2['mean_grp2']          
       
        # Get Z-normed params for T-test 
        df2['zstd_grp1'] = df2.loc[:,zgrp1].std(axis=1, ddof=1)
        df2['zstd_grp2'] = df2.loc[:,zgrp2].std(axis=1, ddof=1)
        df2['znobs_grp1'] = df2.loc[:,zgrp1].count(axis=1) # Counts clipped zeros
        df2['znobs_grp2'] = df2.loc[:,zgrp2].count(axis=1)
        
        #q95grp1 = df2.zstd_grp1.quantile(q=0.95)
        #q95grp2 = df2.zstd_grp2.quantile(q=0.95)
        
        # Current heristic is to average the other side with a quantile STD replacement
        zstds1 = 0.5 * (df2.zstd_grp1.quantile(q=0.75) + df2['zstd_grp2'])
        zstds2 = 0.5 * (df2.zstd_grp2.quantile(q=0.75) + df2['zstd_grp1'])
        
        # Set replacement std fro all-zeros and singlular
        df2['zstd_grp1'] = np.where(azeros1, zstds1, df2['zstd_grp1'])
        df2['zstd_grp2'] = np.where(azeros2, zstds2, df2['zstd_grp2'])
        df2['zstd_grp1'] = np.where(singles1, zstds1, df2['zstd_grp1'])
        df2['zstd_grp2'] = np.where(singles2, zstds2, df2['zstd_grp2'])

        # Singlular counts are artificially set to 2 and STD to q95 for T-test
        df2.loc[singles1, 'znobs_grp1'] = 2
        df2.loc[singles2, 'znobs_grp2'] = 2
                
        pos, neg = util.prpn_log2FCs_over_threshold(df2, f_thresh)
        info(f'  - pair {g1} vs {g2} proteins over theshold: pos; {pos}% neg; {neg}%')
        
        df2 = util.ttest_from_stats_eqvar(df2)
        df2 = df2.drop(columns = zgrp1 + zgrp2 + grp1 + grp2)
       
        col_selection = ['protein_id'] if extra_id_col else []
        
        if 'description' in orig_df:
            df2['description'] = orig_df['description']
            col_selection += ['description']
            
        # Reorder
        col_selection += ['labels','grp1grp2_FC', 'pvals', 'Tpvals', 'nobs_grp1', 'nobs_grp2',
                          'mean_grp1', 'mean_grp2', 'zmean_grp1', 'zmean_grp2', 'zstd_grp1', 'zstd_grp2']
             
        if pep_col:
            col_selection.append('npeps')
         
        # Orig data last
        col_selection += [c for c in df2.columns if c.startswith('inA_')]
        col_selection += [c for c in df2.columns if c.startswith('inB_')]
        
        if intermed_out:
             col_selection += [c for c in df2.columns if c.startswith('log2_A_')]
             col_selection += [c for c in df2.columns if c.startswith('log2_B_')]
 
             if bg_cols:
                 col_selection += [c for c in df2.columns if c.startswith('corrA_')]
                 col_selection += [c for c in df2.columns if c.startswith('corrB_')]

             col_selection += [c for c in df2.columns if c.startswith('normA_')]
             col_selection += [c for c in df2.columns if c.startswith('normB_')]
        
        for col in extra_cols:
            df2[col] = orig_df[col]
            col_selection.append(col)
         
        plotdict[key] = df2[col_selection]

    plots.pvalFC_hists(plotdict, pdf)

    for a, b in pairs:
        pair_name = f'{a}:::{b}'
        
        df =  plotdict[pair_name]
        
        """
        fig, (ax0, ax1, ax2) = plt.subplots(1, 3)
        
        mu1 = np.array(df['zmean_grp1'])
        mu2 = np.array(df['zmean_grp2'])
        sig1 = np.array(df['zstd_grp1'])
        sig2 = np.array(df['zstd_grp2'])
        
        hist1, edges1 = np.histogram(mu1[mu1 > 0], bins=50)
        ax0.plot(edges1[:-1], hist1, color='#0080FF', alpha=0.5)
        hist2, edges2 = np.histogram(mu2[mu2 > 0], bins=50)
        ax0.plot(edges2[:-1], hist2, color='#B0B000', alpha=0.5)
         
        ax1.scatter(mu1, mu2, color='#0080FF', alpha=0.25, s=3)
        ax1.set_xlabel('Mean 1')
        ax1.set_ylabel('Mean 2')
        
        ax1.plot([mu1.min(), mu1.max()], [mu2.min(), mu2.max()], color='#808080', linewidth=0.5)
         
        x_vals = np.maximum(mu1, mu2)
        y_vals = np.minimum(sig1, sig2)
        ax2.scatter(x_vals, y_vals, color='#FF0000', alpha=0.5, s=3)
        ax2.set_xlabel('Mean')
        ax2.set_ylabel('STD')
        
        plt.show()
        """

        plots.volcano(pdf, pair_name, plotdict[pair_name], f_thresh, p_thresh, min_peps,
                      colors, quant_scale, split_x, hq_only, hit_labels, markers, lw=0.25, ls='--',
                      marker_text_col=None)
   
    if (do_pp or do_ff) and len(pairs) > 1:
        for i, pair1 in enumerate(pairs[:-1]):
            key1 = ':::'.join(pair1)
            s1 = set(pair1)
            for pair2 in pairs[i+1:]:
                s2 = set(pair2)
                
                if len(pairs) > 2  and (not s1 & s2):
                  continue
                
                key2 = ':::'.join(pair2)

                if do_pp:
                    plots.dual_comparison_plot(pdf, pair1, plotdict[key1], pair2, plotdict[key2],
                                            p_thresh, f_thresh, min_peps, markers)
                if do_ff:
                    plots.dual_comparison_plot(pdf, pair1, plotdict[key1], pair2, plotdict[key2],
                                               p_thresh, f_thresh, min_peps, markers, is_pp=False)

    if table_path:
        save_volcano_table(pairs, plotdict, table_path, f_thresh, p_thresh, min_peps)
    
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
    arg_parse.add_argument(metavar='INPUT_FILE', dest='in',
                           help='The input file path to read data from; the DIA results text file (protein-level, not peptide-level)')
    
    arg_parse.add_argument('-e', '--experiment-table', dest="e", metavar='FILE_PATH', default=None, 
                           help='The experimental design file; this is a comma- or tab-separated text file containing a table, ' \
                                'relating the names of data samples to their experimental groups/categories; labels like "Control", "CondA" etc. ' \
                                'In the experimental design table the sample names (matching the input data file) should be in the first column. ' \
                                'Optionally, a column of short/new names for the samples may be specified; these must be unique for each row.' \
                                'Subsequently, the groups (or categories) to which each samples belong are specified in further columns. ' \
                                'Each grouping column defines a distinct set of comparisons, i.e. samples can be compared according to different variables.' \
                                'A column wiith only a single group is assumed to contain a grouping of background data, as used with the -b option.')

    arg_parse.add_argument('-b', '--background-table', dest="b", metavar='FILE_PATH', default=None, 
                           help='A file specifying background/baseline "subtractions" to perform; so comparisons can be made relative to a set of background data, which may differ for different conditions. ' \
                                'This file is a comma- or tab-separated table with two columns that pair each experimental sample column (to be analysed) with it\'s ' \
                                'background group. These names must be present in the experimental design table.')
   
    arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', default=DEFAULT_SOFTWARE,
                           help=f"The name of the software used to process the data present in the input file. Available choices: {', '.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}")
    
    default_idx = ' '.join(DEFAULT_INDEX_COLS)
    arg_parse.add_argument('-i', '--index-columns', dest="i", metavar='INDEX_COLUMN', default=None, nargs='+',
                           help=f'The names, or numbers starting from 1, of one or more input columns used to uniquely index the data rows. If not specified, column names will be sought from the default list: {default_idx}.')

    arg_parse.add_argument('-n', '--min-peps', dest="n", metavar='PEP_COUNT_COLUMN', default=DEFAULT_MIN_PEPS, type=int, 
                           help=f'The minimum number of theoretical peptides required for a protein to be considered as a hit. '\
                                  'No selection will be done unless the peptide column is specified via the -nc option.'\
                                  'Proteins that fail will be reported separately.')

    arg_parse.add_argument('-nc', '--num-peps-column', dest="nc", metavar='PEP_COUNT_COLUMN', default=None,
                           help=f'The name, or number starting from 1, of the input data column containing peptide counts. Proteins may be filtered according to this column using the -n option.')

    arg_parse.add_argument('-r', '--ref-groups', dest="r", metavar='GROUP_NAME', nargs='+',  default=None,
                           help='An optional list of names specifiying which group(s) of samples/columns are considered reference groups. '
                                'Typicaly this would include any control group. These names should match the group/category names in the experiment design file (-e). ' \
                                'When reference groups are present other relevant groups will be compared only to these references, and ' \
                                'not among themselves. If there is no reference group then all groups will be compared to all relavent others. ')
    
    arg_parse.add_argument('-m', '--marker-ids', dest="m", metavar='MARKER_ID', nargs='+',  default=None,
                           help='An optional list of marker IDs/accessions to label on plots; these must be space separated and match values in the index columns (-i).')
 
    arg_parse.add_argument('-g', '--graphics-pdf', dest="g", metavar='PDF_FILE_PATH', default=None,
                           help=f"Optional path to save graphs as a PDF file. If not specified graphics will be plotted to screen.")
 
    arg_parse.add_argument('-o', '--out-table', dest="o", metavar='TABLE_FILE_PATH', default=None,
                           help=f"Optional save file path for volcano plot results. The output format is determined by the file extension; " \
                                "'.xlsx' or '.xls' saves as Excel spreadsheets, '.csv' saves as comma-separated text, and anything else as tab-separated text. " \
                                "For simple textual formats, different comparison/volcano reults are saved in separate files, labelled according to the groups compared.")

    arg_parse.add_argument('-k', '--ignore-contaminants',dest="k",  default=False, action='store_true',
                           help='Whether to ignore keratin and other contaminants. If NOT set contaminents will be removed.')
    
    arg_parse.add_argument('-c', '--columns', dest="c", metavar='COLUMN_SELECTION', nargs='+', default=None,
                           help='If not already specified within an experimental design file (-e), this option defines column/samples which should be grouped and compared. ' \
                                'This option is designed for use within automated pipelines; for regular users using the -e option should be easier. ' \
                                'Column/sample groups are specified in the form "1,2,3:4,5,6:7,8,9  3,4:5,6:10,11", with numbers ' \
                                'starting at 1, identifying the data column. Columns within the same group are joined with commas "," and '\
                                'two or more groups to be compared with one another are joined with colons ":". Several sets of comparisons may be ' \
                                'specified, separated by spaces. Columns numbers may appear in multiple groups/comparisons.') 
 
    arg_parse.add_argument('-z', '--z-norm-fc', dest="z", action='store_true',
                           help=f"When set, fold changes will be calculated from Z-normalised values.")

    arg_parse.add_argument('-d', '--descriptions', dest="d", action='store_true',
                           help=f"If present, cary over any protein description information (inc. from FASTA headers) to output tables.")

    arg_parse.add_argument('-v', '--intermed-values-out', dest="v", action='store_true',
                           help=f"If set, write all intermediate log-transformed values, normalized values etc. to output tables.")

    arg_parse.add_argument('-l', '--label-columns', dest="l", nargs='+', default=None,
                           help='The column names from the input to use as labels for data points; defaults to use gene names, if available, or else protein IDs')

    arg_parse.add_argument('-lf', '--label-file', dest="lf", default=None,
                           help="A separate file from which to take protein/row labels, i.e. where the main input file doesn't contain an appropriate label column. " \
                                "This file should be a comma- or tab-separated text file with two columns; the first column contains one or more IDs, delimited by semicolons, " \
                                "(typically protein accessions or unique identifiers to index data rows) " \
                                "and the second column giving a label that maps to those IDs. Using this option overrides the -l option.")

    arg_parse.add_argument('-x', '--extra-columns', dest="x", nargs='+', default=None,
                           help='Extra column names from the input to carry across (unchanged) to the output tables')

    arg_parse.add_argument('-pp', '--pvalue-dual', dest="pp", action='store_true',
                           help='When two or more separate comparisons are done, create dual plots of p-values from different sources.')
    
    arg_parse.add_argument('-ff', '--fold0change-dual', dest="ff", action='store_true',
                           help='When two or more separate comparisons are done, create dual plot of fold-changes from different sources.')

    arg_parse.add_argument('-qs', '--quantile-scaling', dest="qs", action='store_true',
                           help='Use quantile scaling of abundance to determine spot size in the volcno plots, otherwise the (maximum) mean abundance is used directly.')
                           
    arg_parse.add_argument('-lg', '--log-status', dest="lg", action='store_true',
                           help=f"When set, writes status information to a log file; the log file path will be based on the input path.")
    
    arg_parse.add_argument('-q', '--quiet-mode', dest="q", action='store_true',
                           help='Proceed quietly, without user interaction; implicitly accept defaults to user queries')
    
    arg_parse.add_argument('-f', '--min-fold-change', dest="f", metavar='MIN_RATIO', type=float, default=DEFALUT_FC,
                           help=f'Minimum fold change for potentially significant hits. Default: {DEFALUT_FC}')
   
    arg_parse.add_argument('-p', '--max-p-value', dest="p", metavar='MAX_PVALUE', type=float, default=DEFALUT_MAXP,
                           help=f'Maximum threshold p-value for selecting significant hits. Default: {DEFALUT_MAXP}')
    
    arg_parse.add_argument('--pos-color', dest="pos-color", metavar='COLOR', default=DEFAULT_POS_COLOR,
                           help=f'Optional color specification (used by matplotlib) for positive hits on volcano plots, e.g. "blue" or "#0000FF". Default: {DEFAULT_POS_COLOR}')
    
    arg_parse.add_argument('--neg-color', dest="neg-color", metavar='COLOR', default=DEFAULT_NEG_COLOR,
                           help=f'Optional color specification (used by matplotlib) for negative hits on volcano plots, e.g. "red" or "#FF0000". Default: {DEFAULT_NEG_COLOR}')
                           
    arg_parse.add_argument('--low-color', dest="low-color", metavar='COLOR', default=DEFAULT_LOW_COLOR,
                           help=f'Optional color specification (used by matplotlib) for insignificant points on volcano plots, e.g. "grey" or "#808080". Default: {DEFAULT_LOW_COLOR}')

    arg_parse.add_argument('-xs', '--split-x', dest="xs", metavar='XGAP_WIDTH', type=float, default=0.0,
                           help='Optionally split the X-axis (fold change) on volcano plots if the largest gap between data points on ' \
                                'the postive and/or negative side is larger than the value (in axis units) specified here; useful if ' \
                                'you expect a wide range of fold-changes')

    arg_parse.add_argument('--no-labels', dest="no-labels", action='store_true',
                           help='If set, suppress the labelling of significant hits in the volcano plots')
    
    arg_parse.add_argument('-hq', '--hq-only', dest="hq-only", action='store_true',
                           help='If set, plot only the high-quality points on volcano plots. Otherwise low quality points ' \
                                 '(thous with s single data point in one group) will be plotted, albeit differently. Data points with only ' \
                                 'one value in both compared groups are never plotted.')
    
    args = vars(arg_parse.parse_args(argv))
 
    in_path = args['in']
    pdf_path = args['g']
    table_path = args['o']
    software = args['s']
    remove_contaminents = not args['k']
    idx_cols = args['i']
    ref_groups = args['r']
    markers = args['m']
    columns = args['c']
    exp_path = args['e']
    bg_path = args['b']
    f_thresh = args['f']
    p_thresh = args['p']
    quiet = args['q']
    log = args['lg']
    znorm_fc = args['z']
    min_peps = args['n']
    pep_col = args['nc']
    quant_scale = args['qs']
    take_descripts = args['d']
    label_cols = args['l']
    label_file = args['lf']
    extra_cols = args['x']
    intermed_out = args['v']
    do_pp = args['pp']
    do_ff = args['ff']
    
    split_x = args['xs']
    
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
    
    if label_file and label_cols:
        warn('Both "-lf" (label file) and "-l" (label column) options were specified; the latter will be ignored.')
        label_cols = None
    
    if f_thresh <= 1.0:
        fail('Minimum fold change threshold (-f) must be greeater than 1.0')
    
    f_thresh = np.log2(f_thresh)   
        
    if p_thresh > 0.5:
        fail('Maximum p-value threshold must be < 0.5')
    
    if p_thresh < 0.0:
        fail('Maximum p-value threshold must be positive')
    
    _check_color(pos_color)
    _check_color(neg_color)
    _check_color(low_color)
    
    colors = (pos_color, low_color, neg_color)

    lava(in_path, exp_path, bg_path, software, pdf_path, table_path, idx_cols, ref_groups, markers,
         columns, min_peps, pep_col, take_descripts, label_file, label_cols, extra_cols,
         remove_contaminents, f_thresh, p_thresh, quiet, colors, split_x,
         hit_labels, hq_only, znorm_fc, quant_scale, do_pp, do_ff, intermed_out)


# Excel sheet name 31 char limit

if __name__ == '__main__':
    
    main()
    
