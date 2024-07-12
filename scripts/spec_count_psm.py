import sys, re

from math import log
from collections import defaultdict

def spec_count_psm(file_path, out_file_path, prot_file_path=None, 
                   acc_col="Master Protein Accessions", desc_col="Master Protein Descriptions",
                   file_col="File ID", prot_abun_cols='Abundance (F\d+) Sample \w+', prot_acc_col='Accession'):
  
  abun_dict = {}
  
  if prot_file_path:
    print(f'Protein abundances from {prot_file_path}')   
    
    with open(prot_file_path) as file_obj:
      head = file_obj.readline().rstrip('\n').split('\t')
      head = [x.strip('"') for x in head]
      
      cols = []
      file_head_dict = {}
      j = head.index(prot_acc_col)
      
      for i, col in enumerate(head):
        file_num = re.findall(prot_abun_cols, col)
        
        if file_num:
          cols.append((i, file_num[0]))
      
      abun_dict = {}
      for line in file_obj:
         row = [x.strip('"') for x in line.rstrip('\n').split('\t')]
         acc = row[j]
         abundances = {key:row[i].strip() for i, key in cols}
         abun_dict[acc] = abundances
          
    
  print(f'Reading spectra from {file_path}')   
  with open(file_path) as file_obj:
    head = file_obj.readline().rstrip('\n').split('\t')
    head = [x.strip('"') for x in head]
    
    col_acc = head.index(acc_col)
    col_desc = head.index(desc_col)
    col_file = head.index(file_col)
    
    ddict = {}
    cdict = {}
    
    for line in file_obj:
      row = [x.strip('"') for x in line.rstrip('\n').split('\t')]
      acc = row[col_acc]
      
      if not acc:
        continue
      
      if acc.startswith('cRAP'):
        continue
      
      desc = row[col_desc]
      file_id = row[col_file]
      
      if file_id not in cdict:
        cdict[file_id] = defaultdict(int)
      
      ddict[acc] = desc
      cdict[file_id][acc] += 1
    
  
  file_ids = sorted(cdict)
  prot_ids = sorted(ddict)
  n_lines = 0
  n_missing_prots = 0
  
  with open(out_file_path, 'w') as out_file_obj:
    
    head = ['prot_ids','gene_name',]
    head += file_ids
    
    if abun_dict:
      head += [f'Ab:{x}' for x in file_ids]
      head += [f'log_2:{x}' for x in file_ids]
    
    head += ['protein_description',]
    out_file_obj.write('\t'.join(head) + '\n')
    
    for pid in prot_ids:
      
      desc = ddict[pid]
      gene_name = re.findall('GN=(\S+)', desc)
      
      if gene_name:
        gene_name = ';'.join(sorted(set(gene_name)))
      else:
        gene_name = '?'
      
      file_data = '\t'.join(['%d' % cdict[f][pid] for f in file_ids])
      
      if abun_dict:
        if pid not in abun_dict:
          n_missing_prots += 1
          continue
        
        abunds = [abun_dict[pid][f] for f in file_ids]
        abun_data = '\t'.join([x if x else '0' for x in abunds])

        log_abunds = ['%.3f' % log(float(x), 2.0) if x else '0' for x in abunds]
        log_abun_data = '\t'.join(log_abunds)
      
        line = f'{pid}\t{gene_name}\t{file_data}\t{abun_data}\t{log_abun_data}\t{desc}\n'
      
      else:
        line = f'{pid}\t{gene_name}\t{file_data}\t{desc}\n'
  
      out_file_obj.write(line)
      n_lines += 1
      
  if abun_dict:
     print(f'Missing proteins compated to PSM file: {n_missing_prots:,}')   
  
  print(f'Write {n_lines:,} lines to {out_file_path}')   
      
      
if __name__ == '__main__':
   
   if len(sys.argv) == 4:
     in_file_path = sys.argv[1]
     prot_file_path = sys.argv[2]
     out_file_path = sys.argv[3]
   
   else:
     in_file_path = sys.argv[1]
     out_file_path = sys.argv[2]
     prot_file_path = None
   
   spec_count_psm(in_file_path, out_file_path, prot_file_path)


"""
python3 spec_count_psm.py PSMs.txt spec_counts.tsv

python3 spec_count_psm.py PSMs.txt Proteins.txt abundances_spec_counts.tsv


"""
