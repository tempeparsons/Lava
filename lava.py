VALID_SOFTWARE = ('DIANN', 'MaxQuant', 'ProteomeDiscoverer') # Allow "PD", "MQ" etc?
DEFAULT_SOFTWARE = VALID_SOFTWARE[-1]


def lava(in_path, software=VALID_SOFTWARE, remove_contaminents=True):
  
  pass
  

def main(argv=None):

  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]
  
  arg_parse = ArgumentParser(prog='proteome_consec_pep_scan.py', description='Find sequences with the longest subsequences that match consecutive arrays of sub-peptide sequences',
                             epilog=epilog, prefix_chars='-', add_help=True)
  
  # Do we allow multiple (similar) inputs?
  arg_parse.add_argument(metavar='INPUT_FILE', nargs='+', dest='i',
                         help='The file path of the DIA results text file (protein-level, not peptide-level)')

  arg_parse.add_argument('-s', '--software', dest="s", metavar='SOFTWARE', nargs='+', default=None,
                         help='The name of the software used to process the data. Available choices: {' ,'.join(VALID_SOFTWARE)}. Deafult: {DEFAULT_SOFTWARE}')

  arg_parse.add_argument('-k', '--ignore-keratin', default=False, action='store_true',
                         help='Whether to ignore keratin contaminents. If NOT set contaminents will be removed.')
                         
                         
  #arg_parse.add_argument('-rc', '--min-num-consec', default=DEFAULT_N_CONSEC, metavar='PEP_COUNT', type=int, dest="n",
  #                       help=f'Minimum number of consecutive sub-peptide sequences. Default {DEFAULT_N_CONSEC}')

  #arg_parse.add_argument('-o', '--out-file', default=None, metavar='OUT_FILE', dest="o",
  #                       help=f'Optional output file to write results as a tab-separated table')

  args = vars(arg_parse.parse_args(argv))
  
  in_path = args['i']
  software = args['s']
  remove_contaminents = not args['k']
  
  lava(in_path, software, remove_contaminents)
