from datetime import datetime
import sys
sys.path.append('/scripts/')

import pandas

from partitioning_HC import HC

print(datetime.now())
hc = HC('test_a')
hc.allele_mx =  pandas.read_csv('/test_data/alleles.tsv', sep='\t')
hc.run()
print(datetime.now())