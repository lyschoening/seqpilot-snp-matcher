# SeqPilot SNP Matcher

Report for SNP calls exported from JSI Sequence Pilot in TXT format.

Generates a PDF report with two tables, one comparing SNPs in the input samples to the reference, one comparing
individual samples to each other, counting mismatches.



	usage: report.py [-h] reference S [S ...] O
	
	Compare a number of JSI Seq Pilot tables

	positional arguments:
		 reference   Reference SNPs, one per line: rs123456<tab>A
		 S           One or more samples to compare
		 O           Output prefix
	
	optional arguments:
	     -h, --help  show this help message and exit


## Requirements

- Python 2.6+
- Jinja 2
