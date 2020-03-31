#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# makeSNPLOC.py
#========================================================================================================

import argparse
import gzip as gz
import csv

class makeSNPLOC:
	def __init__(self, args):
		self.IdField = args.id
		self.VarField = args.var
		self.ChrField = args.chr
		self.BpField = args.bp
		self.InputFil = args.input
		self.OutFil = args.output
	def run(self):
		if self.InputFil.split(".")[-1] == "gz":
			self.reader = csv.reader(gz.open(self.InputFil, 'rt'), delimiter="\t")
		else:
			self.reader = csv.reader(open(self.InputFil, 'rt'), delimiter="\t")
		self.writer = csv.writer(open(self.OutFil, 'wt'), delimiter="\t")
		self.header = next(self.reader)
		#print(self.header)
		if self.VarField != None:
			idx_var = self.header.index(self.VarField)
			idx_id = self.header.index(self.IdField)
			for row in self.reader:
				_id = row[idx_id]
				_var = row[idx_var]
				_chr, _bp, _ref, _alt = _var.split(":")
				self.writer.writerow([_id, _chr, _bp])
		else:
			idx_id = self.header.index(self.IdField)
			idx_chr = self.header.index(self.ChrField)
			idx_bp = self.header.index(self.BpField)
			for row in self.reader:
				_id = row[idx_id]
				_chr = row[idx_chr].strip('chr')
				_bp = row[idx_bp]
				self.writer.writerow([_id, _chr, _bp])

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('--id', required=True, type=str, help = 'String of SNP ID field')
	parser.add_argument('--chr', required=False, type=str, help = 'String of Chromesome field')
	parser.add_argument('--bp', required=False, type=str, help = 'String of BasePair Position field')
	parser.add_argument('--var', required=False, default=None, type=str, help = 'String of BasePair Position field')
	parser.add_argument('-i', '--input', required=True, type=str, help = 'Input GWAS SumStats File Name')
	parser.add_argument('-o', '--output', type=str, help = 'SNP LOC OutPut Name')
	args = parser.parse_args()
	if args.output == None:
		args.output = args.input.split(".")[0] + ".snploc"
	return args

def main():
	args = GetOptions()
	ins = makeSNPLOC(args)
	ins.run()	

	return

if __name__=='__main__':
	main()
