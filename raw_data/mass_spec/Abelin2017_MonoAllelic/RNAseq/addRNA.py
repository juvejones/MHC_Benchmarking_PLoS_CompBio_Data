import sqlite3
import pandas as pd

conn = sqlite3.connect("HLA-peptidome-monoallelic-bk.db")
conn.text_factory = str
db = conn.cursor()

files = {'A2902': 'GSM2450855_A29_02.isoforms.results.txt',
		 'B5101': 'GSM2450856_B51_01.isoforms.results.txt',
		 'B5401': 'GSM2450857_B54_01.isoforms.results.txt'}

#db.execute('CREATE TABLE TRANSCRIPT (HLAtype text, transcript_id text, gene_id int, length int, effective_length real, tpm real, fpkm real, isopct real)')
#db.execute('CREATE UNIQUE INDEX transcript_idx ON TRANSCRIPT(HLAtype, transcript_id)')

for key, value in files.iteritems():
	allele = key
	df = pd.read_csv(value, header=0, sep="\t", index_col=False)
	print df[:5]
	for idx, row in df.iterrows():
		v = (allele, row[0],row[1],row[2],row[3],row[5],row[6],row[7])
		db.execute('INSERT OR REPLACE INTO TRANSCRIPT VALUES (?,?,?,?,?,?,?,?)', v)

db.execute('SELECT * FROM TRANSCRIPT')
print(len(db.fetchall()))

conn.commit()
conn.close()

