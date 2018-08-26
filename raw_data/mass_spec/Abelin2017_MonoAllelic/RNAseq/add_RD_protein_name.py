import sqlite3
import pandas as pd

conn = sqlite3.connect("HLA-peptidome-monoallelic-bk.db")
conn.text_factory = str
db = conn.cursor()

with open("PA_ucsc_proteome.fasta", 'r') as f:
        lines = f.readlines()
entries = {}
for i, line in enumerate(lines[:(len(lines)-10)]):
        if str(line.strip()).startswith('>'):
                entries[line.strip()] = ''
                linesToRead = ""
                for j in range(i+1, i+9):
                        if not lines[j].startswith('>'):
                                linesToRead += lines[j].strip("\n")
                        else:
                                break
                entries[line.strip()] = linesToRead

result = db.execute('SELECT peptide FROM ElutionArray WHERE unique_allele = "RD"').fetchall()
column = [field[0] for field in result]
rd_peptides = column

for i, peptide in enumerate(rd_peptides):
    for key, value in entries.iteritems():
        if peptide in value:
            protein_name = key[1:]
            print protein_name
	    db.execute('UPDATE ElutionArray SET protein_name=?  WHERE peptide=?', (protein_name, peptide))
            break

db.execute('SELECT * FROM ElutionArray WHERE unique_allele == "RD"')
print db.fetchone()
conn.commit()
conn.close()

