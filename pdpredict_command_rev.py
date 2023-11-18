#!/opt/websites/anaconda/envs/hari_37/bin/python
'''#!/usr/bin/env python3.7'''
import cgi, os
from functools import total_ordering
import sys
import time
import shutil
import subprocess
import math
import numpy as np
from os import path
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import uuid
import cgitb;
cgitb.enable()
import timeit
start = timeit.default_timer()

do=input("Enter the type of input: 'pdb-id' or 'pdb-file'")
if do=="pdb-id":
	inn=input("Enter the PDB-ID and protein/DNA/chains (optional) as 1a4t-A-B or 1a4t--B or 1a4t-A")
	inn=inn.split("-")
	if len(inn)>1:
		pdb_id=inn[0]
		chain=inn[1]
	if (len(inn)>2):
		dchain=inn[2]
	else:
		pdb_id=inn[0]
		chain=""
		dchain=""
	pdb_id=pdb_id.lower()
if do=='pdb-file':
	inn=input("Enter the PDB file name in current folder")
	pdb_file=inn
	pdb_id=""
	chain=""
	dchain=""
rna_strand=input("Please select if the RNA is single stranded or double stranded (type ds or ss)")
'''
form = cgi.FieldStorage()
pdb_id=form.getvalue('fname')
chain=form.getvalue('chain')
dchain=form.getvalue('lchain')

if form['pdbf'].filename:          
	fileitem=form['pdbf'].filename
else:
'''
fileitem=b""
#rna_strand=form.getvalue('rna_strand')
pdb_id=pdb_id.lower()
#import sys
sys.stdout.flush()
print ('Content-Type: text/html\r\n')
#print('/r/n')
#os.system("rm -r pd_res_*")
if __name__ == '__main__':
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	#print(pdb_id)
	#os.system("chmod +x pdpredict.py")
	print('<div style="height:400px">')
	print('</div>')
	print('<div style="vertical-align: middle;">')
	print("<center><div style='width:1000px;height:200px;border:solid black; font-size:20px;size = '+4''>")
	print("Feature calculations...")
	print("<br>")
	print("<center><div style='width:1000px;height:100px;border-top:solid black;font-size:15px;size = '+2''>")
	print('Please wait until the calculations are done')
	print('<br>')
	'''
	redirectURL = "/cgi-bin/pdpredict.py"
	print ('Content-Type: text/html\r\n',flush=True)
	print ('<html>',flush=True)
	print ('  <head>',flush=True)
	print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL,flush=True)
	#print ('	<title>You are going to be redirected</title>')
	#print ("Use this url to obtain data: <a href='{}'>Result_page</a>".format(redirectURL)
	#print ('Content-Type: text/html\r\n')
	#print(pdb_id+'/r/n')
	#print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
	#exit()
	g=open("index.txt").readlines()
	for gg in g:
		print(gg,flush=True)
		#gg=gg.rstrip()
		#print ("""{}""").format(gg)
		#print("""print (\"\"\"{}\"\"\")""".format(gg))
	print("calculating features....",flush=True)
	#print(fileitem,flush=True)
	#print(pdb_id,flush=True)
	h=open("footer.txt").readlines()
	for hh in h:
		#hh=hh.rstrip()
		print (hh,flush=True)
		#print("""print (\"\"\"{}\"\"\")""".format(hh))
		#print("\n")
	'''
	#os.system("find -mmin +20 -type d -exec rmdir !('tmp') > files_older_last_deleted 2>&1")
	#find dir_* -mmin +34 -type d -exec rm -r {} \;
	#if chain=="":
	#	print("chain,dchain")
	#pdb_id='1aay'
	#chain='A'
	#dchain='B'
	#fileitem=''
	#rna_strand='ss'
	#print(chain,dchain,pdb_id,fileitem,rna_strand)
	#model='alpha'
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	if pdb_id!="" and chain!="" and dchain!="" and do=='pdb-id' and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=dchain
		method=1
	elif pdb_id!="" and chain!="" and dchain=="" and do=='pdb-id' and rna_strand!="": 
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=""
		method=1
	elif pdb_id!="" and chain=="" and dchain!="" and do=='pdb-id' and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=dchain
		method=1
	elif pdb_id!="" and chain=="" and dchain=="" and do=='pdb-id' and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=""
		print("in")
		method=1
	elif pdb_id==""and do=='pdb-id':
		print ('Content-Type: text/html\r\n')
		#print(fileitem,flush=True)
		print('/r/n')
		print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
		exit()
	elif rna_strand=="":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Input error</title><body>Select DNA strand</body></html>')
		exit()
	elif pdb_id=="" and do=='pdb-file':
		method=2
		#print(fileitem,flush=True)
		#print (fileitem.filename)
		#if form['pdbf'].filename:
		#	fn = os.path.basename(form['pdbf'].filename)
		#	open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	elif pdb_id=="" and chain!="" and dchain=="" and do=='pdb-file':
		#fileitem = form['pdbf']
		chain=chain
		#lchain=form.getvalue('chain')
		method=2
		#print (fileitem.filename)
		#if form['pdbf'].filename:
		#	fn = os.path.basename(form['pdbf'].filename)
		#	open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	else:
		#print ('Content-Type: text/html\r\n')
		#print('/r/n')
		print ('<html><head><title>Input error</title><body>No input found</body></html>')
		exit()
	'''	
	if pdb_id!="" and chain!="" and dchain!="" and fileitem==b"" and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=dchain
		method=1
	elif pdb_id!="" and chain!="" and dchain=="" and fileitem==b""and rna_strand!="": 
		pdb_id_up=pdb_id.upper()
		chain=chain
		dchain=""
		method=1
	elif pdb_id!="" and chain=="" and dchain!="" and fileitem==b""and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=dchain
		method=1
	elif pdb_id!="" and chain=="" and dchain=="" and fileitem==b""and rna_strand!="":
		pdb_id_up=pdb_id.upper()
		chain=""
		dchain=""
		#print("in")
		method=1
	elif pdb_id==""and fileitem=="":
		#print ('Content-Type: text/html\r\n')
		#print(fileitem,flush=True)
		#print('/r/n')
		print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
		exit()
	elif rna_strand=="":
		#print ('Content-Type: text/html\r\n')
		#print('/r/n')
		print ('<html><head><title>Input error</title><body>Select RNA strand</body></html>')
		exit()
	elif pdb_id=="" and fileitem!="":
		method=2
		#print(fileitem,flush=True)
		#print (fileitem.filename)
		if form['pdbf'].filename:
			fn = os.path.basename(form['pdbf'].filename)
			open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	elif pdb_id=="" and chain!="" and dchain=="" and form.getvalue('pdbf')!=b"":
		#fileitem = form['pdbf']
		chain=chain
		lchain=form.getvalue('chain')
		method=2
		#print (fileitem.filename)
		if form['pdbf'].filename:
			fn = os.path.basename(form['pdbf'].filename)
			open('tmp/' + fn, 'wb').write(form['pdbf'].file.read())
	else:
		#print ('Content-Type: text/html\r\n')
		#print('/r/n')
		print ('<html><head><title>Input error</title><body>No input found</body></html>')
		exit()
		model=input("Please enter the structural classification of protein (all-alpha/all-beta/alpha-beta/not-known)")
		if model=="all-alpha":
			model="alpha"
		if model=="all-beta":
			model="beta"
		if model=="alpha-beta":
			model="alphabeta"
			model1=input("Please enter the Functional classification of protein (Enzyme/Regulatory/other)")
			if model1=="Regulatory" or model1=="regulatory":
				model1="reg"
			if model1=="Enzyme" or model1=="Enzyme":
				model1="enz"
			if model1=="other" or model1=="others":
				model1="str"
		if model=="not-known":
			model="none"
	'''
	if rna_strand=='ds':		
		model1=input("Please enter the Functional classification of protein (Enzyme/Regulatory/other)")
		if model1=="Regulatory" or model1=="regulatory":
			model1="reg"
		if model1=="Enzyme" or model1=="enzyme":
			model1="enz"
		if model1=="other" or model1=="others":
			model1="str"
		model=""
		
	if rna_strand=='ss':
		model1=""
		model=input("Please enter the Functional classification of protein (Enzyme/Regulatory/other)")
		if model=="Regulatory" or model=="regulatory":
			model="reg"
		if model=="other" or model=="others":
			model="str"
	#model=""
	#model1=""
	model2=rna_strand
	#model=""	#model2='ss' 
	'''
	if rna_strand=='ds':
		model1=form.getvalue('complex_func') # model='beta' 
		model=form.getvalue('complex_type')
		#print(model) #model1='reg'
		#print(model1)
	if form.getvalue('RNA_strand')=='ss':
		model1=""
		model=""
	#model=""
	#model1=""
	model2=form.getvalue('RNA_strand')
	#model=""	#model2='ss' 
	#model1=""
	#print(model2)
	'''
	'''
	if model=="None" and model2!="ss":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the structural classfication </title><body>Please select the structural classfication of the protein</body></html>')
		exit()
	if model1=="None" and model2!="ss":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the Functional classfication </title><body>Please select the Functional classfication of the protein </body></html>')
		exit()
	if model2=="None":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Please select the  dna strand</title><body>please select type of DNA</body></html>')
		exit()
	'''
	if model1=='reg':
		modeltype="Regulatory"
	if model1=='enz':
		modeltype="Enzyme"
	elif model1=='str':
		modeltype="Others"
	if model2=='ds':
		modeltype2="Double strand"
		modeltype1=''
	elif model2=='ss':
		modeltype2="Single strand"
		modeltype=''
	if model=='reg':
		modeltype1="Regulatory"
	elif model=='str':
		modeltype1="Others"
	
#########################################################################################################	Single strand	 ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
#########################################################################################################	Single strand	  ###########################################################################################################
	
	import shutil
	#import wget
	import glob
	import Bio.PDB as bpdb
	from Bio.PDB import is_aa
	from Bio.PDB import PDBParser, PDBIO, Select
	import urllib
	import os
	import numpy as np
	import re
	import pandas as pd
	import math
	from Bio import PDB
	import warnings
	from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
	pdb_id_up=pdb_id
	if method==1:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname="pr_res_"+uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.system("chmod -R 777 {}".format(path))
		os.system("chmod -777 bin")
		os.system("chmod +x bin")
		
		os.chdir(path)
		#print(os.listdir())
		#print(path)
		#shutil.copyfile("../1aay.pdb", "1aay.pdb")
		#os.system("wget 'https://files.rcsb.org/download/{}.pdb'".format(pdb_id_up))
		#print(pdb_id_up)
		os.system("wget 'https://files.rcsb.org/download/{}.pdb'".format(pdb_id_up))
		try:
			os.system("wget 'https://files.rcsb.org/download/{}.pdb1'".format(pdb_id_up))
			shutil.copyfile(pdb_id_up+".pdb1", pdb_id_up+".pdb")
		except:
			print("n")
		#filename=wget.download('https://files.rcsb.org/download/{}.pdb1'.format(pdb_id_up))
		#filename1=wget.download('https://files.rcsb.org/download/{}.pdb'.format(pdb_id_up))
		
	elif method==2:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname="pr_res_"+uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.chdir(path)
		os.system("mv ../tmp/"+fn+" "+fn)
		pdb_id_up=fn.split(".")[0]
		os.system("scp "+fn+" "+pdb_id_up+".pdb1")
		'''Process PDB'''
	class ProtSelect(Select):
		warnings.simplefilter('ignore', PDBConstructionWarning)
		warnings.simplefilter('ignore', FutureWarning)
		def accept_residue(self, residue):
			if not is_aa(residue, standard=True):
				res = residue.id[0]
				if not res == "W":
					return True 
			else:
				return False
	class ProtSelect1(Select):
		def accept_residue(self, residue):
			warnings.simplefilter('ignore', PDBConstructionWarning)
			warnings.simplefilter('ignore', FutureWarning)
			if is_aa(residue, standard=True):
				return True 
			else:
				return False
	
	print('Page will be automatically redirected to result page once the prediction is done.')
	print('<br>')
	print('Job ID:'+randname)
	
	parser = PDBParser()
	structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb")
	modell = structure[0]
	io = bpdb.PDBIO()
	io.set_structure(modell)
	modd=0
	io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
	io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
	with open('dna_'+pdb_id_up+'.pdb') as filedna:
		with open('apo_'+pdb_id_up+'.pdb') as protfile:
			if len(filedna.readlines())<5 or len(protfile.readlines())<5:
				for modd in range(1,len(structure)):
					modell = structure[modd]
					io = bpdb.PDBIO()
					io.set_structure(modell)
					io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
					io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
					if not len(filedna.readlines())<5 or not len(protfile.readlines())<5:
						break
	
	if dchain!="":
		dchains=dchain.split(",")
		class ChainSelect(Select):
			def __init__(self, schain):
			    self.schain = schain
			def accept_schain(self, schain):
			    if schain.get_id() in self.schain:
			    	return 1
			    else:
			    	return 0
		p = PDBParser(PERMISSIVE=1)       
		structure = p.get_structure('dna'+pdb_id_up+'.pdb', 'dna_'+pdb_id_up+'.pdb')                               
		io_w_no_h = PDBIO()               
		io_w_no_h.set_structure(structure)
		io_w_no_h.save('dna_'+pdb_id_up+'.pdb', ChainSelect(dchains))	
	if chain!="":
		chains=chain.split(",")
		class ChainSelect(Select):
			def __init__(self, schain):
			    self.schain = schain
			def accept_schain(self, schain):
			    if schain.get_id() in self.schain:
			    	return 1
			    else:
			    	return 0
		p = PDBParser(PERMISSIVE=1)       
		structure = p.get_structure('apo_'+pdb_id_up+'.pdb', 'apo_'+pdb_id_up+'.pdb')
	                                
		io_w_no_h = PDBIO()               
		io_w_no_h.set_structure(structure)
		io_w_no_h.save('apo_'+pdb_id_up+'.pdb', ChainSelect(chains))		
	if dchain!="" or chain!="":	
		chains=chain.split(",")
		dchains=dchain.split(",")
		allchain=chains+dchains
		print(allchain)
		class ChainSelect(Select):
			def __init__(self, schain):
			    self.schain = schain
			def accept_schain(self, schain):
			    if schain.get_id() in self.schain:
			    	return 1
			    else:
			    	return 0
		p = PDBParser(PERMISSIVE=1)       
		structure = p.get_structure(pdb_id_up+'.pdb', pdb_id_up+'.pdb')
	                                
		io_w_no_h = PDBIO()               
		io_w_no_h.set_structure(structure)
		io_w_no_h.save(pdb_id_up+'.pdb', ChainSelect(allchain))

	shutil.copyfile("../foldx", "foldx")
	shutil.copyfile("../rotabase.txt", "rotabase.txt")
	shutil.copyfile("../naccess", "naccess")
	shutil.copyfile("../bin","bin")
	shutil.copyfile("../style4.css", "style4.css")
	shutil.copyfile("../clean_pdb.py", "clean_pdb.py")
	#os.system(r"python3 clean_pdb.py {}.pdb".format(pdb_id_up))
	shutil.copyfile("../index.txt", "index.txt")
	shutil.copyfile("../footer.txt", "footer.txt")
	shutil.copyfile("../dssp","dssp")
	shutil.copyfile("../rna_4vdrch.csv", "rna_4vdrch.csv")
	shutil.copyfile("../aa_20vdrch.csv", "aa_20vdrch.csv")
	shutil.copyfile("../potential_res.csv","potential_res.csv")
	shutil.copyfile("../hbplus","hbplus")
	os.system("chmod -R +x {}".format("foldx"))
	os.system("chmod -R +x {}".format("hbplus"))
	os.system("chmod -R +x {}".format("naccess"))
	os.system("chmod -777 {}".format("naccess"))
	os.system("chmod -R +x {}".format("bin"))
	os.system("chmod -R +x {}".format("dssp"))
	os.system("chmod -R 777 {}".format(path))
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	print(model2)
	print(model)
	if model2=='ss' and model=='reg':	
		print("</div></center></div>")	
		print("<html>")
		print("<div style='width:100px;height:1px;overflow:hidden;'>")
		with open("result.txt","w") as resultout:
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[modd]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['A','G','C','U']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						#return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						#return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									print("")
									if (distance<= 3.5):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
										print("")
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
						#df_inter=pd.DataFrame(int_df1)
					else:
						df_inter=pd.DataFrame(int_df1)
					print("")
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						print("")
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
						print("")
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print("")
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict

			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('rna_4vdrch.csv')
			oc_count=[]
			bind=[]
			nc_count=[]	
			so_count=[]
			tot=[]
			dna_n_atoms=[]
			dna_p_atoms=[]
			ele1=[]
			nonpolar,polar,charged=[],[],[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
					sys.stdout.flush()
					print("")
					if len(df_inter.columns)>3:
						atoms_involve = inst2.f8_interaction_type(df_inter)
						print("")
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						
						if 'OC' in df_at.index:
							oc_count.append(df_at.loc[['OC']].values)
						tot.append(df_at['count'].sum())	
						if 'NC' in df_at.index:
							nc_count.append(df_at.loc[['NC']].values)
						if 'SO' in df_at.index:
							so_count.append(df_at.loc[['SO']].values)
						for ao in df_at.index:
							if ao[1]=='N':
								print(ao)
								dna_n_atoms.append(int(df_at.loc[[ao]].values))#+df_at.loc[['NC']].values+df_at.loc[['NN']].values+df_at.loc[['NP']].values)
						for ao in df_at.index:
							if ao[1]=='P':
								print(ao)
								dna_p_atoms.append(int(df_at.loc[[ao]].values))#+df_at.loc[['NC']].values+df_at.loc[['NN']].values+df_at.loc[['NP']].values)
							
						for p0 in df_inter['prot_resno'].tolist():
							bind.append(i+"_"+str(int(p0)))
						aa={'non-polar':['GLY','ALA','VAL','LEU','MET','ILE','PHE','TYR','TRP'],'polar':['SER','THR','CYS','PRO','ASN','GLN'],'charged':['LYS','ARG','HIS','ASP','GLU']} 
						df_c=df_inter[['prot_res','prot_resno']].drop_duplicates(ignore_index=True)
						df_c.drop_duplicates(ignore_index=True)
						df_c=df_c.replace(aa['non-polar'],'non-polar')
						df_c=df_c.replace(aa['polar'],'polar')
						df_c=df_c.replace(aa['charged'],'charged')
						print("")
						df_new=df_c.groupby(df_c.iloc[:,0])["prot_res"].count()
						#print(df_new.values)
						#x=df_new.iloc[0:3].to_list()
						if "non-polar" in df_new.index:
							nonpolar.append(df_new.loc[['non-polar']].values)
						if "charged" in df_new.index:
							charged.append(df_new.loc[['charged']].values)
						if "polar" in df_new.index:
							polar.append(df_new.loc[['polar']].values)
						#charged.append(x[0])
						#nonpolar.append(x[1])
						#polar.append(x[2])
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						ele1.append(energy_df['ele_energy'].sum())
			dna_p_atoms=sum(dna_p_atoms)
			dna_n_atoms=sum(dna_n_atoms)			
			bind=list(set(bind))
			if len(bind)==0:
				print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
				print("Binding site residues are not found. Kindly check the protein/RNA chains!!!")
				exit()
				print("</div>")
			print("</div></center></div>")	
			print("<html>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")					
			if len(nonpolar)>0:
				nonpolar=list(sum(nonpolar))
			else:
				nonpolar=[[0]]
			if len(polar)>0:
				polar=list(sum(polar))
			else:
				polar=[[0]]
			if len(charged)>0:
				charged=list(sum(charged))
			else:
				charged=[[0]]
			sys.stdout.flush()
			print ('Content-Type: text/html\r\n')
			print('/r/n')
			single_tilt,single_twist=[],[]
			method_bind=1			
			import urllib.request 
			k=0
			c=0
			l=1
			try:
				urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
				with open(pdb_id_up+"_bp_step.pars") as file:
					for j in file.readlines():
						if "***local step parameters***" in j:
							k=1
						if k==1:
							if not "#" in j:
								j=j.strip()
								#print(j)
								s=[str for str in j.split(" ") if str.strip()]
								single_tilt.append(float(s[4]))
								
								l=0
					if l==1:
						print(i+"error:single_base_step")
					single_tilt=(sum(single_tilt)/len(single_tilt))
			except:	
					try:
						l=1
						c=-1
						k=0
						hrise=[]
						#shutil.copytree("../x3dna-v2.4/bin", "x3dna-v2.4/bin")
						
						#os.system("source ~/.bashrc")
						#os.chdir("../x3dna-v2.4/bin")
						
						#os.system(">ppoo")
						#os.chdir("../bin")
						#os.system("chmod -R a+X bin")
						#os.system("chmod +X find_pair")
						
						#os.system("./find_pair "+pdb_id_up+".pdb stdout | analyze stdin", shell=True)
						#p4=subprocess.Popen("/var/www/html/bioinfo2/cgi-bin/prapred/x3dna-v2.4/bin/find_pair "+pdb_id_up+".pdb stdout | analyze stdin",shell=True)
						#p4.wait()	
						#p4=subprocess.Popen("./analyze "+pdb_id_up+".bps",shell=True)
						#p4.wait()
						
						with open("bp_step.par") as file:
							for j in file.readlines():
								if "***local step parameters***" in j:
									k=1
								if k==1:
									if not "#" in j:
										j=j.strip()
										#print(j)
										s=[str for str in j.split(" ") if str.strip()]
										single_tilt.append(float(s[4]))
										
										l=0
							if l==1:
								print(i+"error:single_base_step")
							single_tilt=(sum(single_tilt)/len(single_tilt))
						print('other_class')	
						os.chdir("../../"+randname)
					except:
						print("</div>")
						print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
						print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
						print("</div>")
						exit()
			print("</div></center></div>")	
			print("<html>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")	
			os.chdir("../foldx5Linux64.tar_/")
			print(os.listdir)
			p2=subprocess.Popen("./foldx_20231231 --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithRNA=True",  shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			print(fold_dict['Number of Residues'])
			energy_Ionisation=fold_dict['energy Ionisation']
			entropy_sidechain=fold_dict['entropy sidechain']
			Interface_Residues_VdW_Clashing=fold_dict['Interface Residues VdW Clashing']
			Number_of_Residues=fold_dict['Number of Residues']
			os.chdir("../"+randname)
			'''
			
			try:
				import urllib.request
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				hrise=[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "helical" in j:
							k=1
							hrise.append([])
							c=c+1
						if k==1:
							if "ave" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								hrise[c].append(float(s[3]))
								l=0
					if l==1:
						print(i+"error:single_base_step")
					hrise1=sum(hrise[0])/len(hrise[0])
					print('other_class')
			except:
				try:
					import urllib.request
					urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
					l=1
					c=-1
					k=0
					hrise=[]
					with open(pdb_id_up+"_summary.txt") as file:	
						for j in file.readlines():
							if "helical" in j:
								k=1
								hrise.append([])
								c=c+1
							if k==1:
								if "ave" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									hrise[c].append(float(s[3]))
									l=0
						if l==1:
							print(i+"error:single_base_step")
						hrise1=sum(hrise[0])/len(hrise[0])
						print('other_class')					
				except:	
					try:
						l=1
						c=-1
						k=0
						hrise=[]
						#shutil.copytree("../x3dna-v2.4/bin", "x3dna-v2.4/bin")
						os.chdir("../x3dna-v2.4/bin")
						
						p4=subprocess.Popen("./find_pair "+pdb_id_up+" "+pdb_id_up+".bps",shell=True)
						p4.wait()	
						p4=subprocess.Popen("./analyze "+pdb_id_up+".bps",shell=True)
						p4.wait()
						
						with open(pdb_id_up+".out") as file:	
							for j in file.readlines():
								if "helical" in j:
									k=1
									hrise.append([])
									c=c+1
								if k==1:
									if "ave" in j:
										j=j.strip()
										s=[str for str in j.split(" ") if str.strip()]
										hrise[c].append(float(s[3]))
										l=0
							if l==1:
								print(i+"error:single_base_step")
							hrise1=sum(hrise[0])/len(hrise[0])
						print('other_class')	
					except:
						print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
						print("Error in step parameter calcualtion. Please check if the RNA is double stranded")
						print("</div>")
						exit()	
			'''
			#fold_vdw=fold_dict['Van der Waals']
			#print(fold_energy_ion)
			#print(fold_vdw)
			#print(elec)
			'''
			with open(r"apo_bind_"+pdb_id_up+".pdb","w") as file1:
					with open(r"apo_"+pdb_id_up+".pdb") as file:
						for i in file.readlines():
							#print(i[0:26])
							if i[0:4]=="ATOM":
								chnu=i[21:26].split()
								cha=chnu[0]
								number=chnu[1]
								#print(chain)
								#print(number)
								print(cha+"_"+number)
								print(bind[0:3])
								if cha+"_"+number in bind:
									file1.write(i)
			'''
			with open(r"apo_bind_"+pdb_id_up+".pdb","w") as file1:
				with open(r"apo_"+pdb_id_up+".pdb") as file:
					for i in file.readlines():
						if i[0:4]=="ATOM":
							print(i[21:26])
							print("\n")
							try:
								chain=i[21:26].split()[0]
								number=i[21:26].split()[1]
							except:
								chain=i[21].strip()
								number=i[22:26].strip()
							print(number)
							if chain+"_"+number in bind:
								file1.write(i)
			#os.chdir(randname)
			#shutil.copyfile("../hbplus","hbplus")
			#os.system("sudo chmod -R +x {}".format("hbplus"))

			#os.chdir("../")
			#os.chdir(randname)
			p4=subprocess.Popen("./hbplus "+ pdb_id_up+".pdb", shell=True)
			p4.wait()
			k=1
			ss1=0
			with open(pdb_id_up+".hb2") as file:
				for j in file.readlines():
					if "dist" in j:
						k=0
					if k==0:
						#print(i.split())
						#print(j)
						j=j.split()
						if j[5]=="SS":
							ss1=ss1+1
			
			'''
			os.system("chmod -R 777 apo_bind"+pdb_id_up+".pdb")
			#os.system("chmod -R 777 apo_1po6.pdb")
			shutil.copy("apo_bind_"+pdb_id_up+".pdb", "../vossvolvox-master/xyzr/apo_bind_"+pdb_id_up+".pdb")
			os.chdir("../")
			os.chdir("vossvolvox-master/xyzr")
			os.system("chmod +x pdb_to_xyzr")
			p4=subprocess.Popen("./pdb_to_xyzr apo_bind_"+pdb_id_up+".pdb > apo_bind_"+pdb_id_up+".xyzr",shell=True)
			p4.wait()
			os.system("chmod -R 777 apo_bind_"+pdb_id_up+".xyzr")
			shutil.copy("apo_bind_"+pdb_id_up+".xyzr","../apo_bind_"+pdb_id_up+".xyzr")
			os.chdir("../")
			p4=subprocess.Popen("./bin -i apo_bind_"+pdb_id_up+".xyzr > volume_final3",shell=True)
			p4.wait()
			#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
			#p4.wait()
			with open("volume_final3") as file:
				voll=file.readlines()[0].split()
				vol1=voll[2].strip()
				surface1=voll[3].strip()
			
			print(vol1)
			print(surface1)
			os.chdir(path)
			'''
			#os.chdir(randname)
			surface1=1000
			os.chdir("../"+randname)
			#solv=fold_dict['Solvation Hydrophobic']
			#co1=float(sum(co_count))/float(sum(tot))
			os.system("chmod -R 777 apo_bind_"+pdb_id_up+".pdb")
			#os.system("chmod -R 777 apo_1po6.pdb")
			shutil.copy("apo_bind_"+pdb_id_up+".pdb", "../vossvolvox-master/xyzr/apo_bind_"+pdb_id_up+".pdb")
			os.chdir("../")
			os.chdir("vossvolvox-master/xyzr")
			os.system("chmod +x pdb_to_xyzr")
			p4=subprocess.Popen("./pdb_to_xyzr apo_bind_"+pdb_id_up+".pdb > apo_bind_"+pdb_id_up+".xyzr",shell=True)
			p4.wait()
			os.system("chmod -R 777 apo_bind_"+pdb_id_up+".xyzr")
			shutil.copy("apo_bind_"+pdb_id_up+".xyzr","../apo_bind_"+pdb_id_up+".xyzr")
			os.chdir("../")
			p4=subprocess.Popen("./bin -i apo_bind_"+pdb_id_up+".xyzr > volume_final3",shell=True)
			p4.wait()
			#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
			#p4.wait()
			with open("volume_final3") as file:
				voll=file.readlines()[0].split()
				vol1=voll[2].strip()
				surface1=voll[3].strip()
			print(vol1)
			print(surface1)
			os.chdir(path)
			print(os.listdir)
			#print(nn_count)
			if len(polar)>0:
				polar=polar[0]
			#print(nonpolar)
			#NN=[[0]]
			#if len(nn_count)>0:
			#	NN=nn_count[0]
			#print(NN[0][0])
			OC=sum(oc_count)/sum(tot)
			NC=sum(nc_count)/sum(tot)
			SO=sum(so_count)/sum(tot)
			#protein_residue_interface=nonpolar+polar+charged
			#print(protein_residue_interface[0])
			#ss=0
			#n_atoms1=sum(n_atoms)
			ele=sum(ele1)
			from Bio.PDB import PDBParser
			from Bio.PDB.DSSP import DSSP
			p = PDBParser()
			structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
			modelp = structure[modd]
			dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
			#print(dssp)
			from Bio.PDB.DSSP import dssp_dict_from_pdb_file
			dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
			dssp_dict = dssp_tuple[0]
			print(" ")
			dssp_sec=[]
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			print(dssp_dict)
			bind=list(set(bind))
			print(bind)
			for ds in bind:
				dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
			print(dssp_sec)
			dssp=dssp_sec
			he=dssp.count("helix")/len(dssp)
			be=dssp.count("Beta")/len(dssp)
			co=dssp.count("coil")/len(dssp)
			alpha=he*100
			
			#print(NN[0][0],polar,protein_residue_interface[0],back_bone_clash,OP1,surface1,ss1,fold_energy_ion)
			#predval=-5.4669651636618335 *float(ele)+0.34827567021621914 *float(NC)+-0.24907652931359803 *float(polar)+-0.17648013322187248 *float(n_atoms1)+-0.3737410656854529 *float(hrise1)+0.09477647996191321 *float(OC)+0.022299515902453937 *float(interface_residues)+1.0393093258925301 *float(SO)+-7.217379736779106		
			predval=-0.17007889722986055 *float(entropy_sidechain)+-0.7750572118726267 *float(dna_p_atoms)+5.788596848549408 *float(energy_Ionisation)+0.24251937214243677 *float(Interface_Residues_VdW_Clashing)+-0.025537493245634882 *float(dna_n_atoms)+-0.025264276360160137*float(single_tilt)+-0.02202547481601388 *float(alpha)+-0.0033632929277268575 *float(Number_of_Residues)+-6.6902225631777625
			predval="%.2f" % predval
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval))
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
##########################################################################################					Double-all-alpha-reg		  ##########################################################################################
##########################################################################################					Double-all-alpha-reg			   ##########################################################################################
##########################################################################################					Double-all-alpha-reg			   ##########################################################################################
##########################################################################################					Double-all-alpha-reg			   ##########################################################################################
	if model2=='ss' and model=='str':
		with open("result.txt","w") as resultout:
				print("</div></center></div>")	
				print("<html>")
				print("<div style='width:100px;height:1px;overflow:hidden;'>")	
				class Protein_RNA_ineractions:
					def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
						self.pdb_file = pdb_file
						self.prot_chain = prot_chain
						self.rna_chain = rna_chain
						self.pattern ='^ATOM.{16}'
					def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
						warnings.simplefilter('ignore', PDBConstructionWarning)
						warnings.simplefilter('ignore', FutureWarning)
						int_df1 = pd.DataFrame()
						flag = 0
						parser = PDB.PDBParser()
						structure = parser.get_structure("pdb", pdb)
						model = structure[modd]
						prot_chain = model[self.prot_chain] 
						rna_chain = model[self.rna_chain]
						nal1 = ['A','G','C','U']
						pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
						c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
						c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
						checkp = set(c1).intersection(set(nal1))
						checkn = set(c2).intersection(set(pal1))
						if (len(checkp) > 0):
							print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
							flag = 1
							#return int_df1, flag
						if (len(checkn) >= 2):
							print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
							flag = 1
							#return int_df1, flag
						v=1
						for prot_res in prot_chain:
							for prot_atoms in prot_res:
								for rna_res in rna_chain:
									rna_resname = rna_res.resname
									for rna_atoms in rna_res:
										distance = prot_atoms-rna_atoms
										print("")
										if (distance<= 3.5) and not prot_res.resname=='HOH':
											dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
											int_df1 = int_df1.append(dict1,ignore_index=True)
											v=0
											print("")
						if v==0:
							b1 = [x.strip() in nal1 for x in int_df1['na_res']]
							temp1 = int_df1[b1]
							b2 = [x.strip() in pal1 for x in temp1['prot_res']]
							df_inter = temp1[b2]
							df_inter=pd.DataFrame(int_df1)
						else:
							df_inter=pd.DataFrame(int_df1)
						print("")
						return df_inter
					def f7_energy2(self, df_inter, aa_param,na_param):
						df1 = df_inter.copy(deep=True)			   
						vdwrad = 'Vdwradrna'
						vdweps = 'Vdwepsrna'
						total_energy = 0
						df1['prot_atmtype'] ='NA'
						df1['prot_charge'] = 0
						df1['prot_vdwradius'] =0
						df1['prot_vdweps'] =0
						df1['na_atmtype'] ='NA'
						df1['na_charge'] =0
						df1['na_vdwradius'] =0
						df1['na_vdweps'] =0
						df1['Vdw_energy'] =0
						#numcols = len(df1.columns)
						for i in df1.index:
							print("")
							idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
							# print (idx1)
							if (len(idx1) == 0 ):
								print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
							else:
								df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
								df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
								df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
								df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
							idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
							if (len(idx2) == 0 ):
								# print (list(df1.columns))
								print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
							else:
								df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
								df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
								df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
								df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
							eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
							vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
							chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
							vdwr2 = vdwr*vdwr
							vdwr6 = vdwr2*vdwr2*vdwr2
							vdwr12 = vdwr6*vdwr6
							epssqrt = math.sqrt(eps)
							aec12ab = epssqrt*vdwr12
							aec6ab = 2*epssqrt*vdwr6
							rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
							rs = float(1)/rijs
							rij2=rs
							rij6 = rij2*rij2*rij2
							rij12 = rij6*rij6
							v1 = (aec12ab*rij12)-(aec6ab*rij6)
							v2 = chrg*rij2
							v12 =v1+v2
							print (v1,v2,v12)
							total_energy += v12
							if (v12>0 and df1.loc[i,'distance']<=3):
								df1.loc[i,'vdw_energy'] = 0
								df1.loc[i,'ele_energy'] = 0
								df1.loc[i,'total_energy'] =0
							else:
								df1.loc[i,'total_energy'] = v12
								df1.loc[i,'vdw_energy'] = v1
								df1.loc[i,'ele_energy'] = v2
							print("")
						return df1
					def f8_interaction_type (self, df_inter):
						int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
						print("")
						temp_ptr = open('temp.txt','w')
						for index, dfrow in df_inter.iterrows():
							temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
							in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
							if in_t in int_dict:
								int_dict[in_t] += 1
							else:
								if not 'H' in in_t:
									int_dict.update({in_t:1})
						temp_ptr.close()
						return int_dict
					def f9_energy_div (self, energy_df):			  
						energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
						if len(energy_df) ==0:
							return energy_dict
						prot_main_atm = ['N','CA','C']
						na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
						na_sidech_b4 = np.invert(na_mainch_b1)
						prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
						prot_sidech_b3 = np.invert(prot_mainch_b2)
						energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
						energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
						energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
						energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
						energy_dict['total1'] = energy_df['total_energy'].sum()
						return energy_dict

				prot_chain=[]
				dna_chain=[]
				with open("apo_"+pdb_id_up+".pdb") as file:
					for rows in file.readlines():
						if rows[0:4]=='ATOM':
							x=rows[21].strip()
							if not x=='':
								prot_chain.append(x)
				with open("dna_"+pdb_id_up+".pdb") as file:
					for rows in file.readlines():
						if rows[0:4]=='ATOM':
							x=rows[21].strip()
							if not x=='':
								dna_chain.append(x)	
				dna_chain=list(set(dna_chain))
				prot_chain=list(set(prot_chain))
				aa_param = pd.read_csv('aa_20vdrch.csv')
				na_param = pd.read_csv('rna_4vdrch.csv')
				op_count=[]
				bind=[]
				mcmc=[]
				no_count=[]	
				tot=[]
				n_atoms=[]
				ele1=[]
				nonpolar,polar,charged=[],[],[]
				for i in prot_chain:
					for j in dna_chain:
						inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
						df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
						df_inter = pd.DataFrame(df_inter)
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						sys.stdout.flush()
						print("")
						if len(df_inter.columns)>3:
							atoms_involve = inst2.f8_interaction_type(df_inter)
							print("")
							df_at=pd.Series(atoms_involve).to_frame()
							df_at.columns=['count']
							df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
							if 'NO1' in df_at.index:
								no_count.append(df_at.loc[['NO']].values)
							tot.append(df_at['count'].sum())
							for p0 in df_inter['prot_resno'].tolist():
								bind.append(i+"_"+str(int(p0)))
							aa={'non-polar':['GLY','ALA','VAL','LEU','MET','ILE','PHE','TYR','TRP'],'polar':['SER','THR','CYS','PRO','ASN','GLN'],'charged':['LYS','ARG','HIS','ASP','GLU']} 
							df_c=df_inter[['prot_res','prot_resno']].drop_duplicates(ignore_index=True)
							df_c.drop_duplicates(ignore_index=True)
							df_c=df_c.replace(aa['non-polar'],'non-polar')
							df_c=df_c.replace(aa['polar'],'polar')
							df_c=df_c.replace(aa['charged'],'charged')
							print("")
							df_new=df_c.groupby(df_c.iloc[:,0])["prot_res"].count()
							if "charged" in df_new.index:
								charged.append(int(df_new.loc[['charged']].values))
							energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
							energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')	
							energy_dict = inst2.f9_energy_div(energy_df)
							df_en=pd.Series(energy_dict).to_frame()
							df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
							mcmc.append(df_en.loc[['mc_mc']].values)
							ele1.append(energy_df['ele_energy'].sum())
							tot.append(df_at['count'].sum())
				
				mcmc1=float(sum(mcmc))
				NO1=float(sum(no_count))/float(sum(tot))
				ele=sum(ele1)
				with open(r"apo_bind_"+pdb_id_up+".pdb","w") as file1:
					with open(r"apo_"+pdb_id_up+".pdb") as file:
						for i in file.readlines():
							if i[0:4]=="ATOM":
								print(i[21:26])
								print("\n")
								try:
									chain=i[21:26].split()[0]
									number=i[21:26].split()[1]
								except:
									chain=i[21].strip()
									number=i[22:26].strip()
								print(number)
								if chain+"_"+number in bind:
									file1.write(i)
				p4=subprocess.Popen("./hbplus "+ pdb_id_up+".pdb", shell=True)
				p4.wait()
				k=1
				ms1=0
				with open(pdb_id_up+".hb2") as file:
					for j in file.readlines():
						if "dist" in j:
							k=0
						if k==0:
							#print(i.split())
							#print(j)
							j=j.split()
							if j[5]=="MS":
								ms1=ms1+1
				surface1=1000
				os.chdir("../"+randname)
				#solv=fold_dict['Solvation Hydrophobic']
				#co1=float(sum(co_count))/float(sum(tot))
				os.system("chmod -R 777 apo_bind_"+pdb_id_up+".pdb")
				#os.system("chmod -R 777 apo_1po6.pdb")
				shutil.copy("apo_bind_"+pdb_id_up+".pdb", "../vossvolvox-master/xyzr/apo_bind_"+pdb_id_up+".pdb")
				os.chdir("../")
				os.chdir("vossvolvox-master/xyzr")
				os.system("chmod +x pdb_to_xyzr")
				p4=subprocess.Popen("./pdb_to_xyzr apo_bind_"+pdb_id_up+".pdb > apo_bind_"+pdb_id_up+".xyzr",shell=True)
				p4.wait()
				os.system("chmod -R 777 apo_bind_"+pdb_id_up+".xyzr")
				shutil.copy("apo_bind_"+pdb_id_up+".xyzr","../apo_bind_"+pdb_id_up+".xyzr")
				os.chdir("../")
				p4=subprocess.Popen("./bin -i apo_bind_"+pdb_id_up+".xyzr > volume_final3",shell=True)
				p4.wait()
				#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
				#p4.wait()
				with open("volume_final3") as file:
					voll=file.readlines()[0].split()
					vol_bs=voll[2].strip()
					surface_bs=voll[3].strip()
				#print(vol1)
				#print(surface1)
				
				surface1=1000
				os.chdir("../"+randname)
				#solv=fold_dict['Solvation Hydrophobic']
				#co1=float(sum(co_count))/float(sum(tot))
				os.system("chmod -R 777 apo_"+pdb_id_up+".pdb")
				#os.system("chmod -R 777 apo_1po6.pdb")
				shutil.copy("apo_"+pdb_id_up+".pdb", "../vossvolvox-master/xyzr/apo_"+pdb_id_up+".pdb")
				os.chdir("../")
				os.chdir("vossvolvox-master/xyzr")
				os.system("chmod +x pdb_to_xyzr")
				p4=subprocess.Popen("./pdb_to_xyzr apo_"+pdb_id_up+".pdb > apo_"+pdb_id_up+".xyzr",shell=True)
				p4.wait()
				os.system("chmod -R 777 apo_"+pdb_id_up+".xyzr")
				shutil.copy("apo_"+pdb_id_up+".xyzr","../apo_"+pdb_id_up+".xyzr")
				os.chdir("../")
				p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr > volume_final3",shell=True)
				p4.wait()
				#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
				#p4.wait()
				with open("volume_final3") as file:
					voll=file.readlines()[0].split()
					vol=voll[2].strip()
					surface1=voll[3].strip()
				#print(vol1)
				#print(surface1)
				if len(charged)>0:
					charged=sum(charged)
				else:
					charged=[[0]]
				
				predval=-1.4102707323581358 *(mcmc1)+-15.741040601006999 *(NO1)+-0.22933623750146914 *(charged)+-6.408747686469958 *(ele)+-5.820396990846385
				#-25.147730836390743 *float(ON1)+0.0010733817791259064 *float(surface1)+-0.003777867307730284 *float(surface_bs)+-0.20427012318844287 *float(ms1)+0.0005242249248441044 *float(vol_bs)+-7.897981411892616
				predval="%.2f" % predval
				if pdb_id_up=='input':
					resultout.write("User input")
				else:
					resultout.write(pdb_id_up)
				resultout.write("\n")
				#print (chain)
				prot_chain=[]
				dna_chain=[]
				os.chdir("../"+randname)
				with open("apo_"+pdb_id_up+".pdb") as file:
					for rows in file.readlines():
						if rows[0:4]=='ATOM':
							x=rows[21].strip()
							if not x=='':
								prot_chain.append(x)
				with open("dna_"+pdb_id_up+".pdb") as file:
					for rows in file.readlines():
						if rows[0:4]=='ATOM':
							x=rows[21].strip()
							if not x=='':
								dna_chain.append(x)	
				dna_chain=list(set(dna_chain))
				prot_chain=list(set(prot_chain))
				prot_chain=",".join(prot_chain)
				resultout.write(str(prot_chain))
				resultout.write("\n")
				#print (binlig5)
				binlig5=",".join(dna_chain)
				resultout.write(str(binlig5))
				resultout.write("\n")
				resultout.write(str(predval))
				disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
				resultout.write("\n")
				resultout.write(str(disass))
	
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
	if model2=='ds' and model1=='enz':
	#['CO1', 'Solvation Hydrophobic', 'vol_bs', 'mcmc', 'van']
		with open("result.txt","w") as resultout:
			import urllib.request
			
			try:
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.par",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				#roll,shift=[],[]
				Buckle=[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "parameters" in j:
							k=1
							#Buckle.append([])
							
							Buckle.append([])
							c=c+1
						if k==1:
							if not "#" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								Buckle[c].append(float(s[4]))
								#roll[c].append(float(s[11]))
								#shift[c].append(float(s[7]))
								l=0
					if l==1:
						print(i+"error:single_base_step")
					Buckle1=sum(Buckle[0])/len(Buckle[0])
					#roll1=sum(roll[0])/len(roll[0])
			except:
				try:
					l=1
					c=-1
					k=0
					#shutil.copytree("../x3dna-v2.4/bin", "x3dna-v2.4/bin")
					'''
					os.chdir("../x3dna-v2.4/bin")
					shutil.copyfile(pdb_id_up+".pdb", "../../"+pdb_id_up+".pdb")
					print ('Content-Type: text/html\r\n')
					print('/r/n')
					#p4=subprocess.Popen("python3 check_bp.py",shell=True)
					#p4.wait()
					os.system("./x3dna_setup")
					p4=subprocess.Popen("/var/www/html/bioinfo2/cgi-bin/prapred/x3dna-v2.4/bin/x3dna_setup >> pplo",shell=True)
					p4.wait()
					p4=subprocess.Popen("./find_pair 1aay.pdb",shell=True)
					p4.wait()	
					p4=subprocess.Popen("./analyze ../"+pdb_id_up+".bps",shell=True)
					p4.wait()
					'''
					os.chdir("../x3dna-v2.4/bin")
					Buckle=[]
					with open("bp_step.par") as file:	
						for j in file.readlines():
							if "parameters" in j:
								k=1
								#Buckle.append([])
								
								Buckle.append([])
								c=c+1
							if k==1:
								if not "#" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									Buckle[c].append(float(s[4]))
									#roll[c].append(float(s[11]))
									#shift[c].append(float(s[7]))
									l=0
					if l==1:
						print(i+"error:single_base_step")
					Buckle1=sum(Buckle[0])/len(Buckle[0])
					os.chdir("../../"+randname)
					print(os.getcwd())							
				except:
					print("</div>")
					print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
					print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
					print("</div>")
					exit()
	
			print("<html>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[modd]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['A','G','C','U']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						#return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						#return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									print("")
									if (distance<= 3.5):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
										print("")
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
						#df_inter=pd.DataFrame(int_df1)
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			class ProtSelect(Select):
				warnings.simplefilter('ignore', PDBConstructionWarning)
				warnings.simplefilter('ignore', FutureWarning)
				def accept_residue(self, residue):
					if not is_aa(residue, standard=True):
						res = residue.id[0]
						if not res == "W":
							return True 
					else:
						return False
			class ProtSelect1(Select):
				def accept_residue(self, residue):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					if is_aa(residue, standard=True):
						return True 
					else:
						return False
			
			parser = PDBParser()
			structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb1")
			modell = structure[modd]
			io = bpdb.PDBIO()
			io.set_structure(modell)
			io.save('dna_'+pdb_id_up+'.pdb1', ProtSelect())
			io.save('apo_'+pdb_id_up+'.pdb1', ProtSelect1())
			with open("apo_"+pdb_id_up+".pdb1") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb1") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('rna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			mcmc=[]
			bind=[]
			no_count=[]
			cc_count=[]
			tot=[]
			van1=[]
			s_atoms=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb1', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb1',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_df.to_csv(pdb_id_up+"_"+i+"_"+j+'interaction_energy.csv')
						for oo in list(set(energy_df['prot_resno'])):
							bind.append(i+"_"+str(int(oo)))				
						atoms_involve = inst2.f8_interaction_type(df_inter)
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						if 'NO' in df_at.index:
							no_count.append(int(df_at.loc[['NO']].values))
						if 'CC' in df_at.index:
							cc_count.append(int(df_at.loc[['CC']].values))
						for ao in df_at.index:
							if ao[0]=='S':
								print(ao)
								s_atoms.append(df_at.loc[[ao]].values)#+df_at.loc[['NC']].values+df_at.loc[['NN']].values+df_at.loc[['NP']].values)
						tot.append(df_at['count'].sum())
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						mcmc.append(df_en.loc[['mc_mc']].values)
						van1.append(energy_df['vdw_energy'].sum())
			s_atoms1=sum(s_atoms)	
			no_count=sum(no_count)
			cc_count=sum(cc_count)					
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			mcmc1=float(sum(mcmc))
			#os.system("chmod -R +x {}".format(foldx))
			bind=list(set(bind))
			print(bind)
			with open(r"apo_bind_"+pdb_id_up+".pdb","w") as file1:
				with open(r"apo_"+pdb_id_up+".pdb") as file:
					for i in file.readlines():
						if i[0:4]=="ATOM":
							print(i[21:26])
							print("\n")
							try:
								chain=i[21:26].split()[0]
								number=i[21:26].split()[1]
							except:
								chain=i[21].strip()
								number=i[22:26].strip()
							print(number)
							if chain+"_"+number in bind:
								file1.write(i)

			
			os.chdir("../foldx5Linux64.tar_/")
			sys.stdout.flush()
			print ('Content-Type: text/html\r\n')
			print('/r/n')
			
			p2=subprocess.Popen("./foldx_20231231 --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithRNA=True",  shell=True)
			p2.wait()
			
			fold_par=[]
			fold_val=[]
			cn_count=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					#print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
			#os.system("chmod -R +x {}".format(naccess))
			os.chdir("../"+randname)
			torsional_clash=fold_dict['torsional clash']
			#SC1=float(sum(sc_count))/float(sum(tot))
			'''
			os.system("chmod -R 777 apo_bind_"+pdb_id_up+".pdb")
			#os.system("chmod -R 777 apo_1po6.pdb")
			
			shutil.copy("apo_bind_"+pdb_id_up+".pdb", "../vossvolvox-master/xyzr/apo_bind_"+pdb_id_up+".pdb")
			os.chdir("../")
			os.chdir("vossvolvox-master/xyzr")
			os.system("chmod +x pdb_to_xyzr")
			p4=subprocess.Popen("./pdb_to_xyzr apo_bind_"+pdb_id_up+".pdb > apo_bind_"+pdb_id_up+".xyzr",shell=True)
			p4.wait()
			os.system("chmod -R 777 apo_bind_"+pdb_id_up+".xyzr")
			shutil.copy("apo_bind_"+pdb_id_up+".xyzr","../apo_bind_"+pdb_id_up+".xyzr")
			os.chdir("../")
			p4=subprocess.Popen("./bin -i apo_bind_"+pdb_id_up+".xyzr > volume_final3",shell=True)
			p4.wait()
			#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
			#p4.wait()
			with open("volume_final3") as file:
				voll=file.readlines()[0].split()
				vol1=voll[2].strip()
				surface1=voll[3].strip()
			print(vol1)
			print(surface1)
			os.chdir(path)
			#print(co1,solv,vol1,mcmc1,van_val)
			'''
			#urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
			l=1
			c=-1
			k=0
			try:
				import urllib.request
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				Incl,rise=[],[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "helical" in j:
							k=1
							Incl.append([])
							rise.append([])
							c=c+1
						if k==1:
							if "ave" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								Incl[c].append(float(s[4]))
								rise[c].append(float(s[3]))
								l=0
								break
					if l==1:
						print(i+"error:single_base_step")
					Incl1=sum(Incl[0])/len(Incl[0])
					rise1=sum(rise[0])/len(rise[0])
					print('other_class')
			except:
				try:
					import urllib.request
					urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
					l=1
					c=-1
					k=0
					Incl,rise=[],[]
					with open(pdb_id_up+"_summary.txt") as file:	
						for j in file.readlines():
							if "helical" in j:
								k=1
								Incl.append([])
								rise.append([])
								c=c+1
							if k==1:
								if "ave" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									Incl[c].append(float(s[4]))
									rise[c].append(float(s[3]))
									l=0
									break
						if l==1:
							print(i+"error:single_base_step")
						Incl1=sum(Incl[0])/len(Incl[0])
						rise1=sum(rise[0])/len(rise[0])
					print('other_class')
				except:
					try:
						#shutil.copytree("../x3dna-v2.4/bin", "x3dna-v2.4/bin")
						os.chdir("../x3dna-v2.4/bin")
						'''
						p4=subprocess.Popen("./find_pair "+pdb_id_up+".pdb "+pdb_id_up+".bps")
						p4.wait()	
						p4=subprocess.Popen("./analyze "+pdb_id_up+".bps")
						p4.wait()	
						'''
						l=1
						c=-1
						k=0
						Incl,rise=[],[]
						with open(pdb_id_up+".out") as file:	
							for j in file.readlines():
								if "helical" in j:
									k=1
									Incl.append([])
									rise.append([])
									c=c+1
								if k==1:
									if "ave" in j:
										j=j.strip()
										s=[str for str in j.split(" ") if str.strip()]
										Incl[c].append(float(s[4]))
										rise[c].append(float(s[3]))
										l=0
										break
							if l==1:
								print(i+"error:single_base_step")
							Incl1=sum(Incl[0])/len(Incl[0])
							rise1=sum(rise[0])/len(rise[0])
						
					except:
						print("</div>")
						print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
						print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
						print("</div>")
						exit()
					os.chdir("../../"+randname)
			'''
			except:
				print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
				print("Error in step parameter calcualtion. Please check if the RNA is double stranded")
				print("</div>")
				exit()	
			'''
			p3=subprocess.Popen("./naccess "+pdb_id_up+".pdb -h", shell=True)
			p3.wait()
			rsa={}
			pol={'N':'polar','O':'polar','C':'nonpolar','S':'nonpolar'}
			rsa['polar']=0
			rsa['nonpolar']=0
			with open(pdb_id_up+".asa") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						rows=[str for str in rows.split(" ") if str.strip()]
						if rows[4]+"_"+rows[5] in bind:
							if rows[2][0] in pol:	
								rsa[pol[rows[2][0]]]=rsa[pol[rows[2][0]]]+float(rows[-2])
			non_polar_asa=rsa["nonpolar"]
			#print(SC1,s_atoms1,bp_slide,bp_shift,non_polar_asa)
			
			#predval=-223.55348105782042 *(SC1)+1.9539991839704902 *(s_atoms1)+-1.308390501332273 *(bp_slide)+0.6994755185505096 *(bp_shift)+-0.0027643892605309084 *(non_polar_asa)+-9.980486180799497
			predval=0.2375130663618559 *float(Incl1)+-0.1486949148770158 *float(no_count)+0.2666932690410544 *float(torsional_clash)+-0.18556919797182764 *float(cc_count)+-0.10074606464540935 *float(Buckle1)+2.220224892638656 *float(rise1)+-17.42396441340367
			print(Incl1,no_count,torsional_clash,cc_count,Buckle1,rise1)
			#print(ppp)
			predval="%.2f" % predval
			#if pdb_id_up=='input':
			#	resultout.write("User input")
			#else:
			#	resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval))
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
			
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
##########################################################################################					Double-all-beta-nreg		  ##########################################################################################
	if model2=='ds' and model1 =='reg':
		with open("result.txt","w") as resultout:
			single_shift,single_twist=[],[]
			method_bind=1			
			try:
				import urllib.request
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				Tip=[]
				disp=[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "helical" in j:
							k=1
							disp.append([])
							Tip.append([])
							c=c+1
						if k==1:
							if "ave" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								disp[c].append(float(s[2]))
								Tip[c].append(float(s[5]))
								l=0
								break
					if l==1:
						print(i+"error:single_base_step")
					disp=disp[0][0]
					Tip=Tip[0][0]
					print('other_class')
			except:
				try:
					import urllib.request
					urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
					l=1
					c=-1
					k=0
					Tip=[]
					disp=[]
					with open(pdb_id_up+"_summary.txt") as file:	
						for j in file.readlines():
							if "helical" in j:
								k=1
								disp.append([])
								Tip.append([])
								c=c+1
							if k==1:
								if "ave" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									disp[c].append(float(s[2]))
									Tip[c].append(float(s[5]))
									l=0
									break
						if l==1:
							print(i+"error:single_base_step")
						disp=disp[0][0]
						Tip=Tip[0][0]
						print('other_class')
				except:
					try:
						os.chdir("../x3dna-v2.4/bin")
						shutil.copyfile("../../"+randname+"/"+pdb_id_up+".pdb", pdb_id_up+".pdb")
						os.system("chmod 777 "+pdb_id_up+".pdb")
						p4=subprocess.Popen("./find_pair "+pdb_id_up+".pdb "+pdb_id_up+".bps",shell=True)
						p4.wait()	
						p4=subprocess.Popen("./analyze "+pdb_id_up+".bps",shell=True)
						p4.wait()	
						l=1
						c=-1
						k=0
						Tip=[]
						disp=[]
						with open(pdb_id_up+".out") as file:	
							for j in file.readlines():
								if "helical" in j:
									k=1
									disp.append([])
									Tip.append([])
									c=c+1
								if k==1:
									if "ave" in j:
										j=j.strip()
										s=[str for str in j.split(" ") if str.strip()]
										disp[c].append(float(s[2]))
										Tip[c].append(float(s[5]))
										l=0
										break
							if l==1:
								print(i+"error:single_base_step")
							disp=disp[0][0]
							Tip=Tip[0][0]
							print('other_class')
						
					except:
						print("</div>")
						print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
						print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
						print("</div>")
						exit()
					os.chdir("../../"+randname)
			import urllib.request 
			k=0
			c=0
			l=1
			
			try:
				urllib.request.urlretrieve(r"http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.pars",pdb_id_up+"_bp_step.pars")
				with open(pdb_id_up+"_bp_step.pars") as file:
					for j in file.readlines():
						if "***local step parameters***" in j:
							k=1
						if k==1:
							if not "#" in j:
								j=j.strip()
								#print(j)
								s=[str for str in j.split(" ") if str.strip()]
								single_shift.append(float(s[1]))

								l=0
					if l==1:
						print(i+"error:single_base_step")
					single_shift1=sum(single_shift)/len(single_shift)
			except:
					try:
						#shutil.copytree("../x3dna-v2.4/bin", "x3dna-v2.4/bin")
						os.chdir("x3dna-v2.4/bin")
						p4=subprocess.Popen("./find_pair "+pdb_id_up+" "+pdb_id_up+".bps",shell=True)
						p4.wait()	
						p4=subprocess.Popen("./analyze "+pdb_id_up+".bps",shell=True)
						p4.wait()	
						l=1
						c=-1
						k=0
						with open(pdb_id_up+"_bp_step.pars") as file:
							for j in file.readlines():
								if "***local step parameters***" in j:
									k=1
								if k==1:
									if not "#" in j:
										j=j.strip()
										#print(j)
										s=[str for str in j.split(" ") if str.strip()]
										single_shift.append(float(s[1]))

										l=0
							if l==1:
								print(i+"error:single_base_step")
							single_shift1=sum(single_shift)/len(single_shift)
						
					except:
						print("</div>")
						print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
						print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
						print("</div>")
						exit()	
					os.chdir("../../"+randname)
			print("</div></center></div>")
			print("<html>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[modd]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['A','G','C','U']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						#return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						#return int_df1, flag
					v=1
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									print("")
									if (distance<= 3.5):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
										print("")
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame()
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('rna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			tot=[]
			c_atoms=[]
			mcsc=[]
			bind=[]
			cp_count=[]
			oo_count,sn_count,on_count=[],[],[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					if len(df_inter.columns)>3:	
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						print("")
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						print("")
						if 'CP' in df_at.index:
							cp_count.append(int(df_at.loc[['CP']].values))#.values/df_at['count'].sum()
						if 'SN' in df_at.index:
							sn_count.append(df_at.loc[['SN']].values)#df_at['count'].sum()
						if 'ON' in df_at.index:
							on_count.append(df_at.loc[['ON']].values)#/df_at['count'].sum()
						tot.append(df_at['count'].sum())
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						mcsc.append(df_en.loc[['mc_sc']].values)
						for oo in list(set(energy_df['prot_resno'])):
							bind.append(i+"_"+str(int(oo)))
						print("")
						for ao in df_at.index:
							if ao[0]=='C':
								print(ao)
								c_atoms.append(df_at.loc[[ao]].values)#+df_at.loc[['NC']].values+df_at.loc[['NN']].values+df_at.loc[['NP']].values)
						#c_atoms.append(df_at.loc[['CO']].values+df_at.loc[['CC']].values+df_at.loc[['CN']].values+df_at.loc[['CP']].values)
			
			mcsc1=float(sum(mcsc))	
			print(on_count)
			CP=sum(cp_count)
			#OO1=sum(oo_count)/sum(tot)
			#SN1=sum(sn_count)/sum(tot)
			#print(ON1,OO1,SN1,roll1,ss1)
			c_atoms1=sum(c_atoms)
			os.chdir("../foldx5Linux64.tar_/")
			print(os.listdir)
			p2=subprocess.Popen("./foldx_20231231 --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithRNA=True",  shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				print(f)
				with open(f) as file:
					lis=file.readlines()
					print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				print(fold_dict)
			os.chdir("../"+randname)
			Electrostatics=fold_dict['Electrostatics']
			torsional_clash=fold_dict['torsional clash']
			from Bio.PDB import PDBParser
			from Bio.PDB.DSSP import DSSP
			p = PDBParser()
			structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
			modelp = structure[modd]
			dssp = DSSP(modelp, "apo_"+pdb_id_up+".pdb")
			#print(dssp)
			from Bio.PDB.DSSP import dssp_dict_from_pdb_file
			dssp_tuple = dssp_dict_from_pdb_file("apo_"+pdb_id_up+".pdb")
			dssp_dict = dssp_tuple[0]
			print(" ")
			dssp_sec=[]
			data={'H':'helix', 'B': 'Beta','E':'Beta','G':'helix','I':'helix','T':'coil','S':'coil','-':'coil'}
			print(dssp_dict)
			bind=list(set(bind))
			print(bind)
			for ds in bind:
				dssp_sec.append(data[dssp_dict[(ds.split("_")[0], (' ',int(ds.split("_")[1]), ' '))][1]])
			print(dssp_sec)
			dssp=dssp_sec
			he=dssp.count("helix")/len(dssp)
			be=dssp.count("Beta")/len(dssp)
			co=dssp.count("coil")/len(dssp)
			beta=be*100
			
			#print(ON[0][0],c_atoms1[0][0],torsional_clash[0],single_shift1,Electrostatics)
			
			#predval=0.3126714134106675 *float(ON[0][0])+-0.10704044140481556 *float(c_atoms1[0][0])+0.28626576717168284 *float(torsional_clash)+-0.40158231902437125 *float(single_shift1)+1.541749552385311 *float(Electrostatics)+-10.602863827046878
			predval=1.1842872776497844 *(CP)+0.9512270851661783 *(disp)+-0.35182337543487824 *(single_shift1)+-0.02078859232001462 *(beta)+-0.6319821072461116 *(mcsc1)+0.2383250239569817 *(Tip)+-7.848453291949635
			print(CP,disp,single_shift1,beta,mcsc1,Tip)
			print(predval)
			#print(iiiirrri)
			predval="%.2f" % predval
			print(predval)
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval))
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		 ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-reg		  ##########################################################################################
	
	if model2=='ds' and model1 =='str':	
		with open("result.txt","w") as resultout:
			try:
				import urllib.request
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".outs",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				incl=[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "helical" in j:
							k=1
							incl.append([])
							c=c+1
						if k==1:
							if "ave" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								incl[c].append(float(s[4]))
								l=0
								break
					if l==1:
						print(i+"error:single_base_step")
					Incl=incl[0][0]
					print('other_class')
			except:
				try:
					import urllib.request
					urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+".out",pdb_id_up+"_summary.txt")
					l=1
					c=-1
					k=0
					incl=[]
					with open(pdb_id_up+"_summary.txt") as file:	
						for j in file.readlines():
							if "helical" in j:
								k=1
								incl.append([])
								c=c+1
							if k==1:
								if "ave" in j:
									j=j.strip()
									s=[str for str in j.split(" ") if str.strip()]
									incl[c].append(float(s[4]))
									l=0
									break
						if l==1:
							print(i+"error:single_base_step")
						Incl=incl[0][0]
						print('other_class')
				except:
					print("</div>")
					print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
					print("Error in step parameter calcualtion.Error in w3DNA, or Please download the source code at GitHub and use")
					print("</div>")
					exit()	
				os.chdir("../../"+randname)
			print("")
			print("</div></center></div>")
			print("<html>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[modd]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['A','G','C','U']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						#return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						#return int_df1, flag
					v=1
					print("")
					for prot_res in prot_chain:
						for prot_atoms in prot_res:
							for rna_res in rna_chain:
								rna_resname = rna_res.resname
								for rna_atoms in rna_res:
									distance = prot_atoms-rna_atoms
									print("")
									if (distance<= 3.5):
										dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
										int_df1 = int_df1.append(dict1,ignore_index=True)
										v=0
										
					if v==0:
						b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						temp1 = int_df1[b1]
						b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						df_inter = temp1[b2]
					else:
						df_inter=pd.DataFrame(int_df1)
					print(df_inter)
					return df_inter
				def f7_energy2(self, df_inter, aa_param,na_param):
					df1 = df_inter.copy(deep=True)			   
					vdwrad = 'Vdwradrna'
					vdweps = 'Vdwepsrna'
					total_energy = 0
					df1['prot_atmtype'] ='NA'
					df1['prot_charge'] = 0
					df1['prot_vdwradius'] =0
					df1['prot_vdweps'] =0
					df1['na_atmtype'] ='NA'
					df1['na_charge'] =0
					df1['na_vdwradius'] =0
					df1['na_vdweps'] =0
					df1['Vdw_energy'] =0
					#numcols = len(df1.columns)
					for i in df1.index:
						idx1 = aa_param[(str(df1['prot_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot_res'][i]).strip() == aa_param['Res_name'])].index
						# print (idx1)
						if (len(idx1) == 0 ):
							print ("Unavailable in protein parameter: Residue - "+ str(df1['prot_res'][i]).strip()+ " , Atom - "+ str(df1['prot_atm'][i]).strip())
						else:
							df1.loc[i,'prot_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
							df1.loc[i,'prot_charge'] = list(aa_param['Charge'][idx1])[0]
							df1.loc[i,'prot_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
							df1.loc[i,'prot_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
						idx2 = na_param[(str(df1['na_atm'][i]).strip() == na_param['Atm_name']) & (str(df1['na_res'][i]).strip() == na_param['Res_name'])].index
						if (len(idx2) == 0 ):
							# print (list(df1.columns))
							print ("Unavailable in RNA parameter: Residue - "+ str(df1['na_res'][i]).strip()+ " , Atom - "+ str(df1['na_atm'][i]).strip())
						else:
							df1.loc[i,'na_atmtype'] = list(na_param['Atm_type'][idx2])[0]
							df1.loc[i,'na_charge'] = list(na_param['Charge'][idx2])[0]
							df1.loc[i,'na_vdwradius'] = list(na_param[vdwrad][idx2])[0]
							df1.loc[i,'na_vdweps'] = list(na_param[vdweps][idx2])[0]
						eps = float(df1.loc[i,'prot_vdweps']) *float(df1.loc[i,'na_vdweps'])
						vdwr = float(df1.loc[i,'prot_vdwradius']) +float(df1.loc[i,'na_vdwradius'])
						chrg = float(df1.loc[i,'prot_charge']) *float(df1.loc[i,'na_charge'])
						vdwr2 = vdwr*vdwr
						vdwr6 = vdwr2*vdwr2*vdwr2
						vdwr12 = vdwr6*vdwr6
						epssqrt = math.sqrt(eps)
						aec12ab = epssqrt*vdwr12
						aec6ab = 2*epssqrt*vdwr6
						rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
						rs = float(1)/rijs
						rij2=rs
						rij6 = rij2*rij2*rij2
						rij12 = rij6*rij6
						v1 = (aec12ab*rij12)-(aec6ab*rij6)
						v2 = chrg*rij2
						v12 =v1+v2
						# print (v1,v2,v12)
						total_energy += v12
						if (v12>0 and df1.loc[i,'distance']<=3):
							df1.loc[i,'vdw_energy'] = 0
							df1.loc[i,'ele_energy'] = 0
							df1.loc[i,'total_energy'] =0
						else:
							df1.loc[i,'total_energy'] = v12
							df1.loc[i,'vdw_energy'] = v1
							df1.loc[i,'ele_energy'] = v2
						print("")
					return df1
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
				# f9_energy_div convert the energy into main chain and side chain contacts
				#Total energy calculation between main chain and side chains
				def f9_energy_div (self, energy_df):			  
					energy_dict = {'mc_mc':0,'mc_sc':0,'sc_mc':0, 'sc_sc':0,'total1':0}
					if len(energy_df) ==0:
						return energy_dict
					prot_main_atm = ['N','CA','C']
					na_mainch_b1 = np.array([x.rstrip()[-1:] == '\'' or x.rstrip() == 'P' for x in energy_df['na_atm'] ])   
					na_sidech_b4 = np.invert(na_mainch_b1)
					prot_mainch_b2 = np.array([x.rstrip() in prot_main_atm for x in energy_df['prot_atm'] ])
					prot_sidech_b3 = np.invert(prot_mainch_b2)
					energy_dict['mc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['mc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_mainch_b1)]['total_energy'].sum()
					energy_dict['sc_mc'] = energy_df[np.array(prot_mainch_b2) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['sc_sc'] = energy_df[np.array(prot_sidech_b3) & np.array(na_sidech_b4)]['total_energy'].sum()
					energy_dict['total1'] = energy_df['total_energy'].sum()
					return energy_dict
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('rna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			scsc=[]
			catom=0
			bind=[]
			tot=[]
			op_count,nc_count,cn_count,oc_count=[],[],[],[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					print("")
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						print("")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						print("")
						
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						tot.append(df_at['count'].sum())
						#OP', 'NC1', 'Incl', 'CN1', 'OC1']
						if 'OP' in df_at.index:
							op_count.append(int(df_at.loc[['OP']].values))
							
						if 'NC' in df_at.index:
							nc_count.append(int(df_at.loc[['NC']].values))#/df_at['count'].sum())
							
						if 'CN' in df_at.index:
							cn_count.append(int(df_at.loc[['CN']].values))#/df_at['count'].sum())
						if 'OC' in df_at.index:
							oc_count.append(int(df_at.loc[['OC']].values))#/df_at['count'].sum())
						print("")
						energy_df = inst2.f7_energy2(df_inter,aa_param,na_param)
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						scsc.append(df_en.loc[['mc_mc']].values)
						#print(int(df_at.loc[['CN']].values))
			#print(res_bind_count)
			#print(sn_count)
			#print(mcmc)
			scsc1=float(sum(scsc))
			op=sum(op_count)
			NC1=sum(nc_count)/sum(tot)
			cn1=sum(cn_count)/sum(tot)
			oc1=sum(oc_count)/sum(tot)
			os.chdir("../foldx5Linux64.tar_/")
			sys.stdout.flush()
			print ('Content-Type: text/html\r\n')
			print('/r/n')
			
			p2=subprocess.Popen("./foldx_20231231 --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithRNA=True",  shell=True)
			p2.wait()
			
			fold_par=[]
			fold_val=[]
			cn_count=[]
			for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
				fold_dict={}
				with open(f) as file:
					lis=file.readlines()
					#print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
			#os.system("chmod -R +x {}".format(naccess))
			os.chdir("../"+randname)
			Sidechain_Hbond=fold_dict['Sidechain Hbond']
			helix_dipole=fold_dict['helix dipole']
			surface1=1000
			os.chdir("../"+randname)
			#solv=fold_dict['Solvation Hydrophobic']
			#co1=float(sum(co_count))/float(sum(tot))
			#os.system("chmod -R 777 apo_bind_"+pdb_id_up+".pdb")
			#os.system("chmod -R 777 apo_1po6.pdb")
			shutil.copy(pdb_id_up+".pdb", "../vossvolvox-master/xyzr/"+pdb_id_up+".pdb")
			os.chdir("../")
			os.chdir("vossvolvox-master/xyzr")
			os.system("chmod +x pdb_to_xyzr")
			p4=subprocess.Popen("./pdb_to_xyzr "+pdb_id_up+".pdb > "+pdb_id_up+".xyzr",shell=True)
			p4.wait()
			os.system("chmod -R 777 "+pdb_id_up+".xyzr")
			shutil.copy(pdb_id_up+".xyzr","../"+pdb_id_up+".xyzr")
			os.chdir("../")
			p4=subprocess.Popen("./bin -i "+pdb_id_up+".xyzr > volume_final3",shell=True)
			p4.wait()
			#p4=subprocess.Popen("./bin -i apo_"+pdb_id_up+".xyzr> vol_out.txt",shell=True)
			#p4.wait()
			with open("volume_final3") as file:
				voll=file.readlines()[0].split()
				vol=voll[2].strip()
				surface=voll[3].strip()
			os.chdir("../"+randname)
			#predval=0.14505778669322977 *float(op)+32.54575376213373 *float(nc1)+0.15474311201557356 *float(Incl[0])+-19.482560241915714 *float(cn1)+10.53966507800749 *float(oc1)+-12.140663440179148
			predval=-0.7302651215184853 *float(scsc1)+-1.7595942771081356e-05 *float(vol)+4.151051035573872 *float(helix_dipole)+42.57912804864454 *float(NC1)+0.24237181545309472 *float(Sidechain_Hbond)+0.16806118333997772 *float(Incl)+-12.163579352795393

			#print(op,nc1,Incl[0],cn1,oc1)
			print(predval)
			#print(pppp)
			#print(ppppp)
			predval="%.2f" % predval
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			#print (chain)
			prot_chain=",".join(prot_chain)
			resultout.write(str(prot_chain))
			resultout.write("\n")
			#print (binlig5)
			binlig5=",".join(dna_chain)
			resultout.write(str(binlig5))
			resultout.write("\n")
			resultout.write(str(predval))
			disass= '{:.3g}'.format(math.exp(float(predval)/(0.0019*298.15)))
			resultout.write("\n")
			resultout.write(str(disass))
	
	
##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		 ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
##########################################################################################					Double-alpha-beta-nreg		  ##########################################################################################
	
########## common for all the classifications - start


timetaken=timeit.default_timer() - start

with open("result.py", "w") as polyout:
	print("ii")
	polyout.write("""#!/opt/websites/anaconda/envs/hari_37/bin/python\nimport cgi\nimport cgitb; cgitb.enable()\nprint ('Content-Type: text/html\\r\\n')\n""")
	g=open("index.txt").readlines()
	for gg in g:
		gg=gg.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(gg))
		polyout.write("\n")
	f=open("result.txt").readlines()
	# polyout.write ("print ('<center>Time taken for the calculation: {:.2f} seconds</center>')".format(timetaken))
	# polyout.write("\n")
	# polyout.write("print ('<br/>')")
	# polyout.write("\n")
	# polyout.write ("print ('<center>Model selected: {}</center>')".format(modeltype))
	# polyout.write("\n")
	# polyout.write("print ('<br/>')")
	# polyout.write("\n")


	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>Time taken</b></td><td>{:.2f} seconds</td></td>')".format(timetaken))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	if model2!="ss":
		if model=="alpha" or model=="beta":
			polyout.write("print ('<td><b>Model selected</b></td><td>"+modeltype2+"-"+modeltype+"</td>')")
		else:
			polyout.write("print ('<td><b>Model selected</b></td><td>"+modeltype2+"-"+modeltype+"-"+modeltype1+"</td>')")
	if model2=="ss":
		polyout.write("print ('<td><b>Model selected</b></td><td>"+modeltype2+"</td>')")
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")



	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>PDB ID</b></td><td><b>protein Chain(s)</b></td><td><b>RNA chain (s) </b></td><td><b>Predicted &Delta;G (kcal/mol)</b></td><td><b>K<sub>d</sub> (M)</b></td>')")
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	print(f)
	polyout.write("print ('<td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>')".format(f[0].rstrip(),f[1].rstrip(),f[2].rstrip(),f[3].rstrip(),f[4].rstrip()))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	h=open("footer.txt").readlines()
	for hh in h:
		hh=hh.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(hh))
		polyout.write("\n")
os.system("chmod +x result.py")
os.system("ls > remfile")
'''
redirectURL = "/cgi-bin/%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print ('	<title>You are going to be redirected</title>')
'''
redirectURL = "%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print('	<title>You are going to be redirected</title>')
'''
remfile=open("remfile").readlines()
for rf in remfile:
	rf=rf.rstrip()
	if rf =='molecules':
		os.rmdir('{}'.format(rf))
	elif rf =='autodock':
		shutil.rmtree('autodock')
	elif rf not in ['result.txt','result.py',"style4.css"]:
		os.remove('{}'.format(rf))
'''
#print ('	Redirecting... <a href="%s">Click here if you are not redirected</a>' % redirectURL)
print ('</body>')
print ('</html>')

########## common for all the classifications - end
