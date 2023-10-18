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
	'''
	if rna_strand=='ds':
		model=input("Please enter the structural classification of protein (all-alpha/all-beta/alpha-beta)")
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
		#print(model) #model1='reg'
		#print(model1)
	if rna_strand=='ss':
		model1=""
		model=""
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
		modeltype1="Regulatory"
	if model1=='enz':
		modeltype1="Enzyme"
	elif model1=='str':
		modeltype1="Others"
	if model2=='ds':
		modeltype2="Double strand"
	elif model2=='ss':
		modeltype2="Single strand"
		modeltype1=''
		modeltype=''
	if model=='alpha':
		modeltype="All Alpha protein"
	elif model=='beta':
		modeltype="All Beta protein"
	elif model=='alphabeta':
		modeltype="Alpha Beta protein"
	#print(model+model2+model1)
	
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
		os.system("wget 'https://files.rcsb.org/download/{}.pdb1'".format(pdb_id_up))
		#filename=wget.download('https://files.rcsb.org/download/{}.pdb1'.format(pdb_id_up))
		#filename1=wget.download('https://files.rcsb.org/download/{}.pdb'.format(pdb_id_up))
		shutil.copyfile(pdb_id_up+".pdb1", pdb_id_up+".pdb")
	elif method==2:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname="pr_res_"+uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.chdir(path)
		os.system("scp ../"+pdb_file+" input.pdb")
		pdb_id_up=input
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
	io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
	io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
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
	#shutil.copyfile("../bin","bin")
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
	if model2=='ds' and model=='none':
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
				model = structure[0]
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
									#print(int_df1)
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

		if len(dna_chain)==0:
			print("<center><div style='width:1000px;height:100px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
			print("Error in RNA chain, Kindly check the input!!!")			
			
			print("</div>")
			exit()
		if len(prot_chain)==0:
			print("<center><div style='width:1000px;height:100px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
			print("No protein chain found, Kindly check the input!!!")

			print("</div>")
			exit()		
		aa_param = pd.read_csv('aa_20vdrch.csv')
		na_param = pd.read_csv('rna_4vdrch.csv')
		op_count=[]
		bind=[]
		nn_count=[]	
		tot=[]
		nonpolar,polar,charged=[],[],[]
		for i in prot_chain:
			for j in dna_chain:
				inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
				df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
				df_inter = pd.DataFrame(df_inter)
				df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
				if len(df_inter.columns)>3:
					for p0 in df_inter['prot_resno'].tolist():
						bind.append(i+"_"+str(int(p0)))
					sys.stdout.flush()
					print("")
					
		if len(bind)==0:
			print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
			print("Binding site residues are not found. Kindly check the protein/RNA chains!!!")
			exit()
			print("</div>")
		print("</div></center></div>")
		print("<html>")
		print("<div style='width:100px;height:1px;overflow:hidden;'>")
	
		from Bio.PDB import PDBParser
		from Bio.PDB.DSSP import DSSP
		p = PDBParser()
		structure = p.get_structure("apo_"+pdb_id_up+".pdb", "apo_"+pdb_id_up+".pdb")
		modelp = structure[0]
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
		print(he,be,co)
		#print(ghg)
		if (he>=0.40)&(be<=0.05):
			model='alpha'
			modeltype="All Alpha protein"
		elif (be>=0.40)&(he<=0.05):
			model='beta'
			modeltype="All Beta protein"
		else:
			model='alphabeta'
			modeltype="Alpha Beta protein"
	if model2=='ss':	
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
					model = structure[0]
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
			nn_count=[]	
			tot=[]
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
						if 'OP' in df_at.index:
							op_count.append(df_at.loc[['OP']].values)
						tot.append(df_at['count'].sum())	
						if 'NN' in df_at.index:
							nn_count.append(df_at.loc[['NN']].values)
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
			
			fold_energy_ion=fold_dict['Interaction Energy']
			back_bone_clash=fold_dict['backbone clash']
			os.chdir("../"+randname)
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
			
			os.chdir("../")
			os.chdir(randname)
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
			if len(nonpolar)>0:
				nonpolar=nonpolar[0]
			#print(nonpolar)
			NN=[[0]]
			if len(nn_count)>0:
				NN=nn_count[0]
			#print(NN[0][0])
			OP1=sum(op_count)/sum(tot)
			protein_residue_interface=nonpolar+polar+charged
			#print(protein_residue_interface[0])
			#ss=0
			print(NN[0][0],nonpolar,protein_residue_interface[0],back_bone_clash,OP1,surface1,ss1,fold_energy_ion)
			predval=-0.2919358520579966*float(NN[0][0])+0.3977046308700206 *float(nonpolar)+-0.24538896600374907 *float(protein_residue_interface[0])+-0.19355888842504068 *float(back_bone_clash)+82.21219881719684 *float(OP1)+0.0011023851189393974 *float(surface1)+-0.02382303236020783 *float(ss1)+0.04915970407169483 *float(fold_energy_ion)+-7.4948298558072475
		
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
	#print(model2,model,model1)
	#['IntraclashesGroup1', 'Sidechain Hbond', 'Buckle1', 'bp_tilt', 'bp_twist']
	
	if model2=='ds' and model=='alpha':
		with open("result.txt","w") as resultout:
				try:
					import urllib.request
					urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.par",pdb_id_up+"_summary.txt")
					l=1
					c=-1
					k=0
					Buckle,tilt,twist=[],[],[]
					with open(pdb_id_up+"_summary.txt") as file:	
						for j in file.readlines():
							if "parameters" in j:
								k=1
								Buckle.append([])
								tilt.append([])
								twist.append([])
								c=c+1
							if k==1:
								if not "#" in j:
									j=j.strip()
									
									s=[str for str in j.split(" ") if str.strip()]
									Buckle[c].append(float(s[4]))
									tilt[c].append(float(s[10]))
									twist[c].append(float(s[12]))
									l=0
						if l==1:
							print(i+"error:single_base_step")
						Buckle1=sum(Buckle[0])/len(Buckle[0])
						tilt1=sum(tilt[0])/len(tilt[0])
						twist1=sum(twist[0])/len(twist[0])
				except:
					print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
					print("Error in 3DNA. Please check if the RNA is double stranded")
					print("</div>")
					exit()
				print("</div></center></div>")
				print("<html>")
				print("<div style='width:100px;height:1px;overflow:hidden;'>")
				os.chdir("../foldx5Linux64.tar_/")
				sys.stdout.flush()
				print ('Content-Type: text/html\r\n')
				print('/r/n')
				
				p2=subprocess.Popen("./foldx_20231231 --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithRNA=True",  shell=True)
				p2.wait()
				
				fold_par=[]
				fold_val=[]
				fold_dict={}
				fold_dict={}
				for f in glob.glob('Interaction_'+pdb_id_up+'*.fxout'):
					#print(f)
					with open(f) as file:
						lis=file.readlines()
						#print(lis)
						par=lis[-2].split("\t")
						val=lis[-1].split("\t")
						for fo in range(len(val)):
								fold_dict[par[fo]]=val[fo]
					#print(fold_dict)
				
				IntraclashesGroup1=fold_dict['IntraclashesGroup1']
				Sidechain=fold_dict['Sidechain Hbond']
				predval=0.010211045826825409 *float(IntraclashesGroup1)+1.1910719468457243 *float(Sidechain)+-0.25327388624225355 *float(Buckle1)+-0.44561152069592536 *float(tilt1)+0.22704268462013433 *float(twist1)+-15.120669743171025
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
			#****volume_bs and surface_bs*********
#******************************************************************************bind greater than 0.29*****************************************************************	
		
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
##########################################################################################					Double-all-alpha-nreg		  ##########################################################################################
	if model2=='ds' and model=='beta':
		with open("result.txt","w") as resultout:
			import urllib.request
			try:
				urllib.request.urlretrieve("http://web.x3dna.org/data/ndb/_"+pdb_id_up.upper()+"/"+pdb_id_up.upper()+"_bp_step.par",pdb_id_up+"_summary.txt")
				l=1
				c=-1
				k=0
				roll,shift=[],[]
				with open(pdb_id_up+"_summary.txt") as file:	
					for j in file.readlines():
						if "parameters" in j:
							k=1
							#Buckle.append([])
							roll.append([])
							shift.append([])
							c=c+1
						if k==1:
							if not "#" in j:
								j=j.strip()
								s=[str for str in j.split(" ") if str.strip()]
								#Buckle[c].append(float(s[4]))
								roll[c].append(float(s[11]))
								shift[c].append(float(s[7]))
								l=0
					if l==1:
						print(i+"error:single_base_step")
					shift1=sum(shift[0])/len(shift[0])
					roll1=sum(roll[0])/len(roll[0])
			except:
				print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
				print("Error in step parameter calcualtion. Please check if the RNA is double stranded")
				print("</div>")
				exit()	
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
					model = structure[0]
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
			cp_count=[]
			mcmc=[]
			charged,nonpolar,polar=[],[],[]
			bind=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					if not pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv" in os.listdir():
						df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
						df_inter = pd.DataFrame(df_inter)
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
					if len(df_inter.columns)>3:
						aa={'non-polar':['GLY','ALA','VAL','LEU','MET','ILE','PHE','TYR','TRP'],'polar':['SER','THR','CYS','PRO','ASN','GLN'],'charged':['LYS','ARG','HIS','ASP','GLU']}
						df_c=df_inter[['prot_res','prot_resno']].drop_duplicates(ignore_index=True)
						df_c.drop_duplicates(ignore_index=True)
						df_c=df_c.replace(aa['non-polar'],'non-polar')
						df_c=df_c.replace(aa['polar'],'polar')
						df_c=df_c.replace(aa['charged'],'charged')
						print(df_c)
						df_new=df_c.groupby(df_c.iloc[:,0])["prot_res"].count()
						print(df_new)
						x=df_new.iloc[0:3].to_list()
						print(x)				
						for oo in list(set(df_inter['prot_resno'])):
							bind.append(i+"_"+str(oo))
			#os.chdir("../")
			# apo_rsa_dna
			sys.stdout.flush()
			print ('Content-Type: text/html\r\n')
			print('/r/n')
			
			p3=subprocess.Popen("./naccess apo_"+pdb_id_up+".pdb -h", shell=True)
			p3.wait()
			rsa={}
			i="apo_{}.pdb".format(pdb_id_up)
			rsa[i.split(".")[0].split("_")[1]]={}
			
			with open("apo_"+pdb_id_up+".rsa") as file:
				for rows in file.readlines():
					#print(rows)
					if rows[0:3]=='RES':
						if not rows.split(" ")[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows.split(" ")[2]]={}
						rows=[str for str in rows.split(" ") if str.strip()]
						#print(rows)
						if not rows[2].isalpha():
							new_rows=rows[0:2]
							new_rows.append(rows[2][0])
							new_rows.append(rows[2][1:])
							for kkkk in rows[3:]:
								new_rows.append(kkkk)
							rows=new_rows
						#print(rows)						
						if not rows[2] in rsa[i.split(".")[0].split("_")[1]]:
							rsa[i.split(".")[0].split("_")[1]][rows[2]]={}
						rsa[i.split(".")[0].split("_")[1]][rows[2]][rows[3]]=float(rows[4])
						#print(rsa)		
			res={}
			print(rsa)
			new_rsa_1={}
			for i in glob.glob("*atom_interaction.csv"):
				x=i.split(".")[0]
				print(x)
				if x.split("_")[0] in rsa:
					print(x)
					if len(rsa[x.split("_")[0]])>0:
						print("i")
						if not i.split(".")[0].split("_")[1] in res:
							res[i.split(".")[0].split("_")[1]]=[]
						if not i.split(".")[0].split("_")[2] in res:
							res[i.split(".")[0].split("_")[2]]=[]
						df=pd.read_csv(i)
						#print(res)
						#nar=df['prot_resno'].tolist()
						if 'prot_resno' in df.columns:
						 	for kk in range(len(df['prot_resno'].tolist())):
						 		if not df['prot_resno'][kk] in res[i.split(".")[0].split("_")[1]]:
						 			if 'C' in df['prot_atm'][kk] or 'S' in df['prot_atm'][kk]:
						 				res[i.split(".")[0].split("_")[1]].append(df['prot_resno'].tolist()[kk])
						#print(res)
						if not x in new_rsa_1:
						 	new_rsa_1[x]=[]
						for j in res[i.split(".")[0].split("_")[1]]:
							print(j)
							if str(int(j)) in rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[1][0]]:
								new_rsa_1[x].append(rsa[i.split(".")[0].split("_")[0]][i.split(".")[0].split("_")[1]][str(int(j))])
							else:
								new_rsa_1[x].append(0)
						print(new_rsa_1)
			apo_protein=sum(new_rsa_1[x]*4)

			
			p4=subprocess.Popen("./hbplus "+ pdb_id_up+".pdb", shell=True)
			p4.wait()
			mm1,ms1,sm1,ss1=0,0,0,0
			with open(pdb_id_up+".hb2") as file:
				for j in file.readlines():
					if "dist" in j:
						k=0
					if k==0:
						print(i.split())
						print(j)
						j=j.split()
						if j[5]=="SS":
							ss1=ss1+1
			print(apo_protein, shift1, roll1, ss1)
			predval=-0.0008624283376737335 *(apo_protein)+0.8680412808210729 *(shift1)+-0.25811320699421725 *(roll1)+0.09527901445842585 *(ss1)+-4.455640666734628
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
			print(disass)
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
##########################################################################################					Double-all-beta-reg		  ##########################################################################################
	if model2=='ds' and model=='alphabeta' and model1=='enz':
	#['CO1', 'Solvation Hydrophobic', 'vol_bs', 'mcmc', 'van']
		with open("result.txt","w") as resultout:
			import urllib.request 
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
					model = structure[0]
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
			modell = structure[0]
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
			co_count=[]
			tot=[]
			van1=[]
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
						if 'CO' in df_at.index:
							co_count.append(df_at.loc[['CO']].values)
						tot.append(df_at['count'].sum())
						energy_dict = inst2.f9_energy_div(energy_df)
						df_en=pd.Series(energy_dict).to_frame()
						df_en.to_csv(pdb_id_up+"_"+i+"_"+j+"_energy_side.csv")
						mcmc.append(df_en.loc[['mc_mc']].values)
						van1.append(energy_df['vdw_energy'].sum())
			van_val=sum(van1)									
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
							chain=i[21:26].split()[0]
							number=i[21:26].split()[1]
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
			solv=fold_dict['Solvation Hydrophobic']
			co1=float(sum(co_count))/float(sum(tot))
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
			print(co1,solv,vol1,mcmc1,van_val)
			predval=8.391379760832502 *float(co1)+-0.18110998367933648 *float(solv)+-0.0004974099679880273 *float(vol1)+5.430960389988162 *float(mcmc1)+-0.3920346717189805 *float(van_val)+-12.347546956921164
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
	if model2=='ds' and model=='alphabeta' and model1 =='reg':
		with open("result.txt","w") as resultout:
			single_roll,single_twist=[],[]
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
								single_roll.append(float(s[5]))
								
								l=0
					if l==1:
						print(i+"error:single_base_step")
					roll1=sum(single_roll)/len(single_roll)
			except:
				print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
				print("Error in step parameter calcualtion. Please check if the RNA is double stranded")
				print("</div>")
				exit()	
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
					model = structure[0]
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
			scmc=[]
			bind=[]
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
						if 'OO' in df_at.index:
							oo_count.append(df_at.loc[['OO']].values)#.values/df_at['count'].sum()
						if 'SN' in df_at.index:
							sn_count.append(df_at.loc[['SN']].values)#df_at['count'].sum()
						if 'ON' in df_at.index:
							on_count.append(df_at.loc[['ON']].values)#/df_at['count'].sum()
						tot.append(df_at['count'].sum())
						print("")
			
			print(res_bind_count)
			sys.stdout.flush()
			print('Content-Type: text/html\r\n')
			print('/r/n')
			k=-1
			res_bind_uni=list(set(res_bind_count))
			p4=subprocess.Popen("./hbplus "+ pdb_id_up+".pdb", shell=True)
			p4.wait()
			mm1,ms1,sm1,ss1=0,0,0,0
			with open(pdb_id_up+".hb2") as file:
				for j in file.readlines():
					if "dist" in j:
						k=0
					if k==0:
						print(i.split())
						print(j)
						j=j.split()
						if j[5]=="SS":
							ss1=ss1+1
						
			

			ON1=sum(on_count)/sum(tot)
			OO1=sum(oo_count)/sum(tot)
			SN1=sum(sn_count)/sum(tot)
			print(ON1,OO1,SN1,roll1,ss1)
			predval=27.94911646797724 *(ON1)+27.486006745904962 *(OO1)+176.9092612790746 *(SN1)+-0.05924294775064354 *(ss1)+0.055034930067101584 *(roll1)+-13.720978319340816
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
	
	if model2=='ds' and model=='alphabeta' and model1 =='str':	
		with open("result.txt","w") as resultout:
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
					if l==1:
						print(i+"error:single_base_step")
					Incl=sum(incl[0])/len(incl[0])
					print('other_class')
			except:
				print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
				print("Error in step parameter calcualtion. Please check if the RNA is double stranded")
				print("</div>")
				exit()	
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
					model = structure[0]
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
			mcmc=[]
			catom=0
			bind=[]
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
						#OP', 'NC1', 'Incl', 'CN1', 'OC1']
						if 'OP' in df_at.index:
							op_count.append(df_at.loc[['OP']].values)
						if 'NC' in df_at.index:
							nc_count.append(df_at.loc[['NC']].values/df_at['count'].sum())
						if 'CN' in df_at.index:
							cn_count.append(df_at.loc[['CN']].values/df_at['count'].sum())
						if 'OC' in df_at.index:
							oc_count.append(df_at.loc[['OC']].values/df_at['count'].sum())
						print("")
			#print(res_bind_count)
			#print(sn_count)
			#print(mcmc)
			
			op1=sum(op_count)
			nc1=sum(nc_count)
			cn1=sum(cn_count)
			oc1=sum(oc_count)
			predval=0.16870498399171072 *(op1)+49.18021763808354 *(nc1)+0.27721453797715156 *(Incl)+-22.18970702005487 *(cn1)+13.547932133596287 *(oc1)+-14.750443984324237
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
