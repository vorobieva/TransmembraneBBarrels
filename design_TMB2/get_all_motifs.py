import glob
import itertools
import shutil
import re
import os

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

for pdb_file in glob.glob("*.pdb"):
	pdb_name = pdb_file.split(".")[0]
	nets = {}
	Y68 = []
	Y10 = []
	num_nets = len(glob.glob(pdb_name+"/*.cst"))
	for i in range (1,num_nets+1):
		net = pdb_name + "/" + pdb_name + "_0001_network_" + str(i)
		with open(net+".cst", 'r') as cst_file:
			net_name = ""
			cst = ""
			for line in cst_file:
				if " network_" in line:
					vals = line.split()
					net_name = vals[2]
					if "Y_10" in net_name:
						Y10.append((net_name, i))
					elif "Y_68" in net_name:
						Y68.append((net_name, i))
				elif "#AtomPair" in line:
					cst = line[1:]
			nets[i] = (net_name, cst)
	prod = list(itertools.product(Y68, Y10))
	for p in prod:
		n_comb = ",".join(n[0] for n in p)
		comb = [x[1] for x in p]
		c_comb = "".join(nets[y][1] for y in comb)
		nets[max(nets, key=int) + 1] = (n_comb, c_comb)

	for key in nets:
		source_pdb = pdb_name + "/" + pdb_name + "_0001_network_" + str(key) + ".pdb"
		new_path = "../round2/"+pdb_name+"_"+str(key)+"/"
		try:
			os.mkdir(new_path)
		except OSError:
			print ("Creation of the directory %s failed" % new_path)	
	
		muts = nets[key][0].split(",")
		pairs = list(chunks(muts,2))
		m ={}
		for j in range(0, len(pairs)):
			m_name = ",".join(pairs[j])
			for k, d in nets.items():
				if nets[k][0] == m_name:
					s_pdb = pdb_name + "/" + pdb_name + "_0001_network_" + str(k) + ".pdb"
					res_id = re.findall("\d+", m_name)
					for res in res_id:
						coord = ""
						with open(s_pdb, 'r') as template:
							for line in template:
								if "ATOM" in line:
									if line.split()[5] == res:
										coord += line
						m[res] = coord
		
		with open(pdb_file, 'r') as template:
			with open(new_path+pdb_name+"_"+str(key) + ".pdb", 'w') as out_f:
				for l in template:
					if "ATOM" in l and l.split()[5] in m:
						if l.split()[2] == "N":
							out_f.write(m[l.split()[5]])
						else:
							pass
					elif "ATOM" in l:
						out_f.write(l)
				for r in m:
					out_f.write("REMARK PDBinfo-LABEL:   %s HBNet\n" %(r))

		with open(new_path+"cst", "w") as cst_file:
			for pair in pairs:
				p = sorted(pair)
				cst_file.write("Dihedral C %s CA %s CB %s CG %s CIRCULARHARMONIC -1.22 0.30\n" %(p[1].split("_")[2],p[1].split("_")[2],p[1].split("_")[2],p[1].split("_")[2]))
			cst_file.write(nets[key][1])
		
		with open(new_path+"core.resfile", "w") as resfile:
			resfile.write("ALLAA\n\nstart\n\n2 A PIKAA NTEDPGRKQS\n6 A PIKAA STDN\n\n16 A PIKAA DENST\n17 A PIKAA AES\n18 A PIKAA D\n19 A PIKAA G\n\n32 A PIKAA STD\n33 A PIKAA P\n34 A PIKAA DEHTY\n\n44 A PIKAA DENST\n45 A PIKAA AES\n46 A PIKAA D\n47 A PIKAA G\n\n62 A PIKAA STD\n63 A PIKAA P\n64 A PIKAA DEHTY\n\n76 A PIKAA DENST\n77 A PIKAA AES\n78 A PIKAA D\n79 A PIKAA G\n\n81 A POLAR\n\n94 A PIKAA STD\n95 A PIKAA P\n96 A PIKAA DEHTY\n\n106 A PIKAA DENST\n107 A PIKAA AES\n108 A PIKAA D\n109 A PIKAA G\n\n")
			for pair in pairs:
				resfile.write("%s A PIKAA %s\n%s A PIKAA %s\n" %(pair[0].split("_")[2],pair[0].split("_")[1],pair[1].split("_")[2],pair[1].split("_")[1]))	



