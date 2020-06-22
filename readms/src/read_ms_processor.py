#!/usr/bin/python
import sys
from mfe_adv import read_headlst, read_file_respect_comments
import os
import struct
import re

#store each microstate information
class MSRECORD: 
   def __init__(self):
      self.counter = []
      self.stdev=0.0
      self.occ=[]


class read_ms_processor:
	def __init__(self):
		self.conformers = []
		self.ms_states = []
		self.enum_flag = 0

	def read_re_ms(self, file_path='re_ms.dat'):
		'''
		read re_ms.dat file, a txt file of microstate file, into ms_states
		input: re_ms.dat, format version: (/home/cai/mcce3.5_enum_ms)
		output: ms_state, a list
		'''

		#disregard lines with '#' in the input file
		lines = read_file_respect_comments(file_path)
	
		#read first three lines in the input file: residue number, residue names, method for microstate.
		n_spe = int(lines.pop(0).split()[-1])
		spe_list = lines.pop(0).split()
		method = lines.pop(0).split()[-1]
	
		if len(self.ms_states) > 0:
			print("WARNING: adding to non empty microstate list.")

		#read head3.lst, to read information for each conformer
		try: #succeed to read conformer information from head3.lst
			self.conformers = read_headlst()
			crg_flag=1  
		except:  # failed
			crg_flag=0

		#for monte carlo sampling, use self.counter to store microstate counts,
		#for enumerate, use occupancy to store microstate occupancy.
		if (method == "ENUMERATE"): 
			self.enum_flag =1
			print("re_ms.dat is obtained from enumeration MCCE.")
		else: self.enum_flag =0
	
		#every two lines are the information for one microstate
		#line1: conformer id of the microstate
		#line2: energy and occupancy information of the microstate
		count=0   #index of microstate
		for i in range(0, len(lines), 2):
			ms_state =MSRECORD()
			confid_list = lines[i].split()
			energies = lines[i+1].split()
			#print energies
			#print confid_list
			ms_state.id = count
			ms_state.confseq = confid_list
			ms_state.cumE = float(energies[2])
			ms_state.E = float(energies[5])
			ms_state.counter = float(energies[-1])
			ms_state.occ= float(energies[-1])
	
			#assign charge state for the microstate
			if crg_flag==0:
				ms_state.crg=None
				ms_state.confidseq=None
				ms_state.crgseq=None
			else:
				ms_state.crg = 0   #sotre total charge of the microstate
				ms_state.confidseq=[]  #store conformer id list of the microstate
				ms_state.crgseq=""     #store residue charge (string) of the microstate
				for index in ms_state.confseq:
					ms_state.crg += self.conformers[int(index)].crg
					ms_state.crgseq += str(int(self.conformers[int(index)].crg))
					confid = self.conformers[int(index)].id
					ms_state.confidseq.append(confid)
			count += 1
			self.ms_states.append(ms_state)
		return self.ms_states

	def read_ms(self, file_path='ms.dat'):
		'''
		read ms.dat file, a binary file of microstate file, into ms_states
		input: ms.dat, format version: (/home/cai/mcce3.5_enum_ms)
		output: ms_state, a list
		'''

		with open(file_path, mode="rb") as f:
			#read head3.lst, to read information for each conformer
			try: #succeed to read conformer information from head3.lst
				self.conformers = read_headlst()
				crg_flag=1
			except: #failed
				crg_flag=0


			#read number of residue, residue names, method
			byte = f.read(4)
			n_spe = struct.unpack('1i', byte)[0]
			print("There are %i residues." % (n_spe))
			spe_list=""
			for i in range(n_spe):
				byte = f.read(8)
				res_spe= struct.unpack('8s', byte)[0]
				spe_list = str(spe_list) + str(res_spe) + ","
			print("Microstates on "+str(spe_list))
			byte = f.read(9)
			method = struct.unpack('9s', byte)[0]
			print("Microstate is obtained from %s." %(method))
			if (method == b"ENUMERATE"):
				self.enum_flag=1
			else: 
				self.enum_flag=0
		
			#read information for each microstate
			count=0
			byte = f.read(2)
			while byte:
				ms_state =MSRECORD()
				confid_list=[]

				#read first conformer id of the microstate
				ms_state.id = count
				res_confid = struct.unpack('1H', byte)[0]
				confid_list.append(res_confid)
			
				#read the rest of residue conformer id of the microstate 
				for i in range(n_spe-1):  
					byte=f.read(2)
					res_confid = struct.unpack('1H', byte)[0]
					confid_list.append(res_confid)
				ms_state.confseq = confid_list
				byte = f.read(8)
				ms_state.cumE = struct.unpack('1d', byte)[0]
				byte = f.read(8)
				ms_state.E = struct.unpack('1d', byte)[0]
				if (self.enum_flag):
					byte=f.read(8)
					ms_state.occ= struct.unpack('1d', byte)[0]
				else:
					byte=f.read(4)
					ms_state.counter=struct.unpack('1i', byte)[0]
			
				#assign charge state for the microstate
				if crg_flag==0:
					ms_state.crg=None
					ms_state.confidseq=None
					ms_state.crgseq=None
				else:
					ms_state.crg = 0   #store total charge of the microstate
					ms_state.confidseq=[] #store conformer id list of the microstate
					ms_state.crgseq=""    #store residue charge (string) of the microstate
					for index in ms_state.confseq:
						ms_state.crg += self.conformers[int(index)].crg  
						ms_state.crgseq += str(int(self.conformers[int(index)].crg))
						confid = self.conformers[int(index)].id
						ms_state.confidseq.append(confid)
				count += 1
				self.ms_states.append(ms_state)
			
				#read next microstate first conformer id
				byte=f.read(2)
		return self.ms_states

	def read_old_ms(self, file_path='ms.dat', *args, **kwargs):
		'''
		read ms.dat file, a binary file of microstate file, into ms_states
		input: ms.dat, format version: (/home/cai/mcce3.5)
		output: ms_state, a list
				lowest_state (microstate with lowest energy), a dictionary with 
				fixed conformer id list and free conformer id list
		'''
		ms_start=kwargs.get("ms_start",None)
		ms_end=kwargs.get("ms_end", None)

	
		with open(file_path, "rb") as f:
			#read head3.lst, to read information for each conformer
			try: #succeed to read conformer information from head3.lst
				self.conformers = read_headlst()
				crg_flag=1
			except: #failed
				crg_flag=0

			#read number of residue, residue names, method
			byte = f.read(4)
			n_spe = struct.unpack('1i', byte)[0]
			print("There are %i residues." % (n_spe))
			spe_list=""
			for i in range(n_spe):
				byte = f.read(8)
				res_spe= struct.unpack('8s', byte)[0]
				spe_list = str(spe_list) + str(res_spe) + ","
			print("Microstates on "+str(spe_list))
		
			self.enum_flag=0  #only for monte carlo sampling method
		
			#read information for each microstate
			#bytes for each microstates
			ms_bytes=2*n_spe+8+8+4+self.enum_flag*4
			if ms_start and ms_end: 
				f.seek(ms_start*ms_bytes,1)   #move to the offset postion in the file to read microstate
				count = ms_start 
			else:
				count=0
			byte = f.read(2)
			while byte:
				ms_state =MSRECORD()
				confid_list=[]

				#read first conformer id of the microstate
				ms_state.id = count
				res_confid = struct.unpack('1H', byte)[0]
				confid_list.append(res_confid)
			
				#read the rest of residue conformer id of the microstate 
				for i in range(n_spe-1):  
					byte=f.read(2)
					res_confid = struct.unpack('1H', byte)[0]
					confid_list.append(res_confid)
				ms_state.confseq = confid_list
				byte = f.read(8)
				ms_state.cumE = struct.unpack('1d', byte)[0]
				byte = f.read(8)
				ms_state.cumEsq = struct.unpack('1d', byte)[0]

				if (self.enum_flag):
					byte=f.read(8)
					ms_state.occ= struct.unpack('1d', byte)[0]
				else:
					byte=f.read(4)
					ms_state.counter=struct.unpack('1i', byte)[0]
			
				#calculate state energy of the microstate
				ms_state.E = ms_state.cumE/ms_state.counter  
			
				#assign charge state for the microstate
				if crg_flag==0:
					ms_state.crg=None
					ms_state.confidseq=None
					ms_state.crgseq=None
				else:
					ms_state.crg = 0   #store total charge of the microstate
					ms_state.confidseq=[] #store conformer id list of the microstate
					ms_state.crgseq=""    #store residue charge (string) of the microstate
					for index in ms_state.confseq:
						ms_state.crg += self.conformers[int(index)].crg  
						ms_state.crgseq += str(int(self.conformers[int(index)].crg))
						confid = self.conformers[int(index)].id
						ms_state.confidseq.append(confid)
				count += 1
				self.ms_states.append(ms_state)
				#check
				if ms_end and count > ms_end: break

				#read next microstate first conformer id
				byte=f.read(2)
		return self.ms_states


	def read_new_ms(self, file_path ):
		'''
		read ms_out/pH#eH#ms.txt file, a txt file of microstate file, into ms_states
		input: ms_out/pH#eHms.txt, format version: (/home/cai/Github/Stable-MCCE/bin/mcce)
		output: ms_state, a list
				lowest_state (microstate with lowest energy): a dictionary
		'''
		n_lines = 0
		lowest_E=float('inf')
		lowest_state={'fixed_conf_ids':[], 'free_conf_ids':[]}
		with open(file_path, 'r') as in_file:
			line = in_file.readline()
			while line:
				#opt out lines starting with "#" (comments) or blank lines
				while (len(line.rstrip()) ==0 or line[0] == '#'): 
					#print('delete',line)
					line = in_file.readline()

				n_lines += 1
				#print("--",line,"--")

				#read line 1: ph, eh, temp
				if (n_lines == 1):
					tmp= re.split(',|:', line)
					T, pH, eH = float(tmp[1]), float(tmp[3]), float(tmp[5])
					print('Temperature: {:.2f}, pH: {:.2f}, eH: {:.2f}'.format(T, pH, eH))

				#read line 2: method
				elif (n_lines == 2):
					tmp=line.rstrip().split(':')
					method = tmp[1]
					print('Method: {}'.format(method))

				#read line 3: fixed conformer ids
				elif (n_lines == 3):
					fixed_conf_ids=[]
					n_fixed, fixed_conf_id_str= line.split(':')
					n_fixed = int(n_fixed)
					fixed_conf_ids = list(map(int, fixed_conf_id_str.split()))
					print('n_fixed_residue: {:d}, fixed_conf_ids: {}'.format(n_fixed, fixed_conf_ids))
					
					if len(fixed_conf_ids) > n_fixed:
						raise ValueError('number of fixed conf ids must not be large than n_fixed')

					lowest_state['fixed_conf_ids'] = fixed_conf_ids

				#read line 4: conformer ids of each free residue
				elif (n_lines == 4):
					conf_res_inter={}
					n_free_check = 0
					n_free, free_conf_id_str=line.split(':')
					n_free = int(n_free)
					
					for x in free_conf_id_str.strip(';\n').split(';'):
						for y in x.split():
							conf = int(y)
							conf_res_inter[conf] = n_free_check
						n_free_check += 1
					print('n_free: {}, free_conf_residue_iter: {}'.format(n_free, conf_res_inter))

					if (n_free != n_free_check):
						raise ValueError('number of free residues must be equal to n_free value')

				#read inital state for each MC sampling cycle 
				elif (line[0:2] == 'MC'):
					i_MC = int(line.split(':')[1])

					#print microstate after each monte carlo sampling
					if i_MC > 0:
						print('End microstate of {}-th Monte Carlo sampling: {}'.format(i_MC-1, state))

					line = in_file.readline()
					while (len(line.rstrip()) ==0 or line[0] == '#'): 
						#print('delete',line)
						line = in_file.readline()
					n_lines += 1
					tmp, free_conf_id_str = line.rstrip().split(':')
					if (n_free != int(tmp)):
						raise ValueError('n_free value not match to record in fourth line')
					state = list(map(int,free_conf_id_str.split()))
					print("Start microstate of {}-th Monte Carlo sampling: {}".format(i_MC, state))

					if (len(state) != n_free):
						raise ValueError('number of conf in initial state must equal to n_free value')

				#read each flips for every MC sampling
				else:
					energy, count, flips_str=line.split(',')
					energy, count=float(energy), int(count)
					new_confs = list(map(int, flips_str.split()))
					for conf in new_confs:
						state[conf_res_inter[conf]]=conf
						#print(conf)
					#print(state)
					if energy < lowest_E:
						lowest_E = energy
						lowest_state['free_conf_ids'] = state


				line = in_file.readline()
			print('End microstate of {}-th Monte Carlo sampling: {}'.format(i_MC,state))
			print('\nState with lowest energy Emin {:.2f} is:\nFixed_conf_ids: {}\nFree_conf_ids: {}'\
				.format(lowest_E, lowest_state['fixed_conf_ids'], lowest_state['free_conf_ids']))



if __name__ == '__main__':
	if len(sys.argv) < 3:
		print("./read_ms_processor.py old_version/new_version ms_dat_file")
		sys.exit()
	else:
		version=sys.argv[1]
		file=sys.argv[2]

	#read microstates
	read_ms_processor= read_ms_processor()
	if version == 'old_version':
		ms_states = read_ms_processor.read_old_ms(file_path=file)
	elif version ==  'new_version':
		ms_states = read_ms_processor.read_ms(file_path= file)
	elif version == 'Stable-MCCE':
		read_ms_processor.read_new_ms(file_path= file)
	else: 
		ms_states = read_ms_processor.read_re_ms(file_path=file)
	print("Microstate movements: {0}".format(len(ms_states)))
	


