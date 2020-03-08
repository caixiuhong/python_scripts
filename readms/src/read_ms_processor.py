#!/usr/bin/python
import sys
from mfe_adv import read_headlst, read_file_respect_comments
import os
import struct

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

	def read_old_ms(self, file_path='ms.dat'):
		'''
		read ms.dat file, a binary file of microstate file, into ms_states
		input: ms.dat, format version: (/home/cai/mcce3.5)
		output: ms_state, a list
		'''
	
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
			
				#read next microstate first conformer id
				byte=f.read(2)
		return self.ms_states



if __name__ == '__main__':
	if len(sys.argv) < 3:
		print("./readme.py old_version/new_version ms_dat_file")
		sys.exit()
	else:
		version=sys.argv[1]
		file=sys.argv[2]

	#read microstates
	if version_ms_dat == 'old_version':
		read_old_ms(file_path=file)
	elif version_ms_dat == 'new_version':
		read_ms(file_path= file)
	else: 
		read_re_ms(file_path=file)
	


