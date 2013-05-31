"""
The Linac Lattice Factory generates the Linac Accelerator Lattice from the information
inside of the LinacStructureTree instance which in turn was generated by the LinacParser.
The Linac Lattice Factory uses a predefined set of Linac Acc Elements. If you do not have 
the LinacStructureTree instance you can create the Linac Accelerator Lattice directly in 
the script.
"""

import os
import math

# import the linac structure tree with all sequences and nodes, but without drifts
from LinacParser import LinacStructureTree

from LinacAccNodes import BaseLinacNode, LinacNode, LinacMagnetNode, MarkerLinacNode, Drift, Quad, BaseRF_Gap, Bend
from LinacAccNodes import DCorrectorH, DCorrectorV
from LinacAccNodes import RF_Cavity, Sequence

from LinacAccLattice import LinacAccLattice

# import general accelerator elements
from orbit.lattice import AccNode

# import pyORBIT Python utilities classes for objects with names, types, and dictionary parameters
from orbit.utils import orbitFinalize

class LinacLatticeFactory():
	""" 
	The Linac Lattice Factory generates the Linac Accelerator Lattice 
	from the information inside of the LinacStructureTree instance. 
	"""
	def __init__(self, ltree):
		if(isinstance(ltree, LinacStructureTree) != True):
			msg = "The LinacLatticeFactory constructor: you have to specify the LinacStructureTree instance as input!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)	
		self.ltree = ltree
		#We need to compare positions, lengths etc. This is our delta
		self.zeroDistance = 0.00001
		#The maximal length of the drift. It will be devided if it is more than that.
		self.maxDriftLength = 1.
		
	def setMaxDriftLength(self, maxDriftLength = 1.0):
		"""
		Sets the maximal drift length that is used for 
		the purpose of the space charge calculations and diagnostics.
		"""
		self.maxDriftLength = maxDriftLength
		
	def getMaxDriftLength(self):
		"""
		Returns the maximal drift length that is used for the purpose 
		of the space charge calculations and diagnostics.
		"""
		return self.maxDriftLength
	
	def getLinacAccLattice(self,names):
		"""
		Returns the linac accelerator lattice for specified sequence names.
		"""
		if(len(names) < 1):
			msg = "The LinacLatticeFactory method getLinacAccLattice(names): you have to specify the names array!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)
		#let's check that the names in good order  ==start==
		seqencesLocal = self.ltree.getSeqs()
		seqencesLocalNames = []
		for seq in seqencesLocal:
			seqencesLocalNames.append(seq.getName())
		ind_old = -1
		count = 0
		for name in names:
			ind = seqencesLocalNames.index(name)
			if(ind < 0 or (count > 0 and ind != (ind_old + 1))):
				msg = "The LinacLatticeFactory method getLinacAccLattice(names): sequence names array is wrong!"
				msg = msg + os.linesep
				msg = msg + "existing names=" + str(seqencesLocalNames)
				msg = msg + os.linesep
				msg = msg + "sequence names="+str(names)
				orbitFinalize(msg)
			ind_old = ind
			count += 1
		#	let's check that the names in good order  ==stop==			
		ind_start = seqencesLocalNames.index(names[0])
		sequences = self.ltree.getSeqs()[ind_start:ind_start+len(names)]
		#----make linac lattice
		linacAccLattice = LinacAccLattice(self.ltree.getName())
		#There are the folowing possible types of elements in the linac tree:
		#QUAD - quadrupole
		#RFGAP - RF Gap
		#DCH - horizontal dipole corrector
		#DCV - vertical dipole corrector
		#Marker - anything else with the length equals to 0
		#Before putting enerything into the linacAccLattice we will create sequences 
		# with all nodes.
		#----------------------------------------------------------------------
		# The DRIFTS will be generated additionally and put into right places
		#----------------------------------------------------------------------
		accSeqs = []
		accRF_Cavs = [] 
		seqPosition = 0.
		for seq in sequences:
			#print "debug =========================================== seq=",seq.getName()
			accSeq = Sequence(seq.getName())
			accSeq.setLinacAccLattice(linacAccLattice)
			accSeq.setLength(float(seq.getLength()))
			accSeq.setPosition(seqPosition)
			seqPosition = seqPosition + accSeq.getLength()
			accSeqs.append(accSeq)
			#these nodes are not AccNodes. They are from linac parser
			nodes = seq.getNodes()
			#rf_cav_names is an auxilary array with RF Cav. names
			rf_cav_names = []
			#array of nodes that are AccNodes with zero length
			#They can be positioned inside the thick nodes, and this will be done at the end
			#of this constructor
			thinNodes = []
			for node in nodes:
				node.setParam("pos",float(node.getParam("pos")))
				#------------QUAD-----------------
				if(node.getType() == "QUAD"):
					accNode = Quad(node.getName())
					accNode.updateParamsDict(node.getParamsDict())					
					accNode.setParam("dB/dr",node.getParam("field"))
					accNode.setParam("field",node.getParam("field"))
					accNode.setLength(node.getParam("effLength"))
                                                                                
                                        #if(0.5*accNode.getLength() > self.maxDriftLength):                                        
                                        #    accNode.setnParts(2*int(0.5*int(accNode.getLength()/(self.maxDriftLength + 1.0e-15) + 1  > 0)*int(accNode.getLength()/(self.maxDriftLength + 1.0e-15) + 1 ) + 1))

					accSeq.addNode(accNode)
				#------------BEND-----------------
                                elif(node.getType() == "BEND"):
					accNode = Bend(node.getName())                                        
					accNode.updateParamsDict(node.getParamsDict())	                                        
					accNode.setParam("poles",[int(x) for x in eval(node.getParam("poles"))])
					accNode.setParam("kls", [x for x in eval(node.getParam("kls"))])
                                        accNode.setParam("skews",[int(x) for x in eval(node.getParam("skews"))])
                                        accNode.setParam("ea1",node.getParam("ea1"))
                                        accNode.setParam("ea2",node.getParam("ea2"))
                                        accNode.setParam("theta",node.getParam("theta"))
					accNode.setLength(node.getParam("effLength"))
                                                                                
                                        #if(0.5*accNode.getLength() > self.maxDriftLength):                                        
                                        #    accNode.setnParts(2*int(0.5*int(accNode.getLength()/(self.maxDriftLength + 1.0e-15) + 1  > 0)*int(accNode.getLength()/(self.maxDriftLength + 1.0e-15) + 1 ) + 1))

					accSeq.addNode(accNode)
				#------------RF_Gap-----------------	
				elif(node.getType() == "RFGAP"):
					accNode = BaseRF_Gap(node.getName())
					accNode.updateParamsDict(node.getParamsDict())	
					accNode.setParam("gapOffset",node.getParam("gapOffset"))					
					accNode.setLength(node.getParam("gapLength"))
					accNode.setParam("amp",node.getParam("amp"))
					#the parameter from XAL in MeV, we use GeV
					#accNode.setParam("E0TL",1.0e-3*node.getParam("E0TL"))
					accNode.setParam("E0TL",0.001*node.getParam("E0TL"))
					accNode.setParam("length",node.getParam("gapLength"))
					accNode.setParam("gapLength",node.getParam("gapLength"))		
					accNode.setParam("modePhase",node.getParam("modePhase"))
					rf_cav_name = node.getParam("parentCavity")
					if(rf_cav_name not in rf_cav_names):
						accNode.setParam("firstPhase", (math.pi/180.)*accNode.getParam("firstPhase"))
						rf_cav_names.append(rf_cav_name)
						accRF_Cav = RF_Cavity(rf_cav_name)
						accRF_Cavs.append(accRF_Cav)
						accRF_Cav._setDesignPhase(accNode.getParam("firstPhase"))
						accRF_Cav.setPhase(accNode.getParam("firstPhase"))
						accRF_Cav._setDesignAmp(1.)
						accRF_Cav.setAmp(1.)
						accRF_Cav.setFrequency(seq.getParam("rfFrequency"))
					accRF_Cav = accRF_Cavs[len(accRF_Cavs) - 1] 
					accRF_Cav.addRF_GapNode(accNode)
					accSeq.addNode(accNode)
				else:
					if(node.getParam("length") != 0.):
						msg = "The LinacLatticeFactory method getLinacAccLattice(names): there is a strange element!"
						msg = msg + os.linesep
						msg = msg + "name=" + node.getName()
						msg = msg + os.linesep
						msg = msg + "type="+node.getType()
						msg = msg + os.linesep
						msg = msg + "length(should be 0.)="+str(node.getParam("length"))
						orbitFinalize(msg)						
					thinNodes.append(node)
			#insert the drifts ======================start ===========================
			#-----now check the integrety quads and rf_gaps should not overlap
			#-----and create drifts
			copyAccNodes = accSeq.getNodes()[:]
			firstNode = copyAccNodes[0]
			lastNode = copyAccNodes[len(copyAccNodes)-1]
			#insert the drift before the first element if its half length less than its position
			if(math.fabs(firstNode.getLength()/2.0 - firstNode.getParam("pos")) > self.zeroDistance):
				if(firstNode.getLength()/2.0 > firstNode.getParam("pos")):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the first node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + firstNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + firstNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(firstNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(firstNode.getParam("pos"))						
					orbitFinalize(msg)
				else:
					driftLength = firstNode.getParam("pos") - firstNode.getLength()/2.0
					nDrifts = int(driftLength/self.maxDriftLength) + 1
					driftLength = driftLength/nDrifts
					for idrift_reverse in range(nDrifts):
						idrift = (nDrifts-1) - idrift_reverse
						drift = Drift(accSeq.getName()+":"+firstNode.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",0.+drift.getLength()*(idrift+0.5))
						accSeq.addNode(drift, index = 0)
			#insert the drift after the last element if its half length less + position is less then the sequence length
			if(math.fabs(lastNode.getLength()/2.0 + lastNode.getParam("pos") - accSeq.getLength()) > self.zeroDistance):
				if(lastNode.getLength()/2.0 + lastNode.getParam("pos") > accSeq.getLength()):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the last node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + lastNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + lastNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(lastNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(lastNode.getParam("pos"))					
					msg = msg + os.linesep
					msg = msg + "sequence name=" + accSeq.getName()				
					msg = msg + os.linesep
					msg = msg + "sequence length=" + str(accSeq.getLength())			
					orbitFinalize(msg)
				else:
					driftLength = accSeq.getLength() - (lastNode.getParam("pos") + lastNode.getLength()/2.0)
					nDrifts = int(driftLength/self.maxDriftLength) + 1
					driftLength = driftLength/nDrifts
					for idrift in range(nDrifts):
						drift = Drift(accSeq.getName()+":"+lastNode.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",lastNode.getParam("pos")+lastNode.getLength()/2.0 + drift.getLength()*(idrift+0.5))
						accSeq.addNode(drift)	
			#now move on and generate drifts between (i,i+1) nodes from copyAccNodes
			for node_ind in range(len(copyAccNodes)-1):
				accNode0 = copyAccNodes[node_ind]
				accNode1 = copyAccNodes[node_ind+1]
				ind_of_node =  accSeq.getNodes().index(accNode1)
				dist = accNode1.getParam("pos") - accNode1.getLength()/2 - (accNode0.getParam("pos") + accNode0.getLength()/2)
				if(dist < 0.):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): two nodes are overlapping!"
					msg = msg + os.linesep
					msg = msg + "sequence name=" + accSeq.getName()		
					msg = msg + os.linesep					
					msg = msg + "node 0 name=" + accNode0.getName() + " pos="+ str(accNode0.getParam("pos")) + " L="+str(accNode0.getLength())
					msg = msg + os.linesep
					msg = msg + "node 1 name=" + accNode1.getName() + " pos="+ str(accNode1.getParam("pos")) + " L="+str(accNode1.getLength())			
					msg = msg + os.linesep
					orbitFinalize(msg)				
				elif(dist > self.zeroDistance):
					nDrifts = int(dist/self.maxDriftLength) + 1
					driftLength = dist/nDrifts				
					for idrift in range(nDrifts):
						drift = Drift(accSeq.getName()+":"+accNode0.getName()+":"+str(idrift+1)+":drift")
						drift.setLength(driftLength)
						drift.setParam("pos",accNode0.getParam("pos")+accNode0.getLength()*0.5+drift.getLength()*(idrift+0.5))
						accSeq.addNode(drift, index = ind_of_node+idrift)
				else:
					pass
			#insert the drifts ======================stop ===========================		
			#========================================================================
			#Now we will go over all zero length nodes and attach them into the quads
			#or drifts. We cannot put anything inside RF Cavity.
			# zero length elements insertion ========== start ======================
			# if a zero-length element is inside a quad it will be placed inside this 
			# quad
			accQuads = []
			for accNode in accSeq.getNodes():
				if(isinstance(accNode,Quad)): accQuads.append(accNode)
			usedThinNodes = []
			for node in thinNodes:
				position = node.getParam("pos")
				for quad in accQuads:
					pos = quad.getParam("pos")
					L = quad.getLength()
					nParts = quad.getnParts()
					if(abs(position - pos) < self.zeroDistance or (position > pos - L/2.0 and position < pos +L/2.0)):
						accNode = None
						if(node.getType() == "DCV" or node.getType() == "DCH"):
							if(node.getType() == "DCV"): accNode = DCorrectorV(node.getName())
							if(node.getType() == "DCH"): accNode = DCorrectorH(node.getName())
							accNode.setParam("effLength",float(node.getParam("effLength")))
						else:
							accNode = MarkerLinacNode(node.getName())
							accNode.updateParamsDict(node.getParamsDict())
						accNode.setParam("pos",quad.getParam("pos"))
						quad.addChildNode(accNode, place = AccNode.BODY, part_index = (nParts/2) - 1 , place_in_part = AccNode.AFTER)
						usedThinNodes.append(node)
			#remove all assigned zero-length nodes from list of thin nodes
			for node in usedThinNodes:
				thinNodes.remove(node)
			#----------------
			# chop the drifts if the thin element is inside or insert this element into
			# the sequence at the end or the beginning of the drift
			usedThinNodes = []
			for node in thinNodes:
				position = node.getParam("pos")
				driftNode = self.__getDriftThinNode(position,accSeq)
				if(driftNode != None):
					usedThinNodes.append(node)
					pos = driftNode.getParam("pos")
					L = driftNode.getLength()
					ind_insertion = accSeq.getNodes().index(driftNode)
					accNode = None
					if(node.getType() == "DCV" or node.getType() == "DCH"):
						if(node.getType() == "DCV"): accNode = DCorrectorV(node.getName())
						if(node.getType() == "DCH"): accNode = DCorrectorH(node.getName())
						accNode.setParam("effLength",float(node.getParam("effLength")))
					else:
						accNode = MarkerLinacNode(node.getName())
						accNode.updateParamsDict(node.getParamsDict())
					accNode.setParam("pos",position)					
					if(abs(position - (pos - L/2.0)) < self.zeroDistance):
						#insert before the drift
						accSeq.addNode(accNode, index = ind_insertion)
					elif(abs(position - (pos + L/2.0)) < self.zeroDistance):
						#insert after the drift
						accSeq.addNode(accNode, index = ind_insertion+1)	
					else:
						#we replace this drift with two new
						drift_node_name = driftNode.getName()
						ind_name_pos = drift_node_name.find(":drift")
						drift_node_name = drift_node_name[0:ind_name_pos]
						drift0 = Drift(drift_node_name+":1:drift")
						drift0.setLength(position - (pos - L/2.0))
						drift0.setParam("pos",(pos - L/2.0) + drift0.getLength()/2.0)
						drift1 = Drift(drift_node_name+":2:drift")
						drift1.setLength((pos + L/2.0) - position)
						drift1.setParam("pos",(pos + L/2.0) - drift1.getLength()/2.0)
						accSeq.getNodes().remove(driftNode)
						accSeq.addNode(drift0, index = ind_insertion)
						accSeq.addNode(accNode, index = ind_insertion + 1)	
						accSeq.addNode(drift1, index = ind_insertion + 2)
			#remove all assigned zero-length nodes from list of thin nodes
			for node in usedThinNodes:
				thinNodes.remove(node)			
			if(len(thinNodes) != 0):
				print "==========WARNING!!!!==============="
				print "The seqence =",accSeq.getName()," has nodes that are not assigned to the lattice:"
				for node in thinNodes:
					print "unused node =",node.getName()," pos=",node.getParam("pos")
			# add all AccNodes to the linac lattice
			L_total = 0.
			for accNode in accSeq.getNodes():
				pos = accNode.getParam("pos")
				L = accNode.getLength()
				L_total += L
				linacAccLattice.addNode(accNode)
			linacAccLattice.addSequence(accSeq)
		# zero length elements insertion ========== stop ======================	
		for accRF_Cav in accRF_Cavs:
		 linacAccLattice.addRF_Cavity(accRF_Cav)
		return linacAccLattice
	
	def __getDriftThinNode(self,position,accSeq):
		"""
		This method will return None or the drift AccNode in accSeq which cover this
		position.
		"""
		resNode = None
		ind_start = 0
		ind_stop = len(accSeq.getNodes()) - 1
		while(ind_stop - ind_start > 1):
			ind = (ind_stop + ind_start)/2
			accNode = accSeq.getNodes()[ind]
			pos = accNode.getParam("pos") - accNode.getLength()/2.0
			if(position > pos):
				ind_start = ind
			else:
				ind_stop = ind
		if(isinstance(accSeq.getNodes()[ind_start],Drift)):
			resNode = accSeq.getNodes()[ind_start]
		#check if the last node is the guy
		if(resNode == None):
			accNode = accSeq.getNodes()[len(accSeq.getNodes())-1]
			if(isinstance(accNode,Drift)):
				pos = accNode.getParam("pos") - accNode.getLength()/2.0
				if(pos <= position):
					resNode = accNode
		return resNode


