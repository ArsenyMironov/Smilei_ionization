from .Diagnostic import Diagnostic
from .._Utils import *

class Screen(Diagnostic):
	"""Class for loading a Screen diagnostic"""
	
	def _init(self, diagNumber=None, timesteps=None, subset=None, sum=None, data_log=False, **kwargs):
		
		if diagNumber is None:
			self._error += "Printing available screens:\n"
			self._error += "---------------------------\n"
			diags = self.getDiags()
			for diagNumber in diags:
				self._error += self._printInfo(self._getInfo(diagNumber))
			if len(diags)==0:
				self._error += "      No screens found"
			return
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check the requested diags are ok
		if type(diagNumber) is int:
			if diagNumber<0:
				self._error = "Argument 'diagNumber' cannot be a negative integer."
				return
			self.operation = '#' + str(diagNumber)
		elif type(diagNumber) is str:
			self.operation = diagNumber
		else:
			self._error = "Argument 'diagNumber' must be and integer or a string."
			return
		
		# Get list of requested diags
		self._myinfo = {}
		self._diags = sorted(set([ int(d[1:]) for d in self._re.findall('#\d+',self.operation) ]))
		for d in self._diags:
			try:
				info = self._getInfo(d)
				if info is False: raise
				self._myinfo.update({ d:info })
			except:
				self._error = "Screen #"+str(d)+" invalid"
				return
		try:
			exec(self._re.sub('#\d+','1.',self.operation))
		except ZeroDivisionError: pass
		except:
			self._error = "Cannot understand operation '"+self.operation+"'"
			return
		# Verify that all requested diags all have the same shape
		self._axes = {}
		self._naxes = {}
		for d in self._diags:
			self._axes .update ({ d:self._myinfo[d]["axes"] })
			self._naxes.update ({ d:len(self._axes[d]) })
			if self._naxes[d] != self._naxes[self._diags[0]]:
				self._error = "All diagnostics in operation '"+self.operation+"' must have as many axes." \
					+ " Diagnotic #"+str(d)+" has "+str(self._naxes[d])+" axes and #"+ \
					str(self._diags[0])+" has "+str(self._naxes[self._diags[0]])+" axes"
				return
			for a in self._axes[d]:
				if self._axes[d] != self._axes[self._diags[0]]:
					self._error = "In operation '"+self.operation+"', diagnostics #"+str(d)+" and #"\
						+str(self._diags[0])+" must have the same shape."
					return
		
		self._axes  = self._axes [self._diags[0]]
		self._naxes = self._naxes[self._diags[0]]
		
		# Check subset
		if subset is None: subset = {}
		elif type(subset) is not dict:
			self._error = "Argument `subset` must be a dictionary"
			return
		
		# Check sum
		if sum is None: sum = {}
		elif type(sum) is not dict:
			self._error = "Argument 'sum' must be a dictionary"
			return
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# Get available timesteps
		self.times = {}
		self._times = {}
		self._indexOfTime  = {}
		self._h5items = {}
		for d in self._diags:
			# Gather data from all timesteps, and the list of timesteps
			items = {}
			for path in self._results_path:
				f = self._h5py.File(path+self._os.sep+'Screen'+str(d)+'.h5')
				items.update( dict(f) )
			items = sorted(items.items())
			self._h5items[d] = [it[1] for it in items]
			self.times[d] = self._np.array([ int(it[0].strip("timestep")) for it in items ])
			self._times[d] = self._np.copy(self.times[d])
			# fill the "_indexOfTime" dictionary with indices to the data arrays
			self._indexOfTime.update({ d:{} })
			for i,t in enumerate(self.times[d]):
				self._indexOfTime[d].update({ t : i })
			# If timesteps is None, then keep all timesteps, otherwise, select timesteps
			if timesteps is not None:
				try:
					self.times[d] = self._selectTimesteps(timesteps, self.times[d])
				except:
					self._error = "Argument 'timesteps' must be one or two non-negative integers"
					return
			# Verify that timesteps are the same for all diagnostics
			if (self.times[d] != self.times[self._diags[0]]).any() :
				self._error = "All screens in operation '"+self.operation+"' must have the same timesteps."\
					+" Screen #"+str(d)+" has "+str(len(self.times[d]))+ " timesteps and #"\
					+str(self._diags[0])+" has "+str(len(self.times[self._diags[0]]))+ " timesteps"
				return
		# Now we need to keep only one array of timesteps because they should be all the same
		self.times  = self.times [self._diags[0]]
		self._times = self._times[self._diags[0]]
		
		# Need at least one timestep
		if self.times.size < 1:
			self._error = "Timesteps not found"
			return
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values for all diags
		plot_diff = []
		cell_volume = self._cell_length.prod()
		coeff = 1.
		unitsa = [0,0,0,0]
		spatialaxes = {"x":False, "y":False, "z":False}
		self._finalShape = [[]]*self._naxes
		self._sums = [False]*self._naxes
		self._selection = [self._np.s_[:]]*self._naxes
		hasComposite = False
		
		for iaxis in range(self._naxes):
			axis = self._axes[iaxis]
			axistype = axis["type"]
			
			# Find the vector of values along the axis
			if axis["log"]:
				edges = self._np.linspace(self._np.log10(axis["min"]), self._np.log10(axis["max"]), axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				edges = 10.**edges
				centers = 10.**(centers[:-1])
			else:
				edges = self._np.linspace(axis["min"], axis["max"], axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				centers = centers[:-1]
			axis.update({ "edges"   : edges   })
			axis.update({ "centers" : centers })
			
			# Find some quantities depending on the axis type
			overall_min = "-inf"; overall_max = "inf"
			axis_units = ""
			if   axistype in ["x","y","z","moving_x"]:
				axis_units = "L_r"
				spatialaxes[axistype[-1]] = True
			elif axistype in ["a","b"]:
				axis_units = "L_r"
				hasComposite = True
			elif axistype == "theta" and self._ndim==2:
				axis_units = "rad"
				overall_min = "-3.141592653589793"
				overall_max = "3.141592653589793"
			elif axistype == "theta" and self._ndim==3:
				axis_units = "rad"
				overall_min = "0"
				overall_max = "3.141592653589793"
			elif axistype == "phi":
				axis_units = "rad"
				overall_min = "-3.141592653589793"
				overall_max = " 3.141592653589793"
			elif axistype[:9] == "composite":
				axis_units = "L_r"
				hasComposite = True
				axistype = axistype[10:]
				axis["type"] = axistype
			elif axistype in ["px","py","pz","p"]:
				axis_units = "P_r"
			elif axistype in ["vx","vy","vz","v"]:
				axis_units = "V_r"
			elif axistype in ["vperp2"]:
				axis_units = "V_r**2"
				overall_min = "0"
			elif axistype == "gamma":
				overall_min = "1"
			elif axistype == "ekin":
				axis_units = "K_r"
				overall_min = "0"
			elif axistype == "charge":
				axis_units = "Q_r"
				overall_min = "0"
			else:
				self._error = "axis type "+axistype+" not implemented"
				return None
			
			# if this axis has to be summed, then select the bounds
			if axistype in sum:
				if axistype in subset:
					self._error = "`subset` not possible on the same axes as `sum`"
					return
				
				self._sums[iaxis] = True
				
				try:
					axis["sumInfo"], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectRange(sum[axistype], centers, axistype, axis_units, "sum")
				except:
					return
				
				if axistype in ["x","y","z","moving_x"]:
					first_edge = edges[self._selection[iaxis].start or 0]
					last_edge  = edges[(self._selection[iaxis].stop or len(centers))]
					coeff /= last_edge - first_edge
			
			# if not summed
			else:
				# If taking a subset of this axis
				if axistype in subset:
					try:
						axis["subsetInfo"], self._selection[iaxis], self._finalShape[iaxis] \
							= self._selectSubset(subset[axistype], centers, axistype, axis_units, "subset")
					except:
						return
				# If subset has more than 1 point (or no subset), use this axis in the plot
				if type(self._selection[iaxis]) is slice:
					self._type   .append(axistype)
					self._shape  .append(axis["size"])
					self._centers.append(centers[self._selection[iaxis]])
					self._log    .append(axis["log"])
					self._label  .append(axistype)
					self._units  .append(axis_units)
					if axistype == "theta" and self._ndim==3:
						plot_diff.append(self._np.diff(self._np.cos(edges))[self._selection[iaxis]])
					else:
						plot_diff.append(self._np.diff(edges)[self._selection[iaxis]])
					self._finalShape[iaxis] = len(self._centers[-1])
		
		self._selection = tuple(self._selection)
		
		# Build units
		titles = {}
		units = {}
		for d in self._diags:
			titles.update({ d:"??" })
			units.update({ d:"??" })
			val_units = "??"
			output = self._myinfo[d]["output"]
			if   output == "density":
				titles[d] = "Number" + ("" if hasComposite else " density")
				val_units = "1" if hasComposite else "N_r"
			elif output == "charge_density":
				titles[d] = "Charge" + ("" if hasComposite else " density")
				val_units = "Q_r" if hasComposite else "N_r * Q_r"
			elif output=="ekin_density":
				titles[d] = "Energy" + ("" if hasComposite else " density")
				val_units = "K_r" if hasComposite else "N_r * K_r"
			elif output[0] == "j":
				titles[d] = "J"+output[1] + (" x Volume" if hasComposite else "")
				val_units = "J_r/N_r" if hasComposite else "J_r"
			elif output[0]=="p" and output[-8:] == "_density":
				titles[d] = "P"+output[1].strip("_") + ("" if hasComposite else " density")
				val_units = "P_r" if hasComposite else "N_r * P_r"
			elif output[:8]=="pressure":
				titles[d] = "Pressure "+output[-2] + (" x Volume" if hasComposite else "")
				val_units = "K_r" if hasComposite else "N_r * K_r"
			axes_units = [unit for unit in self._units if (hasComposite or unit!="L_r")]
			units[d] = val_units
			if len(axes_units)>0: units[d] += " / ( " + " * ".join(axes_units) + " )"
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for d in self._diags:
			self._vunits = self._vunits.replace("#"+str(d), "( "+units[d]+" )")
			self._title  = self._title .replace("#"+str(d), titles[d])
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if self._ndim>=1 and not spatialaxes["x"]: coeff /= self._ncels[0]*self._cell_length[0]
		if self._ndim>=2 and not spatialaxes["y"]: coeff /= self._ncels[1]*self._cell_length[1]
		if self._ndim==3 and not spatialaxes["z"]: coeff /= self._ncels[2]*self._cell_length[2]
		
		# Calculate the array that represents the bins sizes in order to get units right.
		# This array will be the same size as the plotted array
		if len(plot_diff)==0:
			self._bsize = 1.
		elif len(plot_diff)==1:
			self._bsize = plot_diff[0]
		else:
			self._bsize = self._np.prod( self._np.array( self._np.meshgrid( *plot_diff ) ), axis=0)
			self._bsize = self._bsize.transpose([1,0]+range(2,len(plot_diff)))
		self._bsize = cell_volume / self._bsize
		if not hasComposite: self._bsize *= coeff
		self._bsize = self._np.squeeze(self._bsize)
		
		# Set the directory in case of exporting
		self._exportPrefix = "Screen_"+"-".join([str(d) for d in self._diags])
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Gets info about diagnostic number "diagNumber"
	def _getInfo(self,diagNumber):
		info = {}
		for path in self._results_path:
			# Open file
			try:
				file = path+self._os.sep+'Screen'+str(diagNumber)+'.h5'
				f = self._h5py.File(file, 'r')
			except:
				return False
			# get attributes from file
			attrs = f.attrs.items()
			axes = []
			# Parse each attribute
			for name, value in attrs:
				if (name == "output"): output = bytes.decode(value)
				if (name == "species"):
					species = bytes.decode(value.strip()).split() # get all species numbers
					species = [int(s) for s in species]
				if (name[0:4] == "axis" ):
					n = int(name[4:]) # axis number
					sp = bytes.decode(value).split()
					axistype  = sp[0]
					axismin  = float(sp[1])
					axismax  = float(sp[2])
					axissize = int(sp[3])
					logscale = bool(int(sp[4]))
					edge_inclusive = bool(int(sp[5]))
					coefficients = 0
					if sp[6]!="[]" :
						coefficients = eval(sp[6])
					while len(axes)<n+1: axes.append({}) # extend the array to the adequate size
					axes[n] = {"type":axistype,"min":axismin,"max":axismax,"size":axissize,"log":logscale,"edges_included":edge_inclusive,"coefficients":coefficients}
			f.close()
			# Verify that the info corresponds to the diag in the other paths
			if info == {}:
				info = {"#":diagNumber, "output":output, "species":species, "axes":axes}
			else:
				if output!=info["output"] or axes!=info["axes"]:
					print("Screen #"+str(diagNumber)+" in path '"+path+"' is incompatible with the other ones")
					return False
		return info
	
	
	# Prints the info obtained by the function "getInfo"
	@staticmethod
	def _printInfo(info):
		# 1 - diag number, type and list of species
		species = ""
		for i in range(len(info["species"])): species += str(info["species"][i])+" " # reconstitute species string
		printedInfo = "Screen#"+str(info["#"])+" - "+info["output"]+" of species # "+species+"\n"
		
		# 3 - axes
		for i in range(len(info["axes"])):
			axis = info["axes"][i];
			logscale = "" if not axis["log"] else " [ LOG SCALE ] "
			edges    = "" if not axis["edges_included"] else " [ INCLUDING EDGES ] "
			printedInfo += "    "+axis["type"]+" from "+str(axis["min"])+" to "+str(axis["max"]) \
				   +" in "+str(axis["size"])+" steps "+logscale+edges+"\n"
		return printedInfo
		
	# Method to print info on all included diags
	def _info(self):
		info = ""
		for d in self._diags:
			info += self._printInfo(self._myinfo[d])+"\n"
		if len(self.operation)>2: info += "Operation : "+self.operation+"\n"
		for ax in self._axes:
			if "sumInfo" in ax: info += ax["sumInfo"]+"\n"
			if "subsetInfo" in ax: info += ax["subsetInfo"]+"\n"
		return info
	
	def getDiags(self):
		allDiags = []
		for path in self._results_path:
			files = self._glob(path+self._os.sep+"Screen*.h5")
			diags = [int(self._re.findall(r"Screen([0-9]+)[.]h5$",file)[0]) for file in files]
			if len(allDiags)>0:
				allDiags = [d for d in diags if d in allDiags]
			else:
				allDiags = diags
		return allDiags
	
	# get all available timesteps for a given diagnostic
	def getAvailableTimesteps(self, diagNumber=None):
		# if argument "diagNumber" not provided, return the times calculated in __init__
		if diagNumber is None:
			return self._times
		# Otherwise, get the timesteps specifically available for the single requested diagnostic
		else:
			times = set()
			for path in self._results_path:
				try:
					file = path+self._os.sep+'Screen'+str(diagNumber)+'.h5'
					f = self._h5py.File(file, 'r')
				except:
					print("Cannot open file "+file)
					return self._np.array([])
				times.update( set(f.keys()) )
				f.close()
			times = [int(t.strip("timestep")) for t in times]
			return self._np.array(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Get arrays from all requested diagnostics
		A = {}
		for d in self._diags:
			# find the index of the array corresponding to the requested timestep
			try:
				index = self._indexOfTime[d][t]
			except: 
				print("Timestep "+str(t)+" not found in this screen")
				return []
			# get data
			B = self._np.empty(self._finalShape)
			try:
				self._h5items[d][index].read_direct(B, source_sel=self._selection) # get array
			except:
				B = self._np.squeeze(B)
				self._h5items[d][index].read_direct(B, source_sel=self._selection) # get array
				B = self._np.reshape(B, self._finalShape)
			B[self._np.isnan(B)] = 0.
			# Apply the summing
			for iaxis in range(self._naxes):
				if self._sums[iaxis]:
					B = self._np.sum(B, axis=iaxis, keepdims=True)
			# remove summed axes
			B = self._np.squeeze(B)
			# Divide by the bins size
			B *= self._bsize
			# Append this diag's array for the operation
			A.update({ d:B })
		# Calculate operation
		data_operation = self.operation
		for d in reversed(self._diags):
			data_operation = data_operation.replace("#"+str(d),"A["+str(d)+"]")
		A = eval(data_operation)
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
