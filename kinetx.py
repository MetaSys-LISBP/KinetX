################
### kinetX ###
################
# 08/06/2018
# kinetX.py
# version 0.4
#
# This script processes series of 1D spectra (acquired as a pseudo 2D and preprocessed
# with splitserph) to extract changes of chemical shifts and intensity, e.g. used for
# pH titration.
#
####################
### Changes log  ###
####################
#
#    v0.4 (08/06/2018)
#        - new feature: integration routine
#        - new feature: extraction & processing of rows (EFP(), APKS())
#                       (splitser is not a dependency anymore)
#    v0.3 (25/05/2018)
#        - clean code
#    v0.2 (15/03/2018)
#        - if no peak can be found in a slice, returns nan (instead of an error)
#    v0.1 (22/02/2018)
#        - forked from interact v1.1
#        - ...
#
####################
#
#    Author: Pierre Millard, pierre.millard@insa-toulouse.fr
#    Copyright 2018, INRA
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################


##############################
# import modules & functions #
##############################

import os, sys, re


def multiple_replace(text, drep):
    rx = re.compile('|'.join(map(re.escape, drep)))
    def one_xlat(match):
        return drep[match.group(0)]
    res = rx.sub(one_xlat, text)
    return(res)

def load_plist(f_in):
	d_in = dict()
	f = open(f_in, 'r')
	for l in f:
		if l[0] != "#":
			k_v = l.strip('\n').split('\t')
			d_in[k_v[0]] = [float(k) for k in k_v[1:]]
	f.close()
	return(d_in)

def write_ppres(res_folder, slist, ndims, arg):
	if arg == "w":
	  outpp = "Peak picking results\n\n"
	  outpp += "PeakName\tPeakID\tExpno\tSlice\tF1 (ppm)\tIntensity\tresF1\tintegral\n"
	  outpp += "\n".join("\t".join(str(j) for j in i) for i in slist)
	elif arg == "a":
	  outpp = "\n"+"\n".join("\t".join(str(j) for j in i) for i in slist)
	fout = os.path.join(res_folder, "results.txt").replace("\\", "/")
	f = open(fout, arg)
	f.write(outpp)
	f.close()
	if arg == "a":
		outpp = "Peak picking results\n\n"+\
		"PeakName\tPeakID\tExpno\tSlice\tF1 (ppm)\tIntensity\tresF1\tintegral"+\
		outpp
	return(outpp)

def peak_picking1D(F1m, F1p):
  # update peak picking window & pick the most intense peak
  PUTPAR("F1P", str(F1p))
  PUTPAR("F2P", str(F1m))
  fullrange = putil.DataChecks.getNMRDataOfSelectedFramePrintMsg().getFullPhysicalRange()
  newRange = fullrange
  newRange[0].setStart(F1p)
  newRange[0].setEnd(F1m)
  newRange[0].setUnit("ppm")
  XCMD(".zx", 1, newRange)
  # run peak picking silently
  XCMD("pps")

def updateXML(icurdata, pk, po, ndims):
  pk_list_f = os.path.join(icurdata[3], icurdata[0], icurdata[1], 'pdata', icurdata[2], 'peaklist.xml').replace("\\", "/")
  searched_field = '<Peak1D F1="%0.6f"' % (po)
  upd_pk_list_f = []
  f = open(pk_list_f, 'r')
  for l in f.readlines():
  	if searched_field in l and "annotation=" not in l:
  	  upd_l = searched_field + ' annotation="' + pk + '"'
  	  l = l.replace(searched_field, upd_l)
  	upd_pk_list_f.append(l)
  f.close()
  f = open(pk_list_f, 'w')
  f.write(''.join(upd_pk_list_f))
  f.close()
    
def parse_intrng(f_intrng, ppm):
	try:
		f=open(f_intrng, 'r')
		for l in f.readlines():
		  lint = l.split(" ")
		  lint = [i for i in lint if i!=""]
		  try:
		  	if (float(lint[2]) < float(ppm) < float(lint[1])):
		  	  vint = lint[3].strip("\n")
		  except:
		  	pass
		f.close()
		return(vint)
	except:
	  return(None)

def createDic1D():
	try:
	  XCMD("dpl")
	  F1p = float(GETPAR("F1P", axis = 2))
	  F1m = float(GETPAR("F2P", axis = 2))
	  result = INPUT_DIALOG("kinetX", "Peak name.", ["Annotation = "], ["Ala_227"], [""], ["1"])
	  d_in = {result[0]:[F1m, F1p]}
	  return(d_in)
	except:
	  return(None)

##################
# run processing #
##################

# check if multiple display is active
if SELECTED_WINDOW().isMultipleDisplayActive():
  ERRMSG(message = "Please exit multiple display before running kinetX.", title="Error", details=None, modal=1)
  EXIT()

# check spectrum dimension
if GETPROCDIM() != 2:
  ERRMSG(message = "The spectrum to process must have 2 dimensions.", title="Error", details=None, modal=1)
  EXIT()

# create output directories (res & tmp)
current_dataset = CURDATA()
res_folder = os.path.join(current_dataset[3], current_dataset[0], "res").replace("\\", "/")
if not os.path.exists(res_folder):
	os.mkdir(res_folder)
	arg = "w"
else:
	# update existing results files ?
	if os.path.isfile(os.path.join(res_folder, "results.txt").replace("\\", "/")):
		value = SELECT("kinetX", "A result file already exists.\nOverwrite existing result file or append the new results?", ["Append", "Overwrite", "Cancel"])
		if value == 0:
			arg = "a"
		elif value == 1:
			arg = "w"
		elif value == 2:
			EXIT()
	else:
		arg = "w"

# extract procnos
# get td
try:
	td = GETPAR("TD", 1)
except:
	ERRMSG(message = "An unknown error has occured.\nPlease verify that the 2D spectrum is correct and has been processed.", title="Error", details=None, modal=1)
	EXIT()
# select rows to extract
extract=False
value = SELECT("kinetX", "Extract rows of the 2D spectrum?\n(into procnos > 1000)", ["Yes", "No", "Cancel"])
if value == 0:
	extract = True
	ask_n = INPUT_DIALOG("kinetX", "", ["First row"], ["1"], [""], ["1"])
	ask_n2 = INPUT_DIALOG("kinetX", "", ["Last row"], [td], [""], ["1"])
	try:
		if (int(ask_n[0])<1 or int(ask_n[0])>int(td)-1 or int(ask_n2[0])<1 or int(ask_n[0])>int(ask_n[0])):
			ERRMSG(message = "An unknown error has occured.\nPlease verify that the 2D spectrum is correct and has been processed.", title="Error", details=None, modal=1)
			EXIT()
		else:
			lr = range(int(ask_n[0]), int(ask_n2[0])+1, 1)
	except:
		ERRMSG(message = "An unknown error has occured.\nPlease verify that the 2D spectrum is correct and has been processed.", title="Error", details=None, modal=1)
		EXIT()
elif value == 2:
	EXIT()
# extract rows
if extract==True:
	for i in lr:
		try:
		  RSER(str(i), str(999+i), show="n")
		except:
		  ERRMSG(message = "An unknown error has occured with row " + str(i) + ",\nplease verify the 2D spectrum.", title="Error", details=None, modal=1)

# run peack picking ?
runPP = True

# get current dataset & list expnos
current_dataset = CURDATA()
expnos = os.listdir(os.path.join(current_dataset[3], current_dataset[0]).replace("\\", "/")) 
expnos_filt=[]
for i in expnos:
	if i!="res":
		try:
			ii=int(i)
			if ii>=1000:
				expnos_filt.append(i)
		except:
			pass
expnos = expnos_filt
if len(expnos)<1:
	ERRMSG(message = "No data to process.", title="Error", details=None, modal=1)
	EXIT()
	


##########################################
#### 1D spectra processing routine
##########################################

# run peak picking)
if runPP == True:
  # get remaining arguments (expected to be peak database)
  la = [sys.argv[i] for i in range(1, len(sys.argv)) if sys.argv[i][0:2] != "--"]
  # display window to enter parameters for a single signal
  if len(la) == 0:
        d_in = createDic1D()
        if d_in == None:
          EXIT()
  # load the database of peaks to process
  elif len(la) == 1:
        f_in = la[0]
        if os.path.isfile(f_in):
                try:
                        d_in = load_plist(f_in)
                except:
                        ERRMSG(message = "Error when reading file '" + f_in + "'.", title="Error", details=None, modal=1)
                        EXIT()
        else:
                ERRMSG(message = "File '" + f_in + "' not found.", title="Error", details=None, modal=1)
                EXIT()
  elif len(la) > 1:
        ERRMSG(message = "Too many arguments.", title="Error", details=None, modal=1)
        EXIT()
        # if d_in is empty
        if len(d_in) == 0:
          ERRMSG(message = "No signal to process.", title="Error", details=None, modal=1)
          EXIT()
  # get the number of expnos to process
  ask_n = INPUT_DIALOG("kinetX", "", ["First row"], ["1"], [""], ["1"])
  if ask_n == None:
        EXIT()
  try:
          i_min = int(ask_n[0])-1
  except:
          ERRMSG(message = "Row to process must be a positive integer.", title="Error", details=None, modal=1)
          EXIT()
  if i_min < 0:
          ERRMSG(message = "Row to process must be a positive integer.", title="Error", details=None, modal=1)
          EXIT()
  ask_n = INPUT_DIALOG("kinetX", "", ["Last row"], [str(len(expnos))], [""], ["1"])
  if ask_n == None:
        EXIT()
  try:
          i_max = int(ask_n[0])-1
  except:
          ERRMSG(message = "Row to process must be a positive integer.", title="Error", details=None, modal=1)
          EXIT()
  if i_max < 1:
          ERRMSG(message = "Row to process must be a positive integer.", title="Error", details=None, modal=1)
          EXIT()
  if i_max > len(expnos)-1:
          ERRMSG(message = "Row number too high (max=)"+len(expnos)+".", title="Error", details=None, modal=1)
          EXIT()
  # ask for confirmation
  val = CONFIRM("Ok", "Process the following expnos ?\n" + "\n".join(x for x in expnos[i_min:i_max]))
  if val == None:
        EXIT()
  # run peak picking & spectrum annotation
  icurdata = current_dataset
  icurexpno = icurdata[1]
  slist = []
  # for each peak to process
  for pk,v in d_in.items():
        # get the window boundary
        F1m, F1p = v[0], v[1]
        for expnoi in expnos[i_min:i_max]:
          SHOW_STATUS("Processing peak '" + pk + "' in expno '" + str(expnoi) + "'...")
          # load expno (ask for the procno if several procno exists)
          icurdata[1] = expnoi
          RE(icurdata, show="y")
          EFP()
          # get the slice from the expno name
          slice = int(expnoi[1:])
          try:
                peak_picking1D(F1m, F1p)
                listp = GETPEAKSARRAY()
          except:
                listp = None
          if listp == None:
                slist.append([pk, "none", expnoi, slice] + ["none"]*4)
          else:
                # get the last peak picked, with highest intensity
                spec = GETPROCDATA(F1m, F1p)
                ppm = spec.index(max(spec))
                interv = (F1p-F1m)/(len(spec)-1)
                ppm_pk = (F1p-ppm*interv)
                dpol = [abs(peak.getPositions()[0]-ppm_pk) for peak in listp]
                if len(dpol) != 0:
                        idpk = dpol.index(min(dpol))
                        if dpol[idpk] > 2*interv:
                                slist.append([pk, "none", expnoi, slice] + ["none"]*4)
                        else:
                                peak = listp[idpk]
                                # get chemical shifts
                                # as a reminder, to list all peak attributes: dir(peak)
                                po = peak.getPositions()[0]
                                # get integral
                                XCMD("abs")
                                XCMD(".int")
                                XCMD(".sret")
                                f_intrng = os.path.join(icurdata[3], icurdata[0], icurdata[1], 'pdata', icurdata[2], 'integrals.txt').replace("\\", "/")
                                pint = parse_intrng(f_intrng, po)
                                # get resolution
                                reso = peak.getHalfWidth()
                                # append results
                                slist.append([pk, peak.getPeakID()+1, expnoi, slice, po, peak.getRealIntensity(), reso, pint])
                                # update topspin annotation by modifying the xml file
                                updateXML(icurdata, pk, po, 1)
                else:
                        slist.append([pk, float('nan'), expnoi, slice, float('nan'), float('nan'), float('nan'), float('nan')])
  # go back to the initial dataset
  RE(current_dataset, show="y")
  # save pp results
  outpp = write_ppres(res_folder, slist, 1, arg)
  out_all = outpp + "\n\n"
  # done
  SHOW_STATUS("Peak picking finished")
# reload the initial 2D spectrum
RE(current_dataset, show="y")

SHOW_STATUS("Processing finished")

# display results
VIEWTEXT(title="kinetX report", header="results", text=out_all)
