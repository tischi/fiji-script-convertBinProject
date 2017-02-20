#
# use elastix from Fiji
#
# Author information:
# 
# tischitischer@gmail.com
#
# Input: 
# 
# Computation:
#
# Output:
#
#


from ij.io import OpenDialog
from ij.io import Opener
from fiji.util.gui import GenericDialogPlus
from ij.plugin import ZProjector, RGBStackMerge, SubstackMaker, Concatenator
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin import Duplicator
from ij.process import StackStatistics
from ij.plugin import ImageCalculator
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
import os, os.path, re, sys
from subprocess import Popen, PIPE
from ij.process import ImageConverter
import os, time, shutil, sys, math
from ij.macro import MacroRunner
from ij.gui import Plot
import collections, pickle, platform

from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions

from automic.table import TableModel			# this class stores the data for the table
from automic.table import ManualControlFrame 	# this class visualises TableModel via GUI
from java.io import File
from automic.utils.roi import ROIManipulator2D as ROIManipulator

#
#  Functions
#  

def close_all_image_windows():
  # forcefully closes all open images windows
  ids = WindowManager.getIDList();
  if (ids==None):
    return
  for i in ids:
     imp = WindowManager.getImage(i)
     if (imp!=None):
       win = imp.getWindow()
       if (win!=None):
         imp.changes = False # avoids the "save changes" dialog
         win.close()

def extractChannel(imp, nChannel, nFrame):
  """ Extract a stack for a specific color channel and time frame """
  stack = imp.getImageStack()
  ch = ImageStack(imp.width, imp.height)
  for i in range(1, imp.getNSlices() + 1):
    index = imp.getStackIndex(nChannel, i, nFrame)
    ch.addSlice(str(i), stack.getProcessor(index))
  return ImagePlus("Channel " + str(nChannel), ch)

def measureSumIntensity3D(imp):
  stats = StackStatistics(imp)
  return stats.mean * stats.pixelCount

def autoThreshold(imp, method):
  impout = imp.duplicate() 
  IJ.run(impout, "Auto Threshold", "method=" + method + " white stack use_stack_histogram");
  impout.setTitle("Auto Threshold")
  return impout

def threshold(_img, threshold):
  imp = Duplicator().run(_img)
  #imp.show(); time.sleep(0.2)
  #IJ.setThreshold(imp, mpar['lthr'], mpar['uthr'])
  IJ.setThreshold(imp, threshold, 1000000000)
  IJ.run(imp, "Convert to Mask", "stack")
  imp.setTitle("Threshold");
  #IJ.run(imp, "Divide...", "value=255 stack");
  #IJ.setMinAndMax(imp, 0, 1);
  return imp

def show_standard_error_message():
  IJ.error("There was an error.\n\
Please check the text below the script editor window.\n\
Please toggle between [Show Errors] and [Show Output], as both are relevant.")
    
def rename(old_filepath, new_filepath):
  try:
    os.rename(old_filepath, new_filepath); time.sleep(1)
  except:
    show_standard_error_message()
    print("\n  error during renaming:")
    print("    renaming: "+old_filepath)
    print("        into: "+new_filepath)
    print("  Often this happens because the target file exists already and is write protected; please try again with a new output folder.")
    sys.exit(0)
  

def project(imp, projection):
  IJ.run(imp, "Set Scale...", "distance=1 global");
  
  if projection=='Z':
    IJ.run(imp, "Z Project...", "projection=[Max Intensity]");
    imp_projected = IJ.getImage()
  elif projection=='X':
    IJ.run(imp, "Reslice [/]...", "output=1.000 start=Top avoid");
    imp_tmp = IJ.getImage()
    IJ.run(imp_tmp, "Z Project...", "projection=[Max Intensity]");
    imp_projected = IJ.getImage()
  elif projection=='Y':
    IJ.run(imp, "Reslice [/]...", "output=1.000 start=Left avoid");
    imp_tmp = IJ.getImage()
    IJ.run(imp_tmp, "Z Project...", "projection=[Max Intensity]");
    imp_projected = IJ.getImage()

  return imp_projected

def project_XYZ(imp):
  IJ.run(imp, "XYZ MaxProject", "");
  imp_projected = IJ.getImage()
  return imp_projected

#
#
# Projection
#

def analyze(iDataSet, tbModel, p):
  
  #
  # Init
  #
  IJ.run("Options...", "iterations=1 count=1"); 
  close_all_image_windows()
    
  #
  # Open file
  #
  filepath = tbModel.getFileAbsolutePathString(iDataSet, "Input", "IMG")
  
  if filepath.endswith('.mha'):
    IJ.openImage(filepath)
    imp = IJ.getImage()
    # get rid of the unsigned int problem
    IJ.run(imp, "32-bit", "");
    IJ.setMinAndMax(0, 65535);
    IJ.run(imp, "16-bit", "");

  else:
    imp = IJ.openImage(filepath)

  #
  # Crop after loading
  #  
  if p['crop']:
    x_min = int(p['mask_roi'][0])
    y_min = int(p['mask_roi'][1])
    z_min = int(p['mask_roi'][2])
    x_width = int(p['mask_roi'][3])
    y_width = int(p['mask_roi'][4])
    z_width = int(p['mask_roi'][5])
  
    stack = imp.getImageStack()
    stack_cropped = stack.crop(x_min, y_min, z_min, x_width, y_width, z_width)
    imp = ImagePlus('',stack)
  
  #
  # Convert
  #
  IJ.setMinAndMax(p['map_to_zero'], p['map_to_max']);
  IJ.run(imp, p['bit_depth'], "");

  #
  # Scale (Binning)
  #
  if p['binning']:
    sx = float(1.0/p['binning_x'])
    sy = float(1.0/p['binning_y'])
    sz = float(1.0/p['binning_z'])
    parameters =  "x="+str(sx)+" y="+str(sy)+" z="+str(sz)+" interpolation=Bilinear average process create"
    IJ.run(imp, "Scale...", parameters);
    imp = IJ.getImage()


  #
  # Save converted volume data
  #
  if p['save_volume_data']:
    converted_filename = tbModel.getFileName(iDataSet, "Input", "IMG")+"--converted."+p['output_format']
    tbModel.setFileAbsolutePath(p['output_folder'], converted_filename, iDataSet, 'Output', 'IMG')
    IJ.saveAs(imp, p['output_format'], tbModel.getFileAbsolutePathString(iDataSet, 'Output', 'IMG'))
    
  
  #
  # Projections 
  #
  
  if p['save_xyz_projections']:
    # compute one XYZ projection in one image
    projection = 'XYZ'
    imp_projected = project_XYZ(imp)
    # construct output filename 
    projected_filename = tbModel.getFileName(iDataSet, "Input", "IMG")+"--"+projection+"."+p['output_format']
    # store in table
    tbModel.setFileAbsolutePath(p['output_folder'], projected_filename, iDataSet, projection, "IMG")
    # store on disk
    IJ.saveAs(imp_projected,p['output_format'],tbModel.getFileAbsolutePathString(iDataSet, projection, "IMG"))
    '''
    projections = ['X','Y','Z']
  
    for projection in projections:
      # compute projection  
      imp_projected = project(imp, projection)
      # construct output filename 
      projected_filename = tbModel.getFileName(iDataSet, "Input", "IMG")+"--"+projection+"."+p['output_format']
      # store in table
      tbModel.setFileAbsolutePath(p['output_folder'], projected_filename, iDataSet, projection, "IMG")
      # store on disk
      IJ.saveAs(imp_projected,p['output_format'],tbModel.getFileAbsolutePathString(iDataSet, projection, "IMG"))
    '''

   
         
  return tbModel


'''
def convert_mha_to_tif(path_mha, path_tif):
  print("    reading mha")
  start_time = time.time()
  IJ.openImage(path_mha)
  print("      time elapsed: "+str(round(time.time()-start_time,3)))
  imp = IJ.getImage()
  IJ.run(imp, "32-bit", "");
  IJ.setMinAndMax(0, 65535);
  IJ.run(imp, "16-bit", "");
  output_file = os.path.join(output_folder,original_filename+"--transformed.tif")
  print("    writing tif")
  start_time = time.time()
  IJ.saveAs(imp, "Tiff", output_file)
  print("      time elapsed: "+str(round(time.time()-start_time,3)))
  # optional to save maximum projections
  if(save_as_max):
    if(imp.getStackSize()>1):
      print("    make and save z-max")
      start_time = time.time()
      IJ.run(imp, "Z Project...", "projection=[Max Intensity]")
      IJ.saveAs(IJ.getImage(), "Tiff", output_file+"--z-max.tif")
      print("      time elapsed: "+str(round(time.time()-start_time,3)))
  return output_file
'''

 
#
# Determine input files
#

def get_file_list(foldername, reg_exp):

  print("#\n# Finding files in: "+foldername+"\n#")
  pattern = re.compile(reg_exp)
   
  files = []
  for root, directories, filenames in os.walk(foldername):
	for filename in filenames:
	   print("Checking:", filename)
	   if filename == "Thumbs.db":
	     continue
	   match = re.search(pattern, filename)
	   if (match == None) or (match.group(0) == None):
	     continue
	   files.append(os.path.join(foldername, filename))  
	   print("Accepted:", filename)	   

  return(sorted(files))

#
# GET PARAMETERS
#

def get_parameters(p):
  gd = GenericDialogPlus("Please enter parameters")

  for k in p['expose_to_gui']['value']:
    if p[k]['type'] == 'boolean':
      gd.addCheckbox(k, p[k]['value'])
    elif p[k]['type'] == 'folder':
      gd.addDirectoryField(k, p[k]['value'], 100)	
    elif p[k]['type'] == 'file':
      gd.addFileField(k, p[k]['value'], 100)	
    elif p[k]['type'] == 'string':
      if p[k]['choices']:
        gd.addChoice(k, p[k]['choices'], p[k]['value'])	
      else:
        gd.addStringField(k, p[k]['value'])	 
    elif p[k]['type'] == 'int':
      if p[k]['choices']:
        gd.addChoice(k, p[k]['choices'], p[k]['value'])	
      else:
        gd.addNumericField(k, p[k]['value'], 0)	 
    elif p[k]['type'] == 'float':
      gd.addNumericField(k, p[k]['value'], 2)
  
  gd.showDialog()
  if gd.wasCanceled():
    return

  for k in p['expose_to_gui']['value']:
    if p[k]['type'] == 'boolean':
      p[k]['value'] = gd.getNextBoolean()
    elif p[k]['type'] == 'folder' or p[k]['type'] == 'file':
      p[k]['value'] = gd.getNextString()
    elif p[k]['type'] == 'string':
      if p[k]['choices']:
        p[k]['value'] = gd.getNextChoice()	
      else:
        p[k]['value'] = gd.getNextString()	 
    elif p[k]['type'] == 'int':
      if p[k]['choices']:
        p[k]['value'] = int(gd.getNextChoice())	
      else:
        p[k]['value'] = int(gd.getNextNumber()) 
    elif p[k]['type'] == 'float':
        p[k]['value'] = gd.getNextNumber()
    
  return p

    
if __name__ == '__main__':

  print("#\n# Elastix registration\n#")
    
  #
  # GET PARAMETERS
  #
  print("#\n# Parameters\n#")

  #
  # Load gui parameters
  #

  od = OpenDialog("Select parameter file (press CANCEL if you don't have one)", None)
  f = od.getPath()
  
  if f:
    print('loading parameters from file')
    f = open(f, 'r'); p_gui = pickle.load(f); f.close()
  else:
    print('starting from default parameters')
    # make parameter structure if it has not been loaded
    p_gui = {}
    # exposed to GUI
    p_gui['expose_to_gui'] = {'value': ['input_folder', 'reg_exp', 'output_folder', 'output_format', 
    'bit_depth', 'map_to_zero', 'map_to_max', 'binning',
    'binning_x','binning_y','binning_z','save_xyz_projections','save_binned_volume_data']}
    p_gui['input_folder'] = {'choices': '', 'value': 'C:\\Users\\acquifer\\Desktop\\882-reg3', 'type': 'folder'}
    p_gui['reg_exp'] = {'choices': '', 'value': '.*--transformed.mha$', 'type': 'folder'}
    p_gui['output_folder'] = {'choices': '', 'value': 'C:\\Users\\acquifer\\Desktop\\xyz01', 'type': 'folder'}
    p_gui['output_format'] = {'choices': ['tif'], 'value': 'tif', 'type': 'string'}
    p_gui['bit_depth'] = {'choices': ['8-bit','16-bit'], 'value': '16', 'type': 'string'}
    p_gui['map_to_zero'] = {'choices':'', 'value': 0, 'type': 'int'}
    p_gui['map_to_max'] = {'choices':'', 'value': 65535, 'type': 'int'}
    p_gui['save_binned_volume_data'] = {'choices': '', 'value': True, 'type': 'boolean'}

    p_gui['binning'] = {'choices': '', 'value': False, 'type': 'boolean'}
    p_gui['binning_x'] = {'choices': '', 'value': 1, 'type': 'float'}
    p_gui['binning_y'] = {'choices': '', 'value': 1, 'type': 'float'}
    p_gui['binning_z'] = {'choices': '', 'value': 1, 'type': 'float'}
    
    p_gui['save_xyz_projections'] = {'choices': '', 'value': True, 'type': 'boolean'}
    p_gui['save_volume_data'] = {'choices': '', 'value': True, 'type': 'boolean'}

    
  #
  # Expose parameters to users
  #
  print p_gui
  p_gui = get_parameters(p_gui)
  print p_gui
  
  #
  # Create derived paramters
  #
  
  #
  # Save gui parameters
  #
  f = open(os.path.join(p_gui['output_folder']['value'], 'fiji-convert-bin-project-gui-parameters.txt'), 'w')
  pickle.dump(p_gui, f)
  f.close()
   
  #
  # Reformat gui parameters for actual usage
  # 
  p = {}
  for k in p_gui.keys():
    p[k] = p_gui[k]['value']
    
  p['crop'] = False; # todo: implement this
  
  #
  # DETERMINE INPUT FILES
  #
      
  tbModel = TableModel(p['input_folder'])
  files = get_file_list(p['input_folder'], p['reg_exp'])
  

  #
  # INIT INTERACTIVE TABLE
  #
  
  tbModel.addFileColumns('Input','IMG')
  tbModel.addFileColumns('Output','IMG')
  tbModel.addFileColumns('XYZ','IMG')
  
  '''
  tbModel.addFileColumns('X','IMG')
  tbModel.addFileColumns('Y','IMG')
  tbModel.addFileColumns('Z','IMG')
  '''
  
  #for ch in p["channels"]:
  #  tbModel.addFileColumns('Input_'+ch,'IMG')
  
  #for ch in p["channels"]:
  #  tbModel.addFileColumns('Transformed_'+ch,'IMG')
  
  sorted_files = sorted(files)
  print("#\n# Files to be analyzed\n#")
  #for ch in p["channels"]:
  iDataSet = 0
  for afile in sorted_files:  #if ch in afile.split(os.sep)[-1]: if ch == p["channels"][0]:
    tbModel.addRow()
    print(str(iDataSet)+": "+afile)
    tbModel.setFileAbsolutePath(afile, iDataSet, "Input","IMG")
    iDataSet = iDataSet + 1

  #frame=ManualControlFrame(tbModel)
  #frame.setVisible(True)
  
  #
  # ANALYZE
  #
  print("#\n# Analysis\n#")

  n_files = tbModel.getRowCount()
  
  for iDataSet in range(n_files):
    #print(p['reference_image_index'],i)
    # compute transformation and transform reference channel
    tbModel = analyze(iDataSet, tbModel, p)
    close_all_image_windows()
    
  #
  # clean up
  #
  close_all_image_windows()
  
  #
  # show the interactive table
  #
  frame = ManualControlFrame(tbModel)
  frame.setVisible(True)

  print("done!")

  
  