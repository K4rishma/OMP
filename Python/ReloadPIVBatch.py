import os
from tqdm import tqdm, trange
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory  # python 3
from PathlinePlot import saveVectorStreams


initialdir = 'K:/Results/Phrase/PIV'
root = Tk()
root.withdraw()  # we don't want a full GUI, so keep the root window from appearing
root.wm_attributes('-topmost', 1)  # force on top
# filename = askopenfilename(parent=root, initialdir=dataPath,
#                           title='Select a file')  # show an "Open" dialog box and return the path to the selected file
rootdir = askdirectory(parent=root, initialdir=initialdir,
                       title='Select the lowest level root of all files to be loaded')
subdirs = [name for name in os.listdir(rootdir) if os.path.isdir(os.path.join(rootdir, name))]
subdirs = [name for name in subdirs if name.startswith('2018') or name.startswith('2019')]
# choose which settings file to open
filename = askopenfilename(parent=root, initialdir=os.path.join(rootdir, subdirs[0]),
                           filetypes=(("Matlab files", "*.mat"), ("all files", "*.*")),
                           title='Select which settings file to open')
_, filename = os.path.split(filename)
settings = filename[filename.find('_')+1:]

##
for subdir in tqdm(subdirs, desc='file loop'):
    dirname = os.path.join(rootdir, subdir)
    fname, = [name for name in os.listdir(dirname) if name.endswith(settings)]
    datename = fname[:fname.find('_')]
    fullpath = os.path.join(dirname, fname)
    saveVectorStreams(fullpath, datename)