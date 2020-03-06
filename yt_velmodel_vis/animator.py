import numpy as np
import yt
import os
import imageio as iio

class animator(object):
    def __init__(self,input_dir,ftype='png'):
        self.input_dir=input_dir
        self.ftype=ftype
        self.prepFileList()
        return

    def prepFileList(self):
        '''
        finds all files of type self.ftype, builds sorted list of full paths
        '''
        self.fileList=[]
        for fi in os.listdir(path=self.input_dir):
            if fi.split('.')[-1].lower()==self.ftype:
                self.fileList.append(os.path.join(self.input_dir,fi))

        self.fileList.sort()
        self.N_frames=len(self.fileList)
        print(str(self.N_frames)+' frames found in '+self.input_dir)

    def writeGif(self,outfile,sec_per_frame=1):
        '''
        writes a gif using imageio

        parameters
        ----------
        outfile           the file to write
        sec_per_frame     seconds per frame (inverse of frames per sec)
        '''
        with iio.get_writer(outfile, mode='I',duration=sec_per_frame) as writer:
            for filename in self.fileList:
                writer.append_data(iio.imread(filename))
