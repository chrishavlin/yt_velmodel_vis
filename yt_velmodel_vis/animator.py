import numpy as np
import yt
import os
import imageio as iio
from PIL import Image
import shutil

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


class TFsequence(animator):
    def __init__(self,volBase,TFbase,input_dir,ftype='png',tempDir=None):
        self.input_dir=input_dir
        self.volBase=volBase # volume rendering file base
        self.TFbase=TFbase # transfer function file base
        self.ftype=ftype
        self.prepFileList()

        if tempDir is None:
            self.tempDir=os.path.join(self.input_dir,'animator_temp')

        return

    def prepFileList(self):

        self.fileDict={}
        for fi in os.listdir(path=self.input_dir):
            isTF=self.TFbase in fi
            isVol=self.volBase in fi
            if fi.split('.')[-1].lower()==self.ftype and (isTF or isVol):
                frame=''.join([d for d in filter(str.isdigit,fi)])
                if frame not in self.fileDict.keys():
                    self.fileDict[frame]={}
                if isTF:
                    self.fileDict[frame]['tf']=os.path.join(self.input_dir,fi)
                else:
                    self.fileDict[frame]['vol']=os.path.join(self.input_dir,fi)

        self.frameList=list(self.fileDict.keys())
        self.frameList.sort()
        return

    def buildFullFrames(self):

        # prep the temporary directory for full frames
        if os.path.isdir(self.tempDir) is False:
            try:
                os.mkdir(self.tempDir)
            except:
                raise OSError('Could not write temporary directory '+self.tempDir)


        # loop over files, build the frames
        self.fileList=[]
        for frame in self.frameList:
            try:
                volFi=self.fileDict[frame]['vol']
            except:
                raise keyError('volume render missing for frame '+frame)

            try:
                tfFi=self.fileDict[frame]['tf']
            except:
                raise keyError('transfer function missing for frame '+frame)

            vol=Image.open(volFi)
            TF=Image.open(tfFi)

            vol_w=vol.size[0]
            TF_w=TF.size[0]
            heights=[TF.size[1],vol.size[1]]
            total_width = vol_w+int(TF_w/2.)
            max_height = max(heights)

            new_im = Image.new('RGB', (total_width, max_height))

            new_im.paste(vol, (int(TF_w/3.),0))
            new_im.paste(TF, (0,max_height-TF.size[1]))

            newFi=os.path.join(self.tempDir,frame+'.png')
            new_im.save(newFi)
            self.fileList.append(newFi)

        self.fileList.sort()

        return
    
    def buildGif(self,outputFile,deleteFrames=True):
        self.buildFullFrames()
        self.writeGif(outputFile)
        if deleteFrames:
            shutil.rmtree(self.tempDir)
        return
