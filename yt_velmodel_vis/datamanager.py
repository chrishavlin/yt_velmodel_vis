import os
class filesysDB(object):
    def __init__(self,top_level_dir=None,verbosity=0):
        '''
        simple file system database manager

        Parameters
        ----------
        top_level_dir   the top level directory to search for data files. If
                        None, will check for environment variable YTVELMODELDIR
        '''
        if top_level_dir is None:
            top_level_dir = os.environ.get('YTVELMODELDIR')

        if top_level_dir is None:
            raise(("Could not intialize filesysDB: top_level_dir is None and "
                   "environment variable YTVELMODELDIR is not set."))

        self.db_path=top_level_dir
        self.buildFileDictionary()
        return

    def buildFileDictionary(self):
        self.FileDict = {}
        for root, subdirs, files in os.walk(self.db_path):
            fileList=files # files in this root
            for file in files:
                self.FileDict[file]=os.path.join(root,file)
        return

    def validateFile(self,fname):
        validFile=False
        if os.path.isfile(fname):
            validFile=fname
        else:
            if fname in self.FileDict.keys():
                validFile=self.FileDict[fname]

        return validFile
