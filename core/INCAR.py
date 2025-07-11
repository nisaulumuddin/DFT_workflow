from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

class INCAR:
        def __init__(self,filepath):

                with open(filepath,'r') as f:
                        self.file_data = [i.split() for i in f.read().splitlines()]
                
                self.tags = {}
                for fields in self.file_data:
                        if len(fields) != 0:
                                try:
                                        if len(" ".join(fields).split("!")[0]) != 0:
                                                line = " ".join(fields).split("!")[0]
                                                tagname = line.split("=")[0].split()[0]
                                                tagdat = line.split("=")[1].split()
                                                self.tags[tagname] = tagdat

                                except:
                                        continue

                                #if len(" ".join(fields).split("!")[0]) != 0:
                                #        line = " ".join(fields).split("!")[0]
                                #        tagname = line.split("=")[0].split()[0]
                                #        tagdat = line.split("=")[1].split()
                                #        self.tags[tagname] = tagdat

                

        def delete_tag(self, tag):
                if tag in self.tags:
                        del self.tags[tag]
                        print(f"Tag '{tag}' has been deleted.")
                else:
                        print(f"Tag '{tag}' does not exist.")
                        
        def set_tag(self,tag,value):
                self.tags[tag] = [value]         

        def get_sorted_tags_string(self):
                taglist = '\n'.join([str(tag) + ' = ' + ' '.join([str(i) for i in self.tags[tag]]) for tag in self.tags])
                return taglist

        def write_out(self,filepath):
                if os.path.isfile(filepath):
                        os.remove(filepath)

                with open(filepath,'a') as f:
                        f.write(self.get_sorted_tags_string())
        
                
