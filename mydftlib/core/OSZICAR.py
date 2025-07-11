from __future__ import annotations

class OSZICAR:
        def __init__(self,filepath):
                self.mag = []

                with open(filepath,'r') as f:
                        self.file_data = f.read().splitlines()

                        for num,line in enumerate(self.file_data):
                                if 'mag=' in line:
                                        self.mag.append(float(line.split()[9]))