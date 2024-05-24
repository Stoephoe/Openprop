#Author: Aran Dokoupil
#Adapted from Autodesk Description-Import spline from csv file

import adsk.core, adsk.fusion, traceback
import io
import os

class UiLogger:
    def __init__(self, forceUpdate):  
        app = adsk.core.Application.get()
        ui  = app.userInterface
        palettes = ui.palettes
        self.textPalette = palettes.itemById("TextCommands")
        self.forceUpdate = forceUpdate
        self.textPalette.isVisible = True 
    
    def print(self, text):       
        self.textPalette.writeText(text)
        if (self.forceUpdate):
            adsk.doEvents() 

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui  = app.userInterface
        #logger = UiLogger(True)
        # Get all components in the active design.
        product = app.activeProduct
        design = adsk.fusion.Design.cast(product)
        title = 'Import Spline csv'
        if not design:
            ui.messageBox('No active Fusion design', title)
            return
        
        dlg = ui.createFileDialog()
        dlg.title = 'Open First CSV File'
        dlg.filter = 'Comma Separated Values (*.csv);;All Files (*.*)'
        if dlg.showOpen() != adsk.core.DialogResults.DialogOK :
            return
        



        filename = dlg.filename
        num_lines = -1

        with io.open(filename, 'r', encoding='utf-8-sig') as f:            
            line = f.readline()            
            while line:
                num_lines = num_lines + 1   
                line = f.readline()           

        #logger.print(str(num_lines))
        allFiles = []


        
        for root, dirs, files in os.walk(os.path.dirname(filename)):
            for file in files:
                if file.endswith(".csv"):
                        fileItt = os.path.join(root, file)
                        fileItt = fileItt.replace("\\","/") 
                        allFiles.append(fileItt)                      

        for x in allFiles:  

            with io.open(x, 'r', encoding='utf-8-sig') as f:
                points = adsk.core.ObjectCollection.create()
                line = f.readline()
                data = []
                while line:
                    pntStrArr = line.split(',')
                    for pntStr in pntStrArr:
                        try:
                            data.append(float(pntStr))
                        except:
                            break
                
                    if len(data) >= 3 :
                        point = adsk.core.Point3D.create(data[0], data[1], data[2])
                        points.add(point)
                    line = f.readline()
                    data.clear()            
            if points.count:
                root = design.rootComponent
                sketch = root.sketches.add(root.xYConstructionPlane)
                sketch.sketchCurves.sketchFittedSplines.add(points)
            else:
                ui.messageBox('No valid points', title)            
            
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


