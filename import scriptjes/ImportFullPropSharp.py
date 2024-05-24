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
        
        rail1 = adsk.core.ObjectCollection.create()
        rail2 = adsk.core.ObjectCollection.create()
        for x in allFiles:  

            with io.open(x, 'r', encoding='utf-8-sig') as f:
                ind = 0
                points1 = adsk.core.ObjectCollection.create()
                points2 = adsk.core.ObjectCollection.create()
                line = f.readline()
                data = []
                while line:
                    ind = ind + 1
                    pntStrArr = line.split(',')
                    for pntStr in pntStrArr:
                        try:
                            data.append(float(pntStr))
                        except:
                            break
                
                    if len(data) >= 3:                        
                        point = adsk.core.Point3D.create(data[0], data[1], data[2])
                        if ind <= num_lines/2+1:
                            points1.add(point)
                        if ind >= num_lines/2+1:
                            points2.add(point)
                        if ind == 1:
                            point = adsk.core.Point3D.create(data[0], data[1], data[2])
                            rail1.add(point)
                        if ind == num_lines/2+1:
                            point = adsk.core.Point3D.create(data[0], data[1], data[2])
                            rail2.add(point)
                    line = f.readline()
                    data.clear()            
            if points2.count:
                root = design.rootComponent
                sketch = root.sketches.add(root.xYConstructionPlane)    
                constraints = sketch.geometricConstraints            
                spline1 = sketch.sketchCurves.sketchFittedSplines.add(points1)
                spline2 = sketch.sketchCurves.sketchFittedSplines.add(points2)                   
                # startPoint2 = spline2.fitPoints.item(0)
                # startCurve2 = spline2.getCurvatureHandle(startPoint2)
                constraints.addTangent(spline2, spline1)
                # tan1 = spline1.getTangentHandle()
                # tan2 = spline2.getTangentHandle(spline2.fitPoints.item(0))


                #sketch.geometricConstraints.addTangent(spline1.getCurvatureHandle(,)
                # fitpoint1 = spline.fitPoints.item(0)
                # fitpoint2 = spline.fitPoints.item(1)
                # fitpoint3 = spline.fitPoints.item(num_lines-1)
                # sketch.sketchCurves.sketchLines.addByTwoPoints(fitpoint1,fitpoint2)
                # sketch.sketchCurves.sketchLines.addByTwoPoints(fitpoint1,fitpoint3)
            else:
                ui.messageBox('No valid points', title)    

        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.sketchCurves.sketchFittedSplines.add(rail1)
        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.sketchCurves.sketchFittedSplines.add(rail2)
        
            
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


