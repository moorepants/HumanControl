import zipfile
import os

figures = {'Fig06': 'benchmarkSteerClosed.eps',
           'Fig07': 'benchmarkSteerOpen.eps',
           'Fig08': 'benchmarkSteerHandling.eps',
           'Fig09': 'openBode.eps',
           'Fig10': 'handling.eps',
           'Fig11': 'eigenvalues.eps',
           'Fig12': 'paths.eps',
           'Fig13': 'phiDistance.eps',
           'Fig14': 'deltaDistance.eps',
           'Fig15': 'TdeltaDistance.eps',
           'Fig17': 'benchmarkRollClosed.eps',
           'Fig18': 'benchmarkRollOpen.eps',
           'Fig19': 'rollDistance.eps',
           'Fig20': 'phasePortraits.eps',
           'FigA2a': 'Bike01-browserIns_sub.tiff',
           'FigA2b': 'Bike02-browser_sub.tiff',
           'FigA2c': 'Bike03-pista_sub.tiff',
           'FigA2d': 'Bike04-fisher_sub.tiff',
           'FigA2e': 'Bike05-yellow_sub.tiff',
           'FigA2f': 'Bike06-yellowRev_sub.tiff',
           'FigA1': 'benchmarkBicycle.ps'}

plotDir = 'plots'

ieeeFile = zipfile.ZipFile(os.path.join(plotDir, 'forIEEE.zip'), 'w')

for k, v in figures.items():
    ieeeFile.write(os.path.join(plotDir, v), k + '-' + v)

ieeeFile.close()
