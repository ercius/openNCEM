from io import StringIO
import numpy as np
import pickle
import itertools
import os, sys, shutil

ElamLoaded = False # We only want to load the Elam database once.  Keep track of if that has been done.  False = not yet loaded.
ElamData = []

# For debugging we sometimes want Elam to print out the lines it loaded.
ElamPrint = lambda s: None #print(s)

# Reading the Elam database is time consuming.  This caches the fluorescence data so it is only read once per program
#  load.
ElamFluorescenceLines = {'nodata':1}

def GetElamFluorescenceLines(ElementName, E=None, I=None):
    global ElamLoaded, ElamData, ElamFluorescenceLines

    # Already loaded fluorescence data for this element?
    if ElementName in ElamFluorescenceLines:
        # Return the (previously or newly) cached element info.
        return ElamFluorescenceLines[ElementName]

    # This is really time consuming, so we also cache it as a pickle file.
    try:
        ElementXRayLines = pickle.load(open(os.path.join(os.path.dirname(__file__), 'Elam', ElementName.lower() +'.pickle'), 'rb'))
        # It was pickled!  So cache it in memory and return it.
        ElamFluorescenceLines[ElementName] = ElementXRayLines
        return ElementXRayLines
    except:
        # It wasn't pickled, no matter.  We'll make it now and pickle it at the end of the function.
        pass

    # It's not already loaded, so continue with loading it.

    # Read the entire Elam database into memory.
    if ElamLoaded == False:
        fo = open(os.path.join(os.path.dirname(__file__), 'Elam', 'ElamDB12.txt'), "r+")
        str = fo.read();
        fo.close()
        ElamData = StringIO(str)
        ElamLoaded = True

    # Initialize the lists that will contain the line information.
    LineNamesIUPAC = []
    LineNamesSiegbahn = []
    Geoms = []

    # Make sure we're at the start of the file
    ElamData.seek(0)

    # Zoom forward to the element we are looking for.
    for line in ElamData:
        if ('Element '+ElementName+' ') in line:
            break

    # Pull out the section with the fluorescence lines, which is first.
    for line in ElamData:
        # Ignore the lines that aren't the fluorescence info
        if ('Edge' in line) or ('Lines' in line) or ('CK' in line):
            continue
        # When we get to the Photoionization section we are done.
        if 'Photo' in line:
            break;
        # Otherwise, this is a fluorescence line.  Add it to our output if it is intense enough to improve the fit
        q = line.split()
        if (float(q[3]) > 0.01):
            LineNamesIUPAC.append(q[0])
            LineNamesSiegbahn.append(q[1])
            Geoms.append([float(q[2]), float(q[3])])
            ElamPrint('Added: ' + line)
        else:
            ElamPrint('Ignored (intensity < 1%): ' + line)

    Geoms = np.array(Geoms)

    # We are going to use the Siegbahn names, not IUPAC.
    # Make a dictionary where the key is the line name and the value is an (energy, intensity) tuple.
    Lines = dict(zip(LineNamesSiegbahn, Geoms))

    ElementXRayLines = {}

    for Series in ['K', 'L', 'M', 'N', 'O']:

        #Obtain the name strings for only the lines that have names that start with K (or L, or M, whatever Series we are on).
        ThisSeries = dict([(k, v) for k, v in Lines.items() if k.startswith(Series)])

        # If these series is not present, then don't keep it.
        if len(ThisSeries) == 0:
            continue

        # Make this edge in ElementXRayLines
        ElementXRayLines[Series] = ThisSeries

    # We've built the fluorescence data for this element.  Cache it and pickle it for next time.
    ElamFluorescenceLines[ElementName] = ElementXRayLines
    pickle.dump(ElementXRayLines, open(os.path.join(os.path.dirname(__file__), 'Elam', ElementName.lower() + '.pickle'), 'wb'))

    # And return what we've got.
    return ElementXRayLines

def GetWeightedSum(Line, Series):
    # The key is that the series we want to sum has to have a common prefix.
    # For example, if Line = 'Ka' then we can sum 'Ka1', 'Ka2' but not 'Kbx'
    # Get only the lines that match our specific subseries.
    vals = np.array([Vals for LineName, Vals in Series.items() if Line in LineName])
    # Line weights are summed to one for the entire series.  Since we often want to weigh a subseries, renormalize.
    vals[:,1] /= np.sum(vals[:,1])
    # Now get a weighted sum: sum(position*relative_intensity)
    WeightedEnergy = np.sum(vals[:,0]*vals[:,1])
    return WeightedEnergy

def GetFluorescenceLineEnergy(ElementName, Series='K', Line=None):
    LineData = GetElamFluorescenceLines(ElementName)

    # If the user requests a "line" which is the series, we will return the weighted sum of all lines in the series (below in the if).
    if Line == Series:
        Line = None

    # If the user specifies a specific line, then we will return its energy.
    try:
        if Line in ['Ka', 'Kb', 'La', 'Lb']:
            # We can make weighted sums of special subseries of lines. 
            return GetWeightedSum(Line, LineData[Series])
        elif Line is not None:
            # Or we can return a specific line which is requested.
            # print(LineData)
            return LineData[Series][Line][0]
        else:
            # The user doesn't specify a line which means he wants to know the energy of the whole series.
            return GetWeightedSum(Series, LineData[Series])
    except:
        # If the user asked for a series or line that's not in the database, then we return none.
        return None

if __name__ == '__main__':
    print('------------------------------')
    # print('Fe-K: %s' % GetFluorescenceLineEnergy('Fe'))
    # print('Al-K: %s' % GetFluorescenceLineEnergy('Al'))
    # print('Mg-K: %s' % GetFluorescenceLineEnergy('Mg'))
    # print('Mg-K: %s' % GetFluorescenceLineEnergy('Mg', Series='K'))
    # print('Mg-K: %s' % GetFluorescenceLineEnergy('Mg', Series='K', Line='K'))
    # print('Fe-Q: %s' % GetFluorescenceLineEnergy('Fe', Series='Q'))
    # print('Fe-Ka2: %f' % GetFluorescenceLineEnergy('Fe', Series='K', Line='Ka2'))
    # print('Fe-Ka1: %f' % GetFluorescenceLineEnergy('Fe', Series='K', Line='Ka1'))
    # print('Fe-Ka: %f' % GetFluorescenceLineEnergy('Fe', Series='K', Line='Ka'))
    # print('Fe-Kb: %f' % GetFluorescenceLineEnergy('Fe', Series='K', Line='Kb'))
    # print('Mn-Ka: %f' % GetFluorescenceLineEnergy('Mn', Series='K', Line='Ka'))
    # print('Pt-La: %f' % GetFluorescenceLineEnergy('Pt', Series='L', Line='La'))
    print('S-K: %f' % GetFluorescenceLineEnergy('S', Series='K', Line='K'))
