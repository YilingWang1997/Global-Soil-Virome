import gdal
import numpy as np 
import sys
import pandas as pd

def CS_Transport(label_path):
    reader= pd.read_csv(label_path, header=0, encoding='utf-8')
    data=reader.values.tolist()

    return data


def GetValue():
    file1 = sys.argv[1]
    label_path = r'$input.csv' #Location information included latitude and longitude of all the sample sites
    data = CS_Transport(label_path)
    file1 = gdal.Open(file1) 
    band1 = file1.GetRasterBand(1)
    print(len(data))
    for i in range(len(data)):
        position_X = float(data[i][3])
        position_Y = float(data[i][2])
        p_X = int((position_X+180)/0.01)
        p_Y = int((90-position_Y)/0.01)
        Array1 = band1.ReadAsArray(xoff = p_X,yoff = p_Y,win_xsize = 1,win_ysize = 1)
        print (Array1[0][0])
        

GetValue()
#GetValue()