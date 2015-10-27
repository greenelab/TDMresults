from sys import float_info
import numpy

# Simple function to get min and max values from a list.
def getMinMax(values):
    maxValue = 0.0
    minValue = float_info.max
    for value in values:
        floatVal = float(value)
        if floatVal > maxValue:
            maxValue = floatVal
        if floatVal < minValue:
            minValue = floatVal
    return (minValue, maxValue)

def getPercentiles(values):
    values = [float(val) for val in values]
    return numpy.percentile(values, [25,50,75])
