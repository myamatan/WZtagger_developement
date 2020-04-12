'''
DecorationName: SmoothZContainedMaxSig

MassCutLow: x > 2500 ? 73.293 : sqrt(pow((14121.7)/x-(76.1159),2)+pow((-0.0194686)*x+(68.8269),2))

MassCutHigh: x > 2500 ? 122.981 : sqrt(pow((272.319)/x+(103.876),2)+pow((0.0403225)*x+(-35.1453),2))

D2Cut: x > 2500 ? 1.53059 : (0.737153)+(0.000183033)*x+(-8.80806e-08)*pow(x,2)+(5.67267e-11)*pow(x,3)
'''

from scipy.optimize import curve_fit
import numpy as np

from ROOT import TF1, TCanvas, TH2D

dt = np.loadtxt('../output/pyvalue.txt', delimiter=',')

print dt

pt = dt[:,0]
d2 = dt[:,1]
mu = dt[:,2]
ml = dt[:,3]

def d2_fit(x,a,b,c,d):
    return  a + b*x + c*x**2 + d*x**3
param, cov = curve_fit(d2_fit, pt, d2)

d2_smooth = d2_fit(pt, param[0], param[1], param[2], param[3])
print 'Params : '
print param

print 'Smoothed  D2 cut : '
print d2_smooth

'''
cv = TCanvas()

f1 = TF1('f1','[0]+[1]*x+[2]*x*x+[3]*x*x*x', 200, 2500)
f1.SetParameter(0, param[0])
f1.SetParameter(1, param[1])
f1.SetParameter(2, param[2])
f1.SetParameter(3, param[3])
f1.SetMinimum(0)
f1.SetMaximum(5)
f1.Update()

cv.cd()
h.Draw()
f1.Draw('same')
cv.SaveAs('../output/test.png')
'''

