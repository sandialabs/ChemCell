# histogram of cell density along E Coli x axis

d = dump("dump.ecoli.flim.long.2")
d.tselect.test("$t >= 10000")

d.aselect.test("$type == 5")
h = histo(d)
x,y1 = h.compute('x',50)
print x,y1

d.aselect.all()
d.aselect.test("$type == 4")
h = histo(d)
x,y2 = h.compute('x',50)
print x,y2

n = len(y2)
sum = n*[0]
frac = n*[0]

for i in xrange(len(y2)):
  sum[i] = y1[i] + y2[i]
  frac[i] = 1.0*y1[i] / sum[i]
  
g = matlab()
g.plot(x,frac)
