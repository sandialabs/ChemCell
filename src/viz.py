# viz tmp.dump as ChemCell model

c = cdata("data.er.surf")
d = dump("tmp.dump")
#d.scale()
d.extra(c)

g = gl(d)
g.arad(0,0.02)
g.arad(8,0.074)
#g.arad([1,3,5,7],0.3)
v = vcr(g)
