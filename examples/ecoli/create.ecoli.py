# Pizza.py script to generate E Coli geometry

c = cdata()
c.seed(18485)

# cylindrical E Coli cell
# q params = triangle resolution of cell membrane

c.cap("cap",'x',0.5,0.5,0.5,0.5,2.3)
c.q("cap",10,15)
c.surf("cell","cap")

# 1260 = 36x35 planar array of CheA at anterior end

c.partarray("CheA",1,36,35,0.2,0.2375,0.245,0.1,0.015,0.015)
c.project("CheA","cell",0.5,0.5,0.5,1.0e-6,1)

# 8200 CheY randomly in cytoplasm

c.part("CheY",8200,"cell")

# 1600 CheZ randomly in cytoplasm

c.part("CheZ",1600,"cell")

# 34 FliM in each of 4 rings equally spaced along cell axis

c.partring("FliM1",34,0.86,0.0,0.5,0.0225,'y')
c.project("FliM1","cell",0,1,0,1.0e-6)
c.partring("FliM2",34,1.22,0.4999,1.0,0.0225,'z')
c.project("FliM2","cell",0,0,1,1.0e-6)
c.partring("FliM3",34,1.58,1.0,0.5,0.0225,'y')
c.project("FliM3","cell",0,1,0,1.0e-6)
c.partring("FliM4",34,1.94,0.4999,0.0,0.0225,'z')
c.project("FliM4","cell",0,0,1,1.0e-6)

# dump to data file

c.unselect("cap")
c.write("data.ecoli")
